import json
import tqdm
import warnings
import pandas as pd
import numpy as np
import statsmodels.api as sm
from basic_tools import eps_selected, eps_dict
import seaborn as sns
from plot_tools import plot_odds_ratio
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeClassifier
warnings.filterwarnings('ignore')
scaler = StandardScaler()

SEED = 4
OUTCOME = 'T1D_STRICT'
eps = json.load(open('eps_' + OUTCOME + '.json', 'r'))
covariates = ['PC'+str(i) for i in range(1, 11)] + ['sex', 'age']
eps_sig = ['AUTOIMMUNE_HYPERTHYROIDISM',  'D3_ANAEMIA_B12_DEF', 'E4_HYTHY_AI_STRICT',
           'K11_COELIAC', 'M13_RHEUMA', 'T1D_STRICT']

# import imputed HLA data
hla_data = pd.read_csv('R9_imputed_HLAs_v2.csv')
genes = hla_data.ID.tolist()
hla_data = hla_data.transpose()
hla_data.columns = genes
hla_data = hla_data.iloc[9:, :]

# convert to dosage probability and create a sum table for genes
hla_df = hla_data[['A*01:01']]
gene_sum_df = pd.DataFrame(columns=[0.0, 1.0, 2.0])
for i in tqdm.tqdm(genes):
    genotype = hla_data[i].str.extract(':(.+),(.+),(.+)')
    genotype['dose_str'] = hla_data[i].str.extract('(.+):')
    genotype['dose_int'] = np.select([(genotype.dose_str == './.'), (genotype.dose_str == '0/0'),
                                      (genotype.dose_str == '0/1'), (genotype.dose_str == '1/1')],
                                     [-1.0, 0.0, 1.0, 2.0])
    a = genotype[0].astype(float) * 0 + genotype[1].astype(float) * 1 + genotype[2].astype(float) * 2
    a1, a2 = a[~a.isna()], a[a.isna()]
    a2 = genotype[genotype.index.isin(a2.keys())].dose_int
    hla_df[i] = a1.append(a2)
    gene_sum_df = gene_sum_df.append(round(hla_df[i]).value_counts().T)

gene_sum_df['genotype'] = gene_sum_df.index
gene_sum_df.index = range(len(gene_sum_df))
gene_sum_df = gene_sum_df.rename(columns={0.0: '0.0', 1.0: '1.0', 2.0: '2.0'})
# remove those genes without any alt allele: 187 -> 129
genes = gene_sum_df[(~gene_sum_df['1.0'].isna()) & (~gene_sum_df['2.0'].isna())].genotype.tolist()
hla_df = hla_df[genes]

# import confounding cols
confounding = pd.read_csv('confounding.csv')
confounding = confounding[confounding.PC1 != 0]
confounding = confounding[~confounding.age.isna()]
hla_df['finngen_id'] = hla_df.index
hla_df = hla_df.merge(confounding, 'inner', on='finngen_id')

# add endpoints cols and create a sum table for endpoints
events = pd.read_csv('/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_endpoint_longitudinal_1.0.txt.gz',
                     sep='\t')
for ep in tqdm.tqdm(eps):
    events_sub = events[events.ENDPOINT == ep]
    hla_df[ep] = np.select([(hla_df.finngen_id.isin(events_sub.FINNGENID)), (~hla_df.finngen_id.isin(events_sub.FINNGENID))],
                           [1, 0])

ep_sum_df = pd.DataFrame(hla_df[eps].sum(), columns=['n_cases'])
ep_sum_df['endpoint'] = ep_sum_df.index
ep_sum_df.index = range(len(ep_sum_df))
ep_sum_df['sex'] = [-1] * (len(eps) - 1) + [1]
ep_sum_df['n_cohort'] = [len(hla_df)] * (len(eps) - 1) + [len(hla_df[hla_df.sex == 1])]
ep_sum_df = ep_sum_df[['n_cohort', 'n_cases', 'endpoint', 'sex']]

# build a for-loop for logistic regression
columns = ['ep1', 'ep2', 'coef', 'se', 'pval', 'beta_025', 'beta_975', 'r_squared_0', 'r_squared_1', 'r_squared_delta',
           'n_ep1', 'n_ep2', 'n_both']


def regression(exposure, outcome, hla_dataset):
    y = hla_dataset[outcome]
    y.index = range(len(y))

    x0 = sm.add_constant(hla_dataset[covariates])
    model = sm.Logit(y, x0).fit(disp=0)
    r_squared_0 = float(model.summary2().tables[0].iloc[0, 3])

    x1 = sm.add_constant(hla_dataset[[exposure] + covariates])
    model = sm.Logit(y, x1).fit(disp=0)
    r_squared_1 = float(model.summary2().tables[0].iloc[0, 3])
    r_squared_delta = r_squared_1 - r_squared_0

    res = model.summary2().tables[1]
    beta_025 = res.loc[exposure, '[0.025']
    beta_975 = res.loc[exposure, '0.975]']
    pval = res.loc[exposure, 'P>|z|']
    se = res.loc[exposure, 'Std.Err.']
    coef = res.loc[exposure, 'Coef.']
    return [outcome, coef, se, pval, beta_025, beta_975, r_squared_0, r_squared_1, r_squared_delta]


def model_loop(endpoint, results, outcome_list, exposure, hla_dataset):
    for outcome in outcome_list:
        n_ep1 = len(hla_dataset[hla_dataset[endpoint] == 1])
        n_ep2 = len(hla_dataset[hla_dataset[outcome] == 1])
        n_both = len(hla_dataset[(hla_dataset[endpoint] == 1) & (hla_dataset[outcome] == 1)])
    res = regression(outcome, exposure, hla_dataset)
    if res:
        results = results.append(pd.Series([endpoint] + res + [n_ep1, n_ep2, n_both],
                                           index=results.columns), ignore_index=True)
    return results


# plot beta (95% CI)
def plot_association1(df, ylabel, color, hline, outcome=OUTCOME, plt_len=15):
    """
    :param df - a DataFrame of summary statistics
    :param outcome - a string which indicates the outcome disease name
    :return - a beta plot of all the diseases listed in df
    """
    df = df.sort_values(by='ep1')
    df.index = range(len(df))
    plt.figure(figsize=(plt_len, 5))
    plt.box(False)
    plt.grid()
    for i, row in df.iterrows():
        alpha = 1 if row.pval <= 0.05 / len(df) and row.r_squared_delta >= 0.01 else 0.12
        plt.plot((i, i), (row.beta_025, row.beta_975), 's', color=color, alpha=alpha)
        plt.plot(i, (row.beta_025 + row.beta_975) / 2, 's', color=color, alpha=alpha)
    plt.xticks(range(len(df)), [eps_dict[i] for i in df.ep1.tolist()], rotation=90)
    plt.ylabel(ylabel, size=12)
    plt.ylabel(y=hline, color='black', linestyle='--', linewidth=1)
    plt.grid()
    plt.show()


def plot_association2(df, sig_eps, ylabel, color, hline, outcome=OUTCOME, plt_len=15):
    """
    :param df - a DataFrame of summary statistics
    :param outcome - a string which indicates the outcome disease name
    :return - a beta plot of all the diseases listed in df
    """
    df = df.sort_values(by='ep1')
    df.index = range(len(df))
    plt.figure(figsize=(plt_len, 5))
    plt.box(False)
    plt.grid()
    for i, row in df.iterrows():
        alpha = 1 if row.pval <= 0.05 / len(df) else 0.12
        plt.plot((i, i), (row.beta_025, row.beta_975), 's', color=color, alpha=alpha)
        plt.plot(i, (row.beta_025 + row.beta_975) / 2, 's', color=color, alpha=alpha)
        if row.ep1 in sig_eps:
            plt.annotate('$', (i+.1, row.beta_975+.01), size=12, color='green')
    plt.xticks(range(len(df)), [eps_dict[i] for i in df.ep1.tolist()], rotation=90)
    plt.ylabel(ylabel, size=12)
    plt.ylabel(y=hline, color='black', linestyle='--', linewidth=1)
    plt.grid()
    plt.show()


# create a heatmap for the results
def plot_heatmap(stat_df):
    ep_gene_mat_pval = pd.pivot_table(stat_df, values='pval', index=['genotype'], columns=['endpoint'])
    ep_gene_mat_sig = np.select([(ep_gene_mat_pval > 0.05 / len(genes)), (ep_gene_mat_pval <= 0.05 / len(genes))],
                                ['', '*'])
    ep_gene_mat_sig = pd.DataFrame(ep_gene_mat_sig, columns=ep_gene_mat_pval.columns, index=ep_gene_mat_pval.index)
    ep_gene_mat_z = pd.pivot_table(stat_df, values='z', index=['genotype'], columns=['endpoint'])
    for i in ep_gene_mat_z[ep_gene_mat_z[OUTCOME] < 0].index:
        ep_gene_mat_z.loc[i, :] = ep_gene_mat_z.loc[i, :] * -1
    # convert endpoint codes to disease names
    ep_gene_mat_z.columns = [eps_dict[i] for i in ep_gene_mat_z.columns]
    # create a heatmap for the results
    sns.set(rc={'figure.figsize': (20, 6)})
    ax = sns.heatmap(ep_gene_mat_z, linewidths=.5, center=0, annot=ep_gene_mat_sig, cmap='RdBu', fmt='')
    ax.set_facecolor('#f7f7f7')
    ax.collections[0].colorbar.set_label('z score')


# shuffle and split
hla_df_shuffled = hla_df.sample(frac=1, random_state=SEED)
hla_df_shuffled.index = range(len(hla_df_shuffled))
# set for 10-fold CV
n_folds = np.array_split(hla_df_shuffled, 10)

method = '10fold_weighted_ridge'
res_prs = pd.DataFrame(columns=columns)
# obtain stats for selected endpoints and save the results
for endpoint in tqdm.tqdm(eps_selected):
    test_df_all = pd.DataFrame(columns=['finngen_id', 'prs_'+endpoint])
    ep_counts = hla_df[endpoint].value_counts()
    model_weights = {0: ep_counts[1]/ep_counts[0], 1: 1.}

    for fold in n_folds:
        test_df = fold
        test_df.index = range(len(test_df))
        train_df = hla_df_shuffled[~hla_df_shuffled.finngen_id.isin(test_df.finngen_id)]
        train_df.index = range(len(train_df))

        model = RidgeClassifier(class_weight=model_weights, random_state=SEED)
        model.fit(train_df[genes], train_df[endpoint])
        prs_weights = dict(zip(genes, model.coef_[0]))
        test_df['prs_'+endpoint] = 0
        test_df.index = range(len(test_df))
        for ep, weight in prs_weights.items():
            test_df['prs_'+endpoint] += test_df[ep]*weight
        test_df['prs_'+endpoint] = scaler.fit_transform(np.array(test_df['prs_'+endpoint]).reshape(-1, 1))
        test_df_all = pd.concat([test_df_all, test_df[['finngen_id', 'prs_'+endpoint]]], axis=0)

    hla_df = hla_df.merge(test_df_all, 'inner', on='finngen_id')
    res_prs = model_loop(endpoint, res_prs, [OUTCOME], 'prs_'+endpoint, hla_df)
    res_prs = model_loop(endpoint, res_prs, [endpoint], 'prs_' + endpoint, hla_df)

res_prs.to_csv('hla_prs_res.csv', index=None)

res_step1 = res_prs[res_prs.index.isin([i for i in range(len(eps_selected)*2) if i % 2 == 1])]
res_step1 = res_step1[res_step1.r_squared_delta >= 0.01]
res_step2 = res_prs[res_prs.index.isin([i for i in range(len(eps_selected)*2) if i % 2 == 0])]
res_step2 = res_step2[res_step2.ep1.isin(res_step1.ep1)]
# to make sure r_squared_delta in step1 are approximately equal or larger than the those in step2
print(res_step1.r_squared_delta - res_step2.r_squared_delta.values)


# plot the associations between each pair of endpoint and PRS
plot_association1(res_step1, 'beta (AD ~ HLA-PRS for AD)', '#285085', 0)
plot_association2(res_step2, eps_sig, 'beta (T1D ~ HLA-PRS for AD)', '#285085', 0, plt_len=7)

# plot number of cases
plt.figure(figsize=(7, 1))
plt.box(False)
plt.grid()
plt.bar(res_step2.ep1.tolist(), res_step2.n_ep1.tolist(), color='#285085', width=0.5)
plt.xticks(range(len(res_step2.ep1)), [eps_dict[i] for i in res_step2.ep1.tolist()], rotation=90)
plt.ylabel('# cases', size=12)
plt.yscale('log')
plt.yticks([100, 1000, 10000])
ax = plt.gca()
ax.get_yaxis().set_major_formatter(ScalarFormatter)
plt.grid()
plt.show()
