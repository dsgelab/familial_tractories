import json
import tqdm
import warnings
import pandas as pd
import numpy as np
import statsmodels.api as sm
import seaborn as sns
from plot_tools import plot_odds_ratio
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from sklearn.preprocessing import MinMaxScaler
warnings.filterwarnings('ignore')

OUTCOME = 'T1D_STRICT'
eps = json.load(open('eps_' + OUTCOME + '.json', 'r'))

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
    genotype['dose_str'] = hla_data[i].extract('(.+):')
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
genes = gene_sum_df[(~gene_sum_df[1.0].isna()) & (~gene_sum_df[2.0].isna())].genotype.tolist()
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
    hla_df[ep] = np.select([(hla_df.index.isin(events_sub.FINNGENID)), (~hla_df.index.isin(events_sub.FINNGENID))],
                           [1, 0])

ep_sum_df = pd.DataFrame(hla_df[eps].sum(), columns=['n_cases'])
ep_sum_df['endpoint'] = ep_sum_df.index
ep_sum_df.index = range(len(ep_sum_df))
ep_sum_df['sex'] = [-1] * (len(eps) - 1) + [1]
ep_sum_df['n_cohort'] = [len(hla_df)] * (len(eps) - 1) + [len(hla_df[hla_df.sex == 1])]
ep_sum_df = ep_sum_df[['n_cohort', 'n_cases', 'endpoint', 'sex']]

# build a for-loop for logistic regression
covariates = ['PC'+str(i) for i in range(1, 11)] + ['sex', 'age']
key_cols = ['Coef.', 'Std.Err.', 'z', 'P>|z|', '[0.025', '0.975]', 'endpoint', 'n_cases', 'n_cohort']


def modeling(subset, gene_list):
    # subset column order:
    # 0-n_cohort-int,
    # 1-n_cases-int,
    # 2-outcome-str,
    # 3-sex-int
    x = hla_df[gene_list + covariates]
    x = sm.add_constant(x)
    y = hla_df[[subset[2]]]
    if subset[3] > -1:
        x = x[x.sex == subset[3]]
        x = x.drop(columns=['sex'])
        y = y[y.index.isin(x.index)]
    try:
        model = sm.Logit(y.astype(int), x).fit(disp=0)
        stat = model.summary2().tables[1]
        stat = stat[stat.index.isin(gene_list)]
    except PerfectSeparationError:
        # When event per variable is too low,
        # the logistic model is under the risk of complete separation or quasi⁃complete separation.
        # https://m.medsci.cn/article/show_article.do?id=606319405019
        pd.DataFrame(dict(zip(key_cols[:6], [2.] * 6)), index=gene_list)
    except np.linalg.LinAlgError:
        # singular matrix
        pd.DataFrame(dict(zip(key_cols[:6], [3.] * 6)), index=gene_list)
    except ValueError:
        # ValueError: Pandas data cast to numpy dtype of object. Check input data with np.asarray(data).
        pd.DataFrame(dict(zip(key_cols[:6], [4.] * 6)), index=gene_list)
    except Exception:
        pd.DataFrame(dict(zip(key_cols[:6], [5.] * 6)), index=gene_list)
    stat['endpoint'] = subset[2]
    stat['n_cases'] = subset[1]
    stat['n_cohort'] = subset[0]
    return stat


def create_stat_df(ep_list, gene_list, single=True):
    original_results = pd.DataFrame(columns=key_cols)
    ep_sum = ep_sum_df[ep_sum_df.endpoint.isin(ep_list)]
    for ep_row in tqdm.tqdm(ep_sum.to_numpy()):
        if single:
            for genotype in gene_list:
                original_results = original_results.append(modeling(ep_row, [genotype]))
        else:
            original_results = original_results.append(modeling(ep_row, gene_list))
    original_results['genotype'] = original_results.index
    processed_results = original_results[['genotype']+key_cols[-3:]+key_cols[:4]]#.drop(columns='z')
    rename_dict = {key_cols[0]: 'coef', key_cols[1]: 'se', key_cols[3]: 'pval'}
    processed_results = processed_results.rename(columns=rename_dict)
    processed_results['or_025'] = np.exp(original_results[key_cols[4]])
    processed_results['or_975'] = np.exp(original_results[key_cols[5]])
    processed_results['or_500'] = np.exp(original_results[key_cols[0]])
    processed_results['p_sig'] = np.select([(processed_results.pval > 0.05/len(gene_list)),
                                            (processed_results.pval <= 0.05/len(gene_list))], [0, 1])
    processed_results = processed_results.merge(gene_sum_df, 'left', on='genotype')
    return processed_results


# find out the most significant hla genotypes that associated with T1D
# Method One:
res_df_1 = create_stat_df([OUTCOME], genes)
print(res_df_1[res_df_1.pval == res_df_1.pval.min()])
# Four genes are strongly associated with T1D (pval = 0.0): DQA1*03:01 DRB1*04:01 DQB1*03:02 DRB4*01:03
# 87 genotypes are significant
_ = plot_odds_ratio(res_df_1[res_df_1.p_sig == 1])
# Method Two: use conditional analysis
# in each loop, add the most significant hla genotype as a covariate
# genotypes are significant


# obtain stats for selected endpoints and save the results
res_df = create_stat_df(eps, genes)
res_df.to_csv('hla_res.csv', index=None)

# plot the associations between all the endpoints and each of the 4 most significant genes for T1D
for i in ['DQA1*03:01', 'DRB1*04:01', 'DQB1*03:02', 'DRB4*01:03']:
    _ = plot_odds_ratio(res_df[res_df.genotype == i])

# create a heatmap for the results
ep_gene_mat_sig = pd.pivot_table(res_df, values='p_sig', index=['genotype'], columns=['endpoint'])
# remove all the rows with only zeros
ep_gene_mat_sig = ep_gene_mat_sig[~(ep_gene_mat_sig == 0).all(axis=1)]
# replace all the non-sig cells from 0 to N/A
ep_gene_mat_sig = ep_gene_mat_sig.replace({0.:np.nan})
ep_gene_mat_coef = pd.pivot_table(res_df, values='coef', index=['genotype'], columns=['endpoint'])
# shrink the coef matrix to the size of sig matrix, so can go the next step
ep_gene_mat_coef = ep_gene_mat_coef[ep_gene_mat_coef.index.isin(ep_gene_mat_sig.index)]
# only keep those sig coefficients
ep_gene_mat_coef = ep_gene_mat_sig.multiply(ep_gene_mat_coef)
ep_gene_mat_or = np.exp(ep_gene_mat_coef)
ep_gene_mat_or = ep_gene_mat_or.round(2)
# # remove the rows that only have N/As or negative ORs
# ep_gene_mat_or['to_keep'] = ep_gene_mat_or.apply(lambda row:1 if len((row > 1).unique()) != 1 else 0, axis=1)
# ep_gene_mat_or = ep_gene_mat_or[ep_gene_mat_or.to_keep == 1]
# ep_gene_mat_or = ep_gene_mat_or.drop(columns=['to_keep'])
# use pval matrix to decide the color in the heatmap
ep_gene_mat_pval = pd.pivot_table(res_df, values='pval', index=['genotype'], columns=['endpoint'])
ep_gene_mat_pval = ep_gene_mat_pval[ep_gene_mat_pval.index.isin(ep_gene_mat_or.index)]
for i in ep_gene_mat_pval.columns:
    print(ep_gene_mat_pval[i].min()) # select -300 from this step for delta added to pval == 0.0
# transform the distribution of the pval matrix from (0, 1) to inverse (0, 1)
ep_gene_mat_imp = -np.log10(ep_gene_mat_pval + 10**(-300))
# remove those non-sig ones
ep_gene_mat_imp = ep_gene_mat_imp.multiply(ep_gene_mat_sig)
# transform the distribution of the matrix from (0.995, 1) to normal distributed (0, 1)
scaler = MinMaxScaler()
ep_gene_mat_imp = scaler.fit_transform(ep_gene_mat_imp)
ep_gene_mat_imp = pd.DataFrame(ep_gene_mat_imp, columns=ep_gene_mat_or.columns, index=ep_gene_mat_or.index)
# create a mask matrix for only those sig coefficients with positive/negative sign
ep_gene_mat_mask = np.select([(ep_gene_mat_coef > 0), (ep_gene_mat_coef < 0), (ep_gene_mat_coef.isna())],
                             [1, -1, np.nan])
ep_gene_mat_mask = pd.DataFrame(ep_gene_mat_mask, columns=ep_gene_mat_coef.columns, index=ep_gene_mat_coef.index)
ep_gene_mat_mask = ep_gene_mat_mask[ep_gene_mat_mask.index.isin(ep_gene_mat_or.index)]
# convert the distribution of the matrix from normal distributed (0, 1) to normal distributed (-1, 1)
ep_gene_mat_imp = ep_gene_mat_imp.multiply(ep_gene_mat_mask)
# keep only those rows with strong positive association with T1D
ep_gene_mat_or_pos = ep_gene_mat_or[ep_gene_mat_or[OUTCOME] > 1]
ep_gene_mat_imp_pos = ep_gene_mat_imp[ep_gene_mat_imp.index.isin(ep_gene_mat_or_pos.index)]
# create a heatmap for the results
sns.set(rc={'figure.figsize': (15, 8)})
ax = sns.heatmap(ep_gene_mat_imp_pos, linewidths=.5, center=0, annot=ep_gene_mat_or_pos, cmap='RdBu', fmt='')
ax.set_facecolor('#f7f7f7')
ax.collections[0].colorbar.set_label('z score')


# HLA PRS analysis
# split the data to training set and test set by 6:4
test_size = 0.4
seed = 4
hla_df_train = hla_df.sample(n=int(test_size * len(hla_df)), random_state=seed)
hla_df_test = hla_df[~hla_df.finngen_id.isin(hla_df_train.finngen_id)]


def filter_genotypes(dataset):
    conditions_to_add = []
    genotypes_to_check = genes
    while len(genotypes_to_check) > 0:
        stat_df = create_stat_df([OUTCOME], genotypes_to_check, dataset, conditions_to_add)
        stat_df = stat_df[stat_df.pval < 0.05/len(genes)]
        if len(stat_df) == 0:
            break
        genotypes_to_check = stat_df.gene.tolist()
        stat_df['filters'] = stat_df.coef/stat_df.se
        most_sig_genotype = stat_df[stat_df.filters == stat_df.filters.max()].gene.tolist()[0]
        conditions_to_add.append(most_sig_genotype)
        print('len_remained:', len(stat_df), '; selected:', conditions_to_add)
    results = create_stat_df([OUTCOME], conditions_to_add, hla_df_train)
    return results


res_train = filter_genotypes(hla_df_train)
res_all = filter_genotypes(hla_df)


