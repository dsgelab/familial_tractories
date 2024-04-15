import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from scipy import stats
scaler = StandardScaler
covariates = ['PC'+str(i) for i in range(1, 11)] + ['sex', 'BL_AGE', 'BL_YEAR']
eps = ['T1D_STRICT', 'E4_HYTHY_AI_STRICT', 'M13_RHEUMA', 'SLE_FG', 'K11_COELIAC', 'L12_PSORIASIS', 'K11_IBD', 'G6_MS']
eps_dict = {
    'T1D_STRICT': 'Type 1 diabetes',
    'E4_HYTHY_AI_STRICT': 'Autoimmune hypothyroidism',
    'M13_RHEUMA': 'Rheumatoid arthritis',
    'SLE_FG': 'Systemic lupus erythematosus',
    'K11_COELIAC': 'Coeliac disease',
    'L12_PSORIASIS': 'Psoriasis',
    'K11_IBD': 'Inflammatory bowel disease',
    'G6_MS': 'Multiple Sclerosis'
}

family = pd.read_csv('family.csv')
# len(family) = 12563
y = family['T1D_STRICT']

events = pd.read_csv('/finngen_R11/phenotype_1.0/data/finngen_R11_minimum_extended_1.0.txt.gz', sep='\t')
events = events[['FINNGENID', 'COHORT']].rename(columns={'FINNGENID':'finngen_id', 'COHORT':'source'})
family = family.merge(events, 'left')
tmp = family[(family.T1D_STRICT != 1)&(family.source.isin(['THL BIOBANK T1D']))].finngen_id
family = family[~family.finngen_id.isin(tmp)]
# len(family) = 10211

tmp = [i for i in family.columns.tolist() if i not in ['finngen_id', 'father_id', 'mother_id', 'source']+covariates+
       list(eps_dict.keys())+['fa_'+ep for ep in eps_dict.keys()]+['mo_'+ep for ep in eps_dict.keys()]]
for col in tmp:
    family[col] = scaler.fit_transform(family[col].values.reshape(-1,1))

r1 = {
    'T1D_STRICT': [0.67, 0.33],
    'E4_HYTHY_AI_STRICT': [0.31, 0.69],
    'M13_RHEUMA': [0.61, 0.39],
    'SLE_FG': [0.48, 0.52],
    'K11_COELIAC': [0.79, 0.21],
    'L12_PSORIASIS': [0.55, 0.45],
    'K11_IBD': [0.29, 0.71],
    'G6_MS': [0.57, 0.43]
}
for ep in eps:
    family['prsr1_'+ep] = (family['hla_'+ep])*(r1[ep][0]) + (family['non_'+ep])*(r1[ep][1])
    family['fa_prsr1_' + ep] = (family['fa_hla_' + ep]) * (r1[ep][0]) + (family['fa_non_' + ep]) * (r1[ep][1])
    family['mo_prsr1_' + ep] = (family['mo_hla_' + ep]) * (r1[ep][0]) + (family['mo_non_' + ep]) * (r1[ep][1])
    family['pa_prsr1_' + ep] = (family['mo_prsr1_' + ep] + family['fa_prsr1_' + ep]) / 2

r3 = {
    'T1D_STRICT': [0.67, 0.33],
    'E4_HYTHY_AI_STRICT': [0.48, 0.52],
    'M13_RHEUMA': [0.76, 0.24],
    'SLE_FG': [0.65, 0.35],
    'K11_COELIAC': [1.00, 0.00],
    'L12_PSORIASIS': [1.00, 0.00],
    'K11_IBD': [1.00, 0.00],
    'G6_MS': [1.00, 0.00]
}
for ep in eps:
    family['improved_'+ep] = (family['hla_'+ep])*(r3[ep][0]) + (family['non_'+ep])*(r3[ep][1])
    family['fa_improved_' + ep] = (family['fa_hla_' + ep]) * (r3[ep][0]) + (family['fa_non_' + ep]) * (r3[ep][1])
    family['mo_improved_' + ep] = (family['mo_hla_' + ep]) * (r3[ep][0]) + (family['mo_non_' + ep]) * (r3[ep][1])
    family['pa_improved_' + ep] = (family['mo_improved_' + ep] + family['fa_improved_' + ep]) / 2

FOLD_NUM = 5
family = family.sample(frac=1, random_state=5)
family.index = range(len(family))
family_n_folds = np.array_split(family, FOLD_NUM)

# different models predicting an AID via PGS for that AID at an individual level
for ep in eps:
    compare_list = [['prscs_'+ep], ['improved_'+ep], ['hla_'+ep, 'non_'+ep], ['hla_'+ep], ['non_'+ep]]
    groups = [[], [], [], [], []]

    for num in range(FOLD_NUM):
        test = family_n_folds[num]
        test.index = range(len(test))
        train = family_n_folds[~family.index.isin(test.index)]
        train.index = range(len(train))

        y_train = train[ep]
        y_test = test[ep]

        for var in compare_list:
            x_train = sm.add_constants(train[var])
            x_test = sm.add_constants(test[var])
            lr = sm.Logit(y_train, x_train).fit(disp=0)
            groups[compare_list.index(var)].append(round(metrics.roc_auc_score(y_test, lr.predict(x_test)), 3))

    print('individual', ep, 'AUC:', np.round(np.array(groups).mean(axis=1), 3))
    print('T-test:', stats.ttest_ind(groups[0], groups[1]).pvalue)

# different models predicting T1D via PGS for an AID at an individual level
for ep in eps:
    compare_list = [['prscs_'+ep], ['improved_'+ep], ['hla_'+ep, 'non_'+ep], ['hla_'+ep], ['non_'+ep]]
    groups = [[], [], [], [], []]

    for num in range(FOLD_NUM):
        test = family_n_folds[num]
        test.index = range(len(test))
        train = family_n_folds[~family.index.isin(test.index)]
        train.index = range(len(train))

        y_train = train['T1D_STRICT']
        y_test = test['T1D_STRICT']

        for var in compare_list:
            x_train = sm.add_constants(train[var])
            x_test = sm.add_constants(test[var])
            lr = sm.Logit(y_train, x_train).fit(disp=0)
            groups[compare_list.index(var)].append(round(metrics.roc_auc_score(y_test, lr.predict(x_test)), 3))

    print('individual', ep, 'AUC:', np.round(np.array(groups).mean(axis=1), 3))
    print('T-test:', stats.ttest_ind(groups[0], groups[1]).pvalue)

# different models predicting T1D in offspring via parental PGS for an AID
for ep in eps:
    compare_list = [['pa_prscs_'+ep], ['pa_improved_'+ep], ['pa_hla_'+ep, 'pa_non_'+ep], ['pa_hla_'+ep], ['pa_non_'+ep]]
    groups = [[], [], [], [], []]

    for num in range(FOLD_NUM):
        test = family_n_folds[num]
        test.index = range(len(test))
        train = family_n_folds[~family.index.isin(test.index)]
        train.index = range(len(train))

        y_train = train['T1D_STRICT']
        y_test = test['T1D_STRICT']

        for var in compare_list:
            x_train = sm.add_constants(train[var])
            x_test = sm.add_constants(test[var])
            lr = sm.Logit(y_train, x_train).fit(disp=0)
            groups[compare_list.index(var)].append(round(metrics.roc_auc_score(y_test, lr.predict(x_test)), 3))

    print('parental', ep, 'AUC:', np.round(np.array(groups).mean(axis=1), 3))
    print('T-test:', stats.ttest_ind(groups[0], groups[1]).pvalue)

# all the 8 parental AID PGSs vs PGS for parental T1D
auc_8aids, auc_t1d = [], []
for num in range(FOLD_NUM):
    test = family_n_folds[num]
    test.index = range(len(test))
    train = family_n_folds[~family.index.isin(test.index)]
    train.index = range(len(train))

    y_train = train['T1D_STRICT']
    y_test = test['T1D_STRICT']
    tmp = ['pa_improved_'+ep for ep in eps]
    x_train = sm.add_constants(train[tmp])
    x_test = sm.add_constants(test[tmp])
    lr = sm.Logit(y_train, x_train).fit(disp=0)
    auc_8aids.append(round(metrics.roc_auc_score(y_test, lr.predict(x_test)), 3))

    x_train = sm.add_constants(train[['pa_improved_T1D_STRICT']])
    x_test = sm.add_constants(test[['pa_improved_T1D_STRICT']])
    lr = sm.Logit(y_train, x_train).fit(disp=0)
    auc_t1d.append(round(metrics.roc_auc_score(y_test, lr.predict(x_test)), 3))

print('8 AIDs, AUC:', round(np.mean(auc_8aids), 3), '  complete res:', auc_8aids)
print('T1D, AUC:', round(np.mean(auc_t1d), 3), '  complete res:', auc_t1d)
print('T-test:', round(stats.ttest_ind(auc_8aids, auc_t1d).pvalue, 3))

