# !conda activate jupyter_env
# !cd familial_analysis
# !python3
import numpy as np
import pandas as pd
import tqdm, re
import statsmodels.api as sm
from basic_tools import load_data
from plot_tools import plot_odds_ratio
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, precision_score, recall_score, \
    precision_recall_curve, roc_curve

first_event_path = '/data/processed_data/endpointer/densified_first_events_DF8_all_endpoints_2021-09-04.txt'
info_path = '/data/processed_data/minimal_phenotype/minimal_phenotype_2022-03-28.csv'
pedigree_path = ''
df_events, df_info = load_data(first_event_path, info_path, pedigree_path)

# select the endpoints you plan to look into
# a list of ADs: https://risteys.finngen.fi/phenocode/AUTOIMMUNE
eps = ['T1D_STRICT', 'M13_RHEUMA', 'M13_RELAPSPOLYCHONDR', 'M13_SJOGREN', 'M13_SYSTSLCE', 'M13_DERMATOPOLY', # E4_DM1
       'M13_WEGENER', 'M13_MICROPOLYANG', 'M13_CHURGSTRAUSS', 'D3_ALLERGPURPURA', 'M13_BEHCET', 'M13_MCTD',
       'M13_HYPERANG', 'SLE_FG',  #'M13_SLE',
       'I9_RHEUFEV', 'G6_MS', 'G6_ADEM', 'G6_DISSOTH', 'G6_NARCOCATA', 'AUTOIMMUNE_HYPERTHYROIDISM',
       'E4_THYROIDITAUTOIM', 'E4_AUTOPOLYFAI', 'E4_HYTHY_AI_STRICT', 'E4_GRAVES_OPHT_STRICT', 'E4_ADDISON',
       'AUTOHEP', 'D3_AIHA_DRUG', 'D3_AIHA_OTHER', 'D3_ITP', 'D3_ANAEMIA_B12_DEF', 'K11_COELIAC', 'K11_IBD',
       'G6_MYASTHENIA', 'G6_OTHDEMYEL', 'G6_MYOMUSCINOTH', 'G6_GUILBAR', 'H7_IRIDOCYC_ANTER',  'CHIRBIL_PRIM',
       'L12_PSORIASIS', 'L12_VITILIGO', 'L12_ALOPECAREATA', 'L12_PEMPHIGOID', 'L12_DERMATHERP',
       'N14_HENOCHSCHONLEIN_NEPHRITIS', 'N14_IGA_NEPHROPATHY', 'T2D', 'GEST_DIABETES']

eps = ['T1D_STRICT', 'E4_THYROIDITAUTOIM', 'K11_COELIAC', 'D3_ANAEMIA_B12_DEF', 'M13_RHEUMA',
       'L12_VITILIGO', 'GRAVES_OPHT', 'K11_CROHN']
demos = ['sex']  # , 'ever_married', 'ISCED97', 'mother_tongue', 'post_code_first', 'number_of_children',
# 'in_social_assistance_registries', 'in_vaccination_registry', 'in_infect_dis_registry', 'in_malformations_registry',
#  'in_cancer_registry', 'ses', 'occupation', 'edulevel', 'edufield']

df = df_info[['FINREGISTRYID', 'ch_year', 'id_mother', 'id_father'] + demos]
df = df[(df.ch_year >= 1960.0) & (df.ch_year < 2000.0)]
df = df[(df.id_mother.isna() != True) & (df.id_father.isna() != True)]

if 'received_social_assistance' in demos:
    df['assisted'] = np.select([(~df.received_social_assistance.isna()), (df.received_social_assistance.isna())],
                               [1, 0])
    demos.remove("received_social_assistance")
    demos += ['assisted']
if 'ISCED97' in demos:
    df = df[~df.ISCED97.isna()]


# np.select faster & easier than df.merge faster than for loop
for i in tqdm.tqdm(range(len(eps))):
    df_events_sub = df_events[df_events.ENDPOINT == eps[i]]
    df['mo_ep' + str(i)] = np.select([
        (df['id_mother'].isin(df_events_sub.FINREGISTRYID)), (~df['id_mother'].isin(df_events_sub.FINREGISTRYID))
    ], [1, 0])
    df['fa_ep' + str(i)] = np.select([
        (df['id_father'].isin(df_events_sub.FINREGISTRYID)), (~df['id_father'].isin(df_events_sub.FINREGISTRYID))
    ], [1, 0])
    df['ch_ep' + str(i)] = np.select([
        (df['FINREGISTRYID'].isin(df_events_sub.FINREGISTRYID)),
        (~df['FINREGISTRYID'].isin(df_events_sub.FINREGISTRYID))
    ], [1, 0])

df.to_csv('data_all_big.csv', index=None)


x = df[[i for i in df.columns if re.match('ch_ep\d', i)]]
n_cases = x.sum()
n_cases_dict = dict(zip(eps, n_cases))
ep_remove = list(dict(filter(lambda item: item[1] < 200, n_cases_dict.items())).keys())

ep_cols = [i for i in df.columns if re.match('\w{2}_ep\d', i) and i not in ep_remove]

co_df = pd.DataFrame(columns=x.columns)
for i in x.columns:
    x_ = x[x[i] == 1]
    sums = np.round((x_.sum() / n_cases[i]) * 100, 1)
    co_df = co_df.append(sums, ignore_index=True)
co_df.columns = eps
co_df.index = eps

# number of parents
for i in range(len(eps)):
    df['pa_ep'+str(i)] = df['mo_ep'+str(i)]+df['fa_ep'+str(i)]
ep_cols = ['pa_ep'+str(i) for i in range(len(eps)) if eps[i] not in ep_remove]


# father and mother
ep_cols = ['mo_ep'+str(i) for i in range(len(eps)) if eps[i] not in ep_remove] + \
          ['fa_ep'+str(i) for i in range(len(eps)) if eps[i] not in ep_remove]

# split dataset into x,y
x_train = df[demos + ep_cols]
x_train = sm.add_constant(x_train)

df['outcome'] = df.ch_ep0
y_train = df['outcome']

# define model
try:
    lr = sm.Logit(y_train, x_train).fit()
except np.linalg.LinAlgError: # Singular matrix
    raise "any col having 0 only?"

stat = lr.summary2().tables[1]
# statistics
stat.columns = ['coef', 'std_err', 'z_stat', 'p_value', 'ci_1', 'ci_2']


def get_stat(ch_ep):
    df['outcome'] = df[ch_ep]
    y_train = df['outcome']
    lr = sm.Logit(y_train, x_train).fit()
    stat = lr.summary2().tables[1]
    stat.columns = ['coef', 'std_err', 'z_stat', 'p_value', 'ci_1', 'ci_2']
    return stat

get_stat('ch_ep0')


plot_odds_ratio(stat, demos, eps, 0)

# train-test split
df_train = df[(df.ch_year >= 1960) & (df.ch_year < 1995)]
df_valid = df[(df.ch_year >= 1995) & (df.ch_year < 2000)]
# split dataset into x,y
x_train = df_train[demos + ep_cols]
y_train = df_train['outcome']
x_valid = df_valid[demos + ep_cols]
y_valid = df_valid['outcome']

# define class weights
w = {0: 1, 1: 99}

# model 1
lr = LogisticRegression(class_weight=w).fit(x_train, y_train)
lr = LogisticRegression().fit(x_train, y_train)

# model 2
clf = RandomForestClassifier(max_depth=5, random_state=0, class_weight=w)
clf.fit(x_train, y_train)

# performance
y_pred = lr.predict(x_valid)
print(f'Accuracy Score: {accuracy_score(y_valid, y_pred)}')
print(f'Confusion Matrix: \n{confusion_matrix(y_valid, y_pred)}')
print(f'Area Under Curve: {roc_auc_score(y_valid, y_pred)}')
print(f'Precision score: {precision_score(y_valid, y_pred)}')
print(f'Recall score: {recall_score(y_valid, y_pred)}')
'''
Accuracy Score: 0.4901959450728078
Confusion Matrix:
[[142232 148673]
 [  1241   1916]]
Area Under Curve: 0.5479173327007315
Recall score: 0.6069052898321191
'''
