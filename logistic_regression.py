import numpy as np
import pandas as pd
import tqdm, re
import statsmodels.api as sm
from basic_tools import load_data
from stat_tools import vif
from plot_tools import plot_odds_ratio
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, precision_score, recall_score, \
    precision_recall_curve, roc_curve

first_event_path = '/data/processed_data/endpointer/wide_first_events_endpoints_2021-09-04_all_densified.txt'
info_path = '/data/processed_data/minimal_phenotype/minimal_phenotype_2022-02-17.csv'
df_events, df_info = load_data(first_event_path, info_path)

# select the endpoints you plan to look into

eps = ['T1D_STRICT', 'E4_THYROIDITAUTOIM', 'K11_COELIAC', 'D3_ANAEMIA_B12_DEF', 'M13_RHEUMA',
       'L12_VITILIGO', 'GRAVES_OPHT', 'K11_CROHN']
demos = ['sex']#, 'ever_married', 'received_social_assistance', 'ISCED97']

df = df_info[['FINREGISTRYID', 'ch_year', 'id_mother', 'id_father']+demos]
df = df[(df.ch_year >= 1960.0) & (df.ch_year < 2000.0)]
df = df[(df.id_mother.isna() != True) & (df.id_father.isna() != True)]

if 'received_social_assistance' in demos:
    df['assisted'] = np.select([(~df.received_social_assistance.isna()), (df.received_social_assistance.isna())], [1, 0])
    demos.remove("received_social_assistance")
    demos += ['assisted']
if 'ISCED97' in demos:
    df = df[~df.ISCED97.isna()]

# METHOD 1
# np.select faster & easier than df.merge faster than for loop
for i in tqdm.tqdm(range(len(eps))):
    df_events_sub = df_events[df_events.ENDPOINT == eps[i]]
    if i == 0:
        df['outcome'] = np.select([
            (df['FINREGISTRYID'].isin(df_events_sub.FINNGENID)), (~df['FINREGISTRYID'].isin(df_events_sub.FINNGENID))
        ], [1, 0])
    df = df.merge(df_events_sub[['FINREGISTRYID', 'YEAR']], 'left')
    df['ch_ep'+str(i)] = df.YEAR.fillna(2022.)
    df = df.drop(columns=['FINREGISTRYID', 'YEAR'])
    df = df.merge(df_events_sub[['FINREGISTRYID', 'YEAR']], 'left', left_on='id_mother', right_on='FINREGISTRYID')
    df['mo_ep'+str(i)] = df.YEAR.fillna(2022.)
    df = df.drop(columns=['FINREGISTRYID', 'YEAR'])
    df = df.merge(df_events_sub[['FINREGISTRYID', 'YEAR']], 'left', left_on='id_father', right_on='FINREGISTRYID')
    df['fa_ep'+str(i)] = df.YEAR.fillna(2022.)
    df = df.drop(columns=['FINREGISTRYID', 'YEAR'])

# ep_cols = ['mo_ep'+str(i) for i in range(len(eps))] +\
#           ['fa_ep'+str(i) for i in range(len(eps))] +\
#           ['ch_ep'+str(i) for i in range(1, len(eps))]
ep_cols = [i for i in df.columns if re.match('\w{2}_ep\d', i) and i != 'ch_ep0']

for ep in ep_cols:
    df[ep] = np.select([(df[ep] < df.ch_ep0), (df[ep] >= df.ch_ep0)], [1, 0])

# METHOD 2
# np.select faster & easier than df.merge faster than for loop
for i in tqdm.tqdm(range(len(eps))):
    df_events_sub = df_events[df_events.ENDPOINT == eps[i]]
    df['mo_ep'+str(i)] = np.select([
        (df['id_mother'].isin(df_events_sub.FINREGISTRYID)),(~df['id_mother'].isin(df_events_sub.FINREGISTRYID))
    ], [1, 0])
    df['fa_ep'+str(i)] = np.select([
        (df['id_father'].isin(df_events_sub.FINREGISTRYID)),(~df['id_father'].isin(df_events_sub.FINREGISTRYID))
    ], [1, 0])
    df['ch_ep'+str(i)] = np.select([
        (df['FINREGISTRYID'].isin(df_events_sub.FINREGISTRYID)), (~df['FINREGISTRYID'].isin(df_events_sub.FINREGISTRYID))
    ], [1, 0])

df.to_csv('data_all.csv', index=None)

ep_cols = [i for i in df.columns if re.match('\w{2}_ep\d', i)]

x = df[[i for i in df.columns if re.match('ch_ep\d', i)]]
n_cases = x.sum()
co_df = pd.DataFrame(columns=x.columns)
for i in x.columns:
    x_ = x[x[i] == 1]
    sums = np.round((x_.sum() / n_cases[i]) * 100, 1)
    co_df = co_df.append(sums, ignore_index=True)
co_df.columns = eps
co_df.index = eps

df = df.rename(columns={'ch_ep1': 'outcome'})
ep_cols = ['mo_ep'+str(i) for i in range(len(eps))] +\
          ['fa_ep'+str(i) for i in range(len(eps))]

# split dataset into x,y
x_train = df[demos+ep_cols]
x_train = sm.add_constant(x_train)
y_train = df['outcome']

# define model
lr = sm.Logit(y_train, x_train).fit()
stat = lr.summary2().tables[1]
# statistics
stat.columns = ['coef', 'std_err', 'z_stat', 'p_value', 'ci_1', 'ci_2']

df = df.rename(columns={'outcome': 'ch_ep7'})
# df = df.rename(columns={'ch_ep7': 'outcome'})
stat

plot_odds_ratio(stat, demos, eps, 0)


# train-test split
df_train = df[(df.ch_year >= 1960) & (df.ch_year < 1995)]
df_valid = df[(df.ch_year >= 1995) & (df.ch_year < 2000)]
# split dataset into x,y
x_train = df_train[demos+ep_cols]
y_train = df_train['outcome']
x_valid = df_valid[demos+ep_cols]
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

