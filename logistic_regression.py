import pandas as pd
import numpy as np
import datetime
import tqdm
import re
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, precision_score, recall_score, \
    precision_recall_curve, roc_curve
from stat_tools import get_summary
from plot_tools import plot_odds_ratio

# select the endpoints you plan to look into
eps = ['T1D_WIDE', 'E4_THYROIDITAUTOIM', 'K11_COELIAC', 'D3_ANAEMIA_B12_DEF', 'K11_CHRONGASTR', 'L12_VITILIGO']

first_event_path = '/data/processed_data/endpointer/wide_first_events_endpoints_2021-09-04_densified.txt'
info_path = '/data/processed_data/minimal_phenotype/minimal_phenotype_file.csv'

# Get first events
df_events = pd.read_csv(first_event_path)

# Get sex and approximate birth date of each indiv
df_info = pd.read_csv(info_path)

df_info['ch_year'] = df_info['date_of_birth'].str.split('-').str[0].astype(float)
df = df_info[['FINREGISTRYID', 'ch_year', 'sex', 'ever_married',
                       'ever_divorced', 'received_social_assistance', 'ISCED97',
                       'eduyears_ISCED97', 'second_highest_income', 'id_mother', 'id_father']]
df = df[(df.ch_year >= 1960.0) & (df.ch_year < 2000.0)]
# 60.17% of individuals moved out of the minimal_phenotype file     7070745 -> 2816468

# print(len(df[df.mo_id.isna() == True]) / len(df))  # 0.08881691537059892
# print(len(df[df.fa_id.isna() == True]) / len(df))  # 0.11149141406896865

df = df[(df.id_mother.isna() != True) & (df.id_father.isna() != True)]
# 11.34% of children moved out of the cohort        2816468 -> 2497094
# incidence is 26052 and prevalence is 1.04%

# np.select faster & easier than df.merge faster than for loop
for i in tqdm.tqdm(range(len(eps))):
    df_events_sub = df_events[df_events.ENDPOINT == eps[i]]
    df['mo_ep'+str(i)] = np.select([
        (df['id_mother'].isin(df_events_sub.FINNGENID)),(~df['id_mother'].isin(df_events_sub.FINNGENID))
    ], [1, 0])
    df['fa_ep'+str(i)] = np.select([
        (df['id_father'].isin(df_events_sub.FINNGENID)),(~df['id_father'].isin(df_events_sub.FINNGENID))
    ], [1, 0])
    if i == 0:
        df['outcome'] = np.select([
            (df['FINREGISTRYID'].isin(df_events_sub.FINNGENID)), (~df['FINREGISTRYID'].isin(df_events_sub.FINNGENID))
        ], [1, 0])

# len(df[df.mother_tongue == 'fi']) / len(df)  # 93.5%
# len(df[(df.mother_tongue == 'fi') & (df.ch_t1d == True)]) / len(df[df.ch_t1d == True])  # 94.5%
#
# len(df[df.ever_married.isna() == True]) / len(df)  # 0%    no NA or NA is included into 0?
# len(df[df.ever_divorced.isna() == True]) / len(df)  # 0%   no NA or NA is included into 0?
# len(df[df.received_social_assistance.isna() == True]) / len(
#     df)  # only 57.86% 1s and NAs   NA is unknown or never received?
# len(df[df.second_highest_income.isna() == True]) / len(df)  # 48.95% missing
#
# len(df[df.ISCED97.isna() == True]) / len(df)  # 10.75% missing
# len(df[(df.ISCED97.isna() == True) & (df.ch_t1d == True)]) / len(df[df.ch_t1d == True])  # 12.81%

df['assisted'] = df.received_social_assistance.isna() == True

demos = ['sex']

# train-test split
df_train = df[(df.ch_year >= 1960) & (df.ch_year < 1995)]
df_valid = df[(df.ch_year >= 1995) & (df.ch_year < 2000)]
# split dataset into x,y
var_list = demos+['mo_ep'+str(i) for i in range(len(eps))]+['fa_ep'+str(i) for i in range(len(eps))]
x_train = df_train[var_list]
y_train = df_train['outcome']
x_valid = df_valid[var_list]
y_valid = df_valid['outcome']

# Control = len(cohort - cases) / len(subcohort & (cohort - cases))

# define class weights
w = {0: 1, 1: 99}
# define model
lr = LogisticRegression(class_weight=w)
# fit it
lr.fit(x_train, y_train)

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

# statistics
summary_df = get_summary(lr, x_train, y_train, var_list)
plot_odds_ratio(df, eps, demos)
