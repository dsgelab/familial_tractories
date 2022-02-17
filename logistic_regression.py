import pandas as pd
import numpy as np
import tqdm
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, precision_score, recall_score, \
    precision_recall_curve, roc_curve
from stat_tools import get_summary # from familial_analysis.stat_tools import get_summary
from plot_tools import plot_odds_ratio # from familial_analysis.plot_tools import plot_odds_ratio
from basic_tools import load_data
import logging

first_event_path = '/data/processed_data/endpointer/wide_first_events_endpoints_2021-09-04_densified.txt'
info_path = '/data/processed_data/minimal_phenotype/minimal_phenotype_file.csv'
df_events, df_info = load_data(first_event_path, info_path)

# select the endpoints you plan to look into
eps = ['T1D_WIDE', 'E4_THYROIDITAUTOIM', 'K11_COELIAC', 'D3_ANAEMIA_B12_DEF', 'K11_CHRONGASTR', 'L12_VITILIGO']
demos = ['sex', 'ever_married', 'received_social_assistance', 'ISCED97']


df = df_info[['FINREGISTRYID', 'ch_year', 'id_mother', 'id_father']+demos]
df = df[(df.ch_year >= 1960.0) & (df.ch_year < 2000.0)]
df = df[(df.id_mother.isna() != True) & (df.id_father.isna() != True)]

# np.select faster & easier than df.merge faster than for loop
for i in tqdm.tqdm(range(len(eps))):
    df_events_sub = df_events[df_events.ENDPOINT == eps[i]]
    df['mo_ep'+str(i)] = np.select([(df['id_mother'].isin(df_events_sub.FINNGENID)),
                                    (~df['id_mother'].isin(df_events_sub.FINNGENID))], [1, 0])
    df['fa_ep'+str(i)] = np.select([(df['id_father'].isin(df_events_sub.FINNGENID)),
                                    (~df['id_father'].isin(df_events_sub.FINNGENID))], [1, 0])
    if i == 0:
        df['outcome'] = np.select([(df['FINREGISTRYID'].isin(df_events_sub.FINNGENID)),
                                   (~df['FINREGISTRYID'].isin(df_events_sub.FINNGENID))], [1, 0])
    else:
        df['ch_ep'+str(i)] = np.select([(df['FINREGISTRYID'].isin(df_events_sub.FINNGENID)),
                                        (~df['FINREGISTRYID'].isin(df_events_sub.FINNGENID))], [1, 0])

if 'received_social_assistance' in demos:
    df['assisted'] = np.select([(~df.received_social_assistance.isna()), (df.received_social_assistance.isna())], [1, 0])
    demos.remove("received_social_assistance")
    demos += ['assisted']
if 'ISCED97' in demos:
    df = df[~df.ISCED97.isna()]

var_list = demos+['mo_ep'+str(i) for i in range(len(eps))]+['fa_ep'+str(i) for i in range(len(eps))]

# split dataset into x,y
x_train = df[var_list]
y_train = df['outcome']

# define class weights
w = {0: 1, 1: 99}
# define model
lr = LogisticRegression(class_weight=w)
# fit it
lr.fit(x_train, y_train)

# statistics
summary_df = get_summary(lr, x_train, y_train, var_list)
plot_odds_ratio(summary_df, eps, demos)


# train-test split
df_train = df[(df.ch_year >= 1960) & (df.ch_year < 1995)]
df_valid = df[(df.ch_year >= 1995) & (df.ch_year < 2000)]
# split dataset into x,y
x_train = df_train[var_list]
y_train = df_train['outcome']
x_valid = df_valid[var_list]
y_valid = df_valid['outcome']

# define class weights
w = {0: 1, 1: 99}

# model 1
lr = LogisticRegression(class_weight=w)
lr.fit(x_train, y_train)

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

