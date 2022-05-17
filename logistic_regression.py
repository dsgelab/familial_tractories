# !conda activate jupyter_env
# !cd familial_analysis
# !python3
import numpy as np
import pandas as pd
import tqdm, re
import statsmodels.api as sm
from basic_tools import eps
from plot_tools import plot_odds_ratio
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, precision_score, recall_score, \
    precision_recall_curve, roc_curve

data = pd.read_csv('data.csv')

demos = ['sex']  # , 'ever_married', 'mother_tongue', 'post_code_first', 'number_of_children',
# 'in_social_assistance_registries', 'in_vaccination_registry', 'in_infect_dis_registry', 'in_malformations_registry',
#  'in_cancer_registry', 'ses', 'occupation', 'edulevel', 'edufield']

df = data[['ID', 'sex', 'ch_year_range', 'fa_year_range', 'mo_year_range', 'sib_number'] + demos]

if 'received_social_assistance' in demos:
    df['assisted'] = np.select([(~df.received_social_assistance.isna()), (df.received_social_assistance.isna())],
                               [1, 0])
    demos.remove("received_social_assistance")
    demos += ['assisted']
if 'ISCED97' in demos:
    df = df[~df.ISCED97.isna()]



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
