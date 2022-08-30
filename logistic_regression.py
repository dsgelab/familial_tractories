# !conda activate jupyter_env
# !cd familial_analysis
# !python3
import tqdm
import json
import numpy as np
import pandas as pd
from statsmodels.discrete import conditional_models
from basic_tools import who_dict
from plot_tools import plot_odds_ratio, plot_crossed_odds_ratio, draw_grouped_bar_plot

OUTCOME = 'M13_RHEUMA'
data = pd.read_csv('data_'+OUTCOME+'.csv')
eps = json.load(open('eps_'+OUTCOME+'.json', 'r'))
rename_dict = {}
for ep in eps:
    ep_index = eps.index(ep)
    ep_col_mo = 'mo_ep'+str(ep_index)
    ep_col_fa = 'fa_ep'+str(ep_index)
    ep_col_pa = 'pa_ep'+str(ep_index)
    rename_dict[ep_col_mo] = 'mo_' + ep
    rename_dict[ep_col_fa] = 'fa_' + ep
    rename_dict[ep_col_pa] = 'pa_' + ep
data = data.rename(columns=rename_dict)
data = data[['sex', 'subclass', 'outcome']+list(rename_dict.values())]

# to_remove ??
# ['E4_GRAVES_OPHT_STRICT', 'I9_RHEUFEV', 'L12_ALOPECAREATA', 'L12_DERMATHERP', 'L12_PEMPHIGOID', 'M13_MCTD']


def c_model(dataset, ep_col_name):
    try:
        n_cases = dataset[ep_col_name].sum()
        lr = conditional_models.ConditionalLogit(endog=dataset.outcome,
                                                 exog=dataset[ep_col_name],
                                                 groups=dataset.subclass).fit(disp=0)
        pval = lr.pvalues[0]
        se = lr.bse[0]
        or_025 = np.exp(lr.conf_int().iloc[0, 0])
        or_975 = np.exp(lr.conf_int().iloc[0, 1])
        subclass_list = []
        for i in range(data.subclass.max() + 1):  # 4 secs for a loop
            if len(data[data.subclass == i][ep_col_name].unique()) > 1:
                subclass_list.append(i)
        n_valid_pair00 = len(dataset[(dataset.subclass.isin(subclass_list)) &
                                     (dataset.outcome == 0) & (dataset[ep_col_name] == 0)])
        n_valid_pair01 = len(dataset[(dataset.subclass.isin(subclass_list)) &
                                     (dataset.outcome == 0) & (dataset[ep_col_name] == 1)])
        n_valid_pair10 = len(dataset[(dataset.subclass.isin(subclass_list)) &
                                     (dataset.outcome == 1) & (dataset[ep_col_name] == 0)])
        n_valid_pair11 = len(dataset[(dataset.subclass.isin(subclass_list)) &
                                     (dataset.outcome == 1) & (dataset[ep_col_name] == 1)])
        return [se, pval, or_025, or_975, len(subclass_list), n_cases,
                n_valid_pair00, n_valid_pair01, n_valid_pair10, n_valid_pair11]
    except ValueError:
        pass


def model_loop(dataset, who, note, res_df, endpoints=eps):
    for endpoint in tqdm.tqdm(endpoints):
        ep_col_name = who_dict[who] + '_' + endpoint
        if ep_col_name in dataset:
            result = c_model(dataset, ep_col_name)
        res_df = res_df.append(pd.Series([endpoint, note, who] + result, index=res_df.columns), ignore_index=True)
    return res_df


res = pd.DataFrame(columns=["endpoint", "note", "who", "se", "pval", "or_025", "or_975", 'n_valid_group', "n_cases",
                            'n_valid_pair00', 'n_valid_pair01', 'n_valid_pair10', 'n_valid_pair11'])
res = model_loop(data, 'parent', 'all', res)
eps_sig = plot_odds_ratio(res[res.who == 'parent'], 'T1D_STRICT')
# comparisons between fathers and mothers
res = model_loop(data, 'father', 'all', res, eps_sig)
res = model_loop(data, 'mother', 'all', res, eps_sig)
# comparisons between sons and daughters
res = model_loop(data[data.sex == 0], 'parent', 'boy', res, eps_sig)
res = model_loop(data[data.sex == 1], 'parent', 'girl', res, eps_sig)
# belows are too complex to understand
# comparisons between fathers and mothers for sons and daughters separately
# res = model_loop(data[data.sex == 0], 'father', 'boy', res, eps_sig)
# res = model_loop(data[data.sex == 1], 'father', 'girl', res, eps_sig)
# res = model_loop(data[data.sex == 0], 'mother', 'boy', res, eps_sig)
# res = model_loop(data[data.sex == 1], 'mother', 'girl', res, eps_sig)
# summary statistics after removing those whose parents were born after 1961
data_sub = data[(data.fa_year_range < 1960) & (data.mo_year_range < 1960)]
res = model_loop(data_sub, 'all_sub', res)
res = model_loop(data_sub[data_sub.sex == 1], 'women_sub', res)
res = model_loop(data_sub[data_sub.sex == 0], 'men_sub', res)
# save the results
res.to_csv('res_'+OUTCOME+'.csv', index=None)

# plot the results
plot_odds_ratio(res[res.note == 'all'], eps, OUTCOME, bar_delta=.05)
plot_odds_ratio(res[res.note == 'women'], eps, OUTCOME, bar_delta=.05)
plot_odds_ratio(res[res.note == 'men'], eps, OUTCOME, bar_delta=.05)
plot_crossed_odds_ratio(res, ('women', 'men'), OUTCOME)

plot_odds_ratio(res[res.note == 'all_sub'], eps, OUTCOME, bar_cap=.05)
plot_crossed_odds_ratio(res, ('all', 'all_sub'), OUTCOME)

draw_grouped_bar_plot(res[res.note == 'all'], 'who', 'n_cases', title='All individuals')
draw_grouped_bar_plot(res[res.note == 'all_sub'], 'who', 'n_cases', title='Individuals whose parents were born < 1962')

