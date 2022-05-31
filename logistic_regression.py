# !conda activate jupyter_env
# !cd familial_analysis
# !python3
import tqdm
import json
import numpy as np
import pandas as pd
from statsmodels.discrete import conditional_models
from basic_tools import who_dict
from plot_tools import plot_odds_ratio, plot_crossed_odds_ratio

OUTCOME = 'T1D_STRICT'
data = pd.read_csv('data_'+OUTCOME+'.csv')
eps = json.load(open('eps_'+OUTCOME+'.json', 'r'))


def c_model(dataset, ep_index, who):
    ep_col_name = who_dict[who]+'_ep'+str(ep_index)
    try:
        endpoint = eps[ep_index]
        n_cases = dataset[ep_col_name].sum()
        lr = conditional_models.ConditionalLogit(endog=dataset.outcome,
                                                 exog=dataset[ep_col_name],
                                                 groups=dataset.subclass).fit(disp=0)
        pval = lr.pvalues[0]
        se = lr.bse[0]
        hr_025 = np.exp(lr.conf_int().iloc[0,0])
        hr_975 = np.exp(lr.conf_int().iloc[0,1])
        return [endpoint, who, se, pval, hr_025, hr_975, n_cases]
    except KeyError:
        return []


def model_loop(dataset, note, res_df):
    for i in tqdm.tqdm(range(len(eps))):
        res_fa = c_model(dataset, i, 'father')
        res_mo = c_model(dataset, i, 'mother')
        if res_fa:
            res_df = res_df.append(pd.Series(res_fa + [note], index=res.columns), ignore_index=True)
        if res_mo:
            res_df = res_df.append(pd.Series(res_mo + [note], index=res.columns), ignore_index=True)
    return res_df


res = pd.DataFrame(columns=["endpoint", "who", "se", "pval", "hr_025", "hr_975", "n_cases", "note"])
res = model_loop(data, 'all', res)
res = model_loop(data[data.sex == 1], 'women', res)
res = model_loop(data[data.sex == 0], 'men', res)
res.to_csv('res_'+OUTCOME+'.csv', index=None)

# plot the results
plot_odds_ratio(res[res.note == 'all'], eps, OUTCOME, bar_delta=.05)
plot_odds_ratio(res[res.note == 'women'], eps, OUTCOME, bar_delta=.05)
plot_odds_ratio(res[res.note == 'men'], eps, OUTCOME, bar_delta=.05)
plot_crossed_odds_ratio(res[res.pval < 0.05/len(eps)], ('women', 'men'), OUTCOME)
