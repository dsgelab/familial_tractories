import pandas as pd
import numpy as np
import datetime
import tqdm
import matplotlib.pyplot as plt
from data_exploration import get_matrix, get_length_heatmap

info_path = '/data/processed_data/minimal_phenotype/minimal_phenotype_file.csv'
all_event_path = '/data/processed_data/endpointer/longitudinal_endpoints_2021_12_20_no_OMITS.txt'

who_dict = {'ch':'child', 'mo':'mother', 'fa':'father'}

# Get all events
df_events = pd.read_csv(all_event_path)
# Get sex and approximate birth date of each indiv
df_info = pd.read_csv(info_path)

a = pd.read_csv('parent_child.csv')

test = df_events.ENDPOINT.value_counts()
# 3181 in total. The last 4 have value fewer than 5, so take the first 3177
ep_freq_df = pd.DataFrame({'ep':test.keys()[:3177].tolist(), 'freq':test.values[:3177]})

plt.rcParams['axes.grid'] = True

def get_prevalence_img(df, ep, base_year_start=1910, base_year_end=2010, base_who='ch'):
    '''
    df - pd.DataFrame
    who - string: 'ch', 'mo', 'fa'
    '''
    years = range(base_year_start, base_year_end)

    def get_data(df, who, base_year_start, base_year_end, base_who):
        coverage = []
        for num in years:  # tqdm.tqdm(range(base_year_start,base_year_end)):
            denominator = len(df[df[base_who + '_year'] == num])
            if denominator == 0:
                coverage.append(0.0)
            else:
                numerator = len(df[(df[who + '_age_start'].isnull() == False) & (df[base_who + '_year'] == num)])
                percent = numerator / denominator
                coverage.append(percent)
        return coverage

    ch_coverage = get_data(df, 'ch', base_year_start, base_year_end, base_who)
    mo_coverage = get_data(df, 'mo', base_year_start, base_year_end, base_who)
    fa_coverage = get_data(df, 'fa', base_year_start, base_year_end, base_who)

    max_coverage = max([max(ch_coverage), max(mo_coverage), max(fa_coverage)])
    max_tick = np.ceil(max_coverage * 10) / 10

    fig, axs = plt.subplots(3, figsize=(16, 9), sharex=True)

    for i, coverage, group in zip([0, 1, 2],
                                  [ch_coverage, mo_coverage, fa_coverage],
                                  ['Children', 'Mothers', 'Fathers']):
        axs[i].scatter(years, coverage)
        axs[i].set_ylabel('Prevalence', size=12)
        axs[i].set_title('Prevalence of ' + group, size=16)
        if (max_tick > 0.2) & (max_tick < 0.5):
            axs[i].set_yticks(np.arange(0.0, max_tick + 0.01, 0.1))
        elif max_tick >= 0.5:
            axs[i].set_yticks(np.arange(0.0, 1.01, 0.2))

    plt.suptitle('Prevalence of ' + ep + ' by birth year', size=20)
    plt.xlabel('Birth year of ' + who_dict[base_who], size=12)
    plt.show()

'''
'F5_DEMENTIA' very late onset
'I9_MI' late onset
'T2D' medium onset
'F5_SCHZPHR' early onset
'''

eps = ['E4_HYTHY_AI_STRICT',
 'T1D',
 'T1D_WIDE',
 'J10_ASTHMA',
 'L12_ATOPIC',
 'F5_DEPRESSIO',
 'F5_ALLANXIOUS',
       'K11_IBD']

ep = eps[7]

start = datetime.datetime.now()

df_events_t2d = df_events[df_events.ENDPOINT == ep]

b = df_info[['FINREGISTRYID']].merge(df_events_t2d[['FINREGISTRYID', 'EVENT_AGE']], 'left', on='FINREGISTRYID')
# b = b.sort_values('AGE')
b['dup_first'] = b.duplicated(subset='FINREGISTRYID', keep='first')
b['dup_last'] = b.duplicated(subset='FINREGISTRYID', keep='last')
age_start = b[b.dup_first == False]
age_end = b[b.dup_last == False]

df = a[['ch_id', 'ch_year', 'mo_id', 'mo_year', 'fa_id', 'fa_year']]

# add the start of the registry for each individual
df = df.merge(age_start[['FINREGISTRYID', 'EVENT_AGE']], how='left', left_on='ch_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE': 'ch_age_start'})
df = df.merge(age_start[['FINREGISTRYID', 'EVENT_AGE']], how='left', left_on='mo_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE': 'mo_age_start'})
df = df.merge(age_start[['FINREGISTRYID', 'EVENT_AGE']], how='left', left_on='fa_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE': 'fa_age_start'})

# add the end of the registry for each individual
df = df.merge(age_end[['FINREGISTRYID', 'EVENT_AGE']], how='left', left_on='ch_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE': 'ch_age_end'})
df = df.merge(age_end[['FINREGISTRYID', 'EVENT_AGE']], how='left', left_on='mo_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE': 'mo_age_end'})
df = df.merge(age_end[['FINREGISTRYID', 'EVENT_AGE']], how='left', left_on='fa_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE': 'fa_age_end'})

# calculate the length of registry coverage for each individual
df['ch_age_delta'] = df.ch_age_end - df.ch_age_start
df['fa_age_delta'] = df.fa_age_end - df.fa_age_start
df['mo_age_delta'] = df.mo_age_end - df.mo_age_start

df = df[['ch_id', 'ch_year', 'ch_age_start', 'ch_age_end', 'ch_age_delta',
         'mo_id', 'mo_year', 'mo_age_start', 'mo_age_end', 'mo_age_delta',
         'fa_id', 'fa_year', 'fa_age_start', 'fa_age_end', 'fa_age_delta']]

end = datetime.datetime.now()
print(end - start)
# 0:08:59.482246

# df = df[(df.ch_year < 2010)&(df.ch_year >= 1960)]

get_prevalence_img(df, ep, 1910, 2010, 'ch')
# get_prevalence_img(df, ep, 1880, 2000, 'mo')
# get_prevalence_img(df, ep, 1880, 2000, 'fa')

test = df[(df.ch_year >= 1910) & (df.ch_year < 2009)].dropna(thresh=7)
get_length_heatmap(get_matrix(test, 'ch'), get_matrix(test, 'mo'), get_matrix(test, 'fa'))