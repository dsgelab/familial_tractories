import pandas as pd
import datetime

from data_exploration import get_coverage_img, get_matrix, get_length_heatmap


info_path = '/data/processed_data/minimal_phenotype/minimal_phenotype_file.csv'
all_event_path = '/data/processed_data/endpointer/longitudinal_endpoints_2021_12_20_no_OMITS.txt'

# Get all events
df_events = pd.read_csv(all_event_path)
# Get sex and approximate birth date of each indiv
df_info = pd.read_csv(info_path)

a = pd.read_csv('parent_child.csv')

test = df_events.ENDPOINT.value_counts()
# 3181 in total. The last 4 have value fewer than 5, so take the first 3177
ep_freq_df = pd.DataFrame({'ep':test.keys()[:3177].tolist(), 'freq':test.values[:3177]})

'''
'F5_DEMENTIA' very late onset
'I9_MI' late onset
'T2D' medium onset
'F5_SCHZPHR' early onset
'''
ep = 'F5_DEMENTIA'

start = datetime.datetime.now()

df_events_t2d = df_events[df_events.ENDPOINT == ep]

b = df_info[['FINREGISTRYID']].merge(df_events_t2d[['FINREGISTRYID','EVENT_AGE']],'left',on='FINREGISTRYID')
# b = b.sort_values('AGE')
b['dup_first'] = b.duplicated(subset='FINREGISTRYID', keep='first')
b['dup_last'] = b.duplicated(subset='FINREGISTRYID', keep='last')
age_start = b[b.dup_first == False]
age_end = b[b.dup_last == False]

df = a[['ch_id', 'ch_year', 'mo_id', 'mo_year', 'fa_id', 'fa_year']]

# add the start of the registry for each individual
df = df.merge(age_start[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='ch_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE':'ch_age_start'})
df = df.merge(age_start[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='mo_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE':'mo_age_start'})
df = df.merge(age_start[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='fa_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE':'fa_age_start'})


# add the end of the registry for each individual
df = df.merge(age_end[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='ch_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE':'ch_age_end'})
df = df.merge(age_end[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='mo_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE':'mo_age_end'})
df = df.merge(age_end[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='fa_id', right_on='FINREGISTRYID')
df = df.rename(columns={'EVENT_AGE':'fa_age_end'})

# calculate the length of registry coverage for each individual
df['ch_age_delta'] = df.ch_age_end - df.ch_age_start
df['fa_age_delta'] = df.fa_age_end - df.fa_age_start
df['mo_age_delta'] = df.mo_age_end - df.mo_age_start

df = df[['ch_id', 'ch_year', 'ch_age_start', 'ch_age_end', 'ch_age_delta',
         'mo_id', 'mo_year', 'mo_age_start', 'mo_age_end', 'mo_age_delta',
         'fa_id', 'fa_year', 'fa_age_start', 'fa_age_end', 'fa_age_delta']]

get_coverage_img(df, 'ch', 1910, 2010, 'ch')
get_coverage_img(df, 'mo', 1910, 2010, 'ch')
get_coverage_img(df, 'fa', 1910, 2010, 'ch')

get_coverage_img(df, 'ch', 1880, 2000, 'mo')
get_coverage_img(df, 'mo', 1880, 2000, 'mo')
get_coverage_img(df, 'fa', 1880, 2000, 'mo')

get_coverage_img(df, 'ch', 1880, 2000, 'fa')
get_coverage_img(df, 'mo', 1880, 2000, 'fa')
get_coverage_img(df, 'fa', 1880, 2000, 'fa')

end = datetime.datetime.now()
print(end - start)
# 0:08:59.482246


test = df[(df.ch_year >= 1910) & (df.ch_year < 2009)].dropna(thresh=7)
get_length_heatmap(get_matrix(test, 'ch'), get_matrix(test, 'mo'), get_matrix(test, 'fa'))