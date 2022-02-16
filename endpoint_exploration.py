import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from data_exploration import pop #get_matrix, get_length_heatmap

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
        axs[i].axvline(x=1960, color='r', linestyle='--')
        axs[i].axvline(x=2000, color='r', linestyle='--')
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


def get_matrix(df, who):
    sub_df = df[[who + '_year', who + '_age_start']]
    sub_df = sub_df[sub_df.isna() == False]  # remove rows without registry records
    sub_df[who + '_age_delta'] = sub_df[who + '_age_start'].fillna(-0.5)

    sub_df[who + '_age_delta_bins'] = pd.cut(x=sub_df[who + '_age_delta'], bins=range(-1, 100),
                                             labels=range(-1, 99))  # max length in df is 97.xx -> 97
    matrix = pd.crosstab(sub_df[who + '_age_delta_bins'], sub_df[who + '_year'])

    # fill up missing rows and cols with 0
    missing_rows = set(range(-1, 100)).difference(matrix.index)
    missing_cols = set(np.arange(1880.0, 2010.0)).difference(matrix.columns)
    t = pd.DataFrame(np.zeros([len(missing_rows), len(matrix.columns)]), index=missing_rows, columns=matrix.columns)
    matrix = pd.concat([matrix, t], axis=0)
    matrix = pd.concat(
        [matrix, pd.DataFrame(np.zeros([101, len(missing_cols)]), index=range(-1, 100), columns=missing_cols)], axis=1)
    matrix = matrix.sort_index()  # sort rows
    matrix = matrix.sort_index(axis=1)  # sort cols

    return matrix


def get_length_heatmap(ch_mat, mo_mat, fa_mat, cut_year=1880.0, cut_count=99):
    '''
    Usage - create heatmaps for length of registry history by birth year
    ch_mat - pd.DataFrame: birth year of child by length of records
    mo_mat - pd.DataFrame: birth year of mother by length of records
    fa_mat - pd.DataFrame: birth year of father by length of records
    cut_year - float: e.g. 1910.0
    cut_count - int: e.g. 60
    '''

    ch_mat = ch_mat.loc[-1:cut_count, cut_year:2009]
    mo_mat = mo_mat.loc[-1:cut_count, cut_year:2009]
    fa_mat = fa_mat.loc[-1:cut_count, cut_year:2009]

    vals = [mo_mat, ch_mat, ch_mat, fa_mat]
    titles = ['Mother', 'Child', 'Child', 'Father']

    fig = plt.figure(figsize=(17, 17))

    grid = AxesGrid(fig, 111,
                    nrows_ncols=(2, 2),
                    axes_pad=0.3,
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single",
                    )

    for val, title, ax in zip(vals, titles, grid):
        im = ax.imshow(val, extent=(cut_year, 2009.0, cut_count, -1))
        ax.set_title(title)
        ax.set_xlabel('Birth year')
        ax.set_ylabel('Age of the first diagnosis')

    grid.cbar_axes[0].colorbar(im)

    for cax in grid.cbar_axes:
        cax.toggle_label(True)

    plt.suptitle('Heatmap of ' + ep + ' by birth year and age of the first diagnosis', size=15)
    plt.show()


ch_mat = get_matrix(df[(df.ch_year >= 1910) & (df.ch_year < 2009)], 'ch').iloc[1:,:]
mo_mat = get_matrix(df[(df.ch_year >= 1910) & (df.ch_year < 2009)], 'mo').iloc[1:,:]
fa_mat = get_matrix(df[(df.ch_year >= 1910) & (df.ch_year < 2009)], 'fa').iloc[1:,:]
get_length_heatmap(ch_mat, mo_mat, fa_mat)


test = df[(df.ch_year >= 1910) & (df.ch_year < 2009)].dropna(thresh=7)
get_length_heatmap(get_matrix(test, 'ch'), get_matrix(test, 'mo'), get_matrix(test, 'fa'))

# to check the incidence rate by age among all the population in Finland
df_all = ch_mat.sum(axis=1)
fig = plt.figure(figsize=(20,6))
plt.plot(df_all.keys(), df_all.values, 'o-')
plt.title('Incidence of '+ep+' by age at the first diagnosis',size=18)
plt.xlabel('Age at the first diagnosis',size=12)
plt.ylabel('Incidence',size=12)
plt.show()

# get all individuals in db for one heatmap
ch_mat = get_matrix(df[(df.ch_year >= 1910)], 'ch').iloc[1:,:]
plt.figure(figsize=(10,7))
plt.imshow(ch_mat, extent=(1920.0,2019.0,100,0))
plt.colorbar()
plt.xlabel('Birth year')
plt.ylabel('Age at the first diagnosis')
plt.suptitle('Heatmap of '+ep+' by birth year and age of the first diagnosis', size=15)
plt.show()


# get prevalence by groups 20, 40, 100
df['ch_age_bins'] = pd.cut(x=df['ch_age_start'], bins=[-1,20,40,100], labels=[20,40,100])
df['ch_year_start'] = np.floor(df.ch_year + df.ch_age_start)

df_test = df[['ch_id','ch_age_bins','ch_year_start']]#.dropna()
df_test = pd.DataFrame(df_test.groupby(['ch_age_bins', 'ch_year_start']).count().to_records())
df_test = df_test[df_test.ch_year_start >= 1960]

df_test = df_test.merge(pop, 'left', left_on='ch_year_start', right_on='year')
df_test['prevalence'] = df_test.ch_id/df_test.counts

plt.figure(figsize=(20,5))
plt.plot(df_test[df_test.ch_age_bins == 20].ch_year_start, df_test[df_test.ch_age_bins == 20].prevalence, 'o-', c='orange', label='0-19')
plt.plot(df_test[df_test.ch_age_bins == 40].ch_year_start, df_test[df_test.ch_age_bins == 40].prevalence, 'o-', c='green', label='20-39')
plt.plot(df_test[df_test.ch_age_bins == 100].ch_year_start, df_test[df_test.ch_age_bins == 100].prevalence, 'o-', c='purple', label='40-99')
plt.ylabel('Prevalence', size=12)
plt.xlabel('Year', size=12)
plt.title('Prevalence of '+ep+' over the years from 1960 to 2019', size=18)
plt.legend(loc='upper left', title='Age group')
plt.show()

# get incidence
plt.figure(figsize=(20,5))
plt.plot(df_test[df_test.ch_age_bins == 20].ch_year_start, df_test[df_test.ch_age_bins == 20].ch_id, 'o-', c='orange')
plt.plot(df_test[df_test.ch_age_bins == 40].ch_year_start, df_test[df_test.ch_age_bins == 40].ch_id, 'o-', c='green')
plt.plot(df_test[df_test.ch_age_bins == 100].ch_year_start, df_test[df_test.ch_age_bins == 100].ch_id, 'o-', c='purple')
plt.ylabel('Incidence', size=12)
plt.xlabel('Year', size=12)
plt.title('Incidence of '+ep+' over the years from 1960 to 2019', size=18)
plt.legend(loc='upper left', title='Age group')
plt.show()

# get overall prevalence
df_test1 = pd.DataFrame(df_test[['ch_year_start','ch_id','count']].groupby(['ch_year_start','count']).sum().to_records())
df_test1['prevalence'] = df_test1.ch_id/df_test1.counts
plt.figure(figsize=(20,5))
plt.plot(df_test1.ch_year_start, df_test1.prevalence, 'o-')
plt.ylabel('Prevalence', size=12)
plt.xlabel('Year', size=12)
plt.title('Prevalence of '+ep+' over the years from 1960 to 2019', size=18)
plt.show()
