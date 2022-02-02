import pandas as pd
import numpy as np

import json
from collections import defaultdict
import datetime
import tqdm
import re

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

first_event_path = '/data/processed_data/endpointer/wide_first_events_endpoints_2021-09-04_densified.txt'
info_path = '/data/processed_data/minimal_phenotype/minimal_phenotype_file.csv'
all_event_path = '/data/processed_data/endpointer/longitudinal_endpoints_2021_12_20_no_OMITS.txt'

who_dict = {'ch':'child', 'mo':'mother', 'fa':'father'}


# Get first events
df_events = pd.read_csv(first_event_path)

# Get sex and approximate birth date of each indiv
df_info = pd.read_csv(info_path)

# Get all events
df_events = pd.read_csv(all_event_path)


# create a table for mother, father and child
start = datetime.datetime.now()
a = df_info[['FINREGISTRYID','date_of_birth','id_mother','id_father']].merge(df_info[['FINREGISTRYID','date_of_birth']],'left',left_on='id_mother',right_on='FINREGISTRYID')
a = a.rename(columns={"FINREGISTRYID_y": "mo_id", "date_of_birth_y": "mo_birth"})
a = a.merge(df_info[['FINREGISTRYID','date_of_birth']],'left',left_on='id_father',right_on='FINREGISTRYID')
a = a.rename(columns={"FINREGISTRYID_x": "ch_id", "date_of_birth_x": "ch_birth","FINREGISTRYID": "fa_id", "date_of_birth": "fa_birth"})
end = datetime.datetime.now()
print(end - start)

def get_coverage_img(df, who):
    '''
    df - pd.DataFrame
    who - string: 'ch', 'mo', 'fa'
    '''

    coverage, count = [], []
    for num in tqdm.tqdm(range(1910,2010)):
        numerator = len(df[(df[who+'_age_start'].isnull() == False)&(df.ch_year == num)])
        denominator = len(df[df.ch_year == num])
        percent = numerator/denominator
        coverage.append(percent)
        count.append(numerator)

    fig, ax = plt.subplots(figsize=(20,8))
    ax.grid()
    scatter = ax.scatter(range(1910,2010), coverage, s=np.log10(count)*10)
    handles, labels = scatter.legend_elements(prop="sizes", alpha=0.7)
    
    # display the right numbers in the count legend
    labels_ = []
    for i in labels:
        num = re.findall('\{([\d\.]+)\}',i)[0]
        labels_.append(re.sub(num,str(int(np.ceil(10**(float(num)/10)))),i))

    legend = ax.legend(handles, labels_, loc="lower right", title="Count")
    ax.add_artist(legend)

    ax.set_ylabel('Registry coverage', size=14)
    ax.set_xlabel('Birth year of child', size=14)
    ax.set_title('Registry coverage by birth year - '+who_dict[who], size=20)
    ax.set_yticks(np.arange(0.0,1.0,0.2))
    plt.show()

get_coverage_img(a, 'ch')
get_coverage_img(a, 'mo')
get_coverage_img(a, 'fa')


df = a[(a.ch_year >= 1960)&(a.ch_year < 2010)][['ch_birth', 'mo_birth', 'fa_birth', 'ch_age',
       'mo_age', 'fa_age', 'ch_year', 'mo_year', 'fa_year']]

df['ch_yr_bins'] = pd.cut(x=df['ch_year'], bins=[1959, 1969, 1979, 1989, 1999, 2009], labels=[1,2,3,4,5])
#                           labels=['1960s','1970s','1980s','1990s','2000s'])

df['fa_age_bins'] = pd.cut(x=df['fa_age'], bins=[-1,5,10,15,20,25,30,35,40,45,50,55,60,150], 
                           labels=[2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,80]).values.add_categories(-20)
df['ch_age_bins'] = pd.cut(x=df['ch_age'], bins=[-1,5,10,15,20,25,30,35,40,45,50,55,60,150], 
                           labels=[2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,80]).values.add_categories(-20)
df['mo_age_bins'] = pd.cut(x=df['mo_age'], bins=[-1,5,10,15,20,25,30,35,40,45,50,55,60,150], 
                           labels=[2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,80]).values.add_categories(-20)
df['fa_age_bins'] = df['fa_age_bins'].fillna(-20)
df['ch_age_bins'] = df['ch_age_bins'].fillna(-20)
df['mo_age_bins'] = df['mo_age_bins'].fillna(-20)

size_labels = ['$\\mathdefault{10}$','$\\mathdefault{100}$','$\\mathdefault{1000}$','$\\mathdefault{10000}$','$\\mathdefault{100000}$']
age_labels = ['$\\mathdefault{1960s}$','$\\mathdefault{1970s}$','$\\mathdefault{1980s}$','$\\mathdefault{1990s}$','$\\mathdefault{2000s}$']

father = pd.DataFrame(df[['fa_year','fa_age_bins','ch_yr_bins','fa_birth']]
                    .groupby(['fa_year','fa_age_bins','ch_yr_bins']).count().to_records())
father = father.rename(columns={'fa_birth':'counts'})
father = father[father.counts != 0]
father['fa_count'] = pd.cut(x=father['counts'], bins=[0,10,100,1000,10000,100000],
                          labels=[10,100,1000,10000,100000])

fig, ax = plt.subplots(figsize=(20,8))
ax.grid()
scatter = ax.scatter(father.fa_year, father.fa_age_bins, c=father.ch_yr_bins, s=np.log(father.fa_count.astype(int))*30, alpha=0.6)

handles1, _ = scatter.legend_elements(prop="colors", alpha=1.0)
legend1 = ax.legend(handles1, age_labels, loc="center left", title="Child age")
ax.add_artist(legend1)

handles2, _ = scatter.legend_elements(prop="sizes", alpha=0.5)#, labels=['1','2','3','4','5'])
legend2 = ax.legend(handles2, size_labels, loc="upper left", title="Count")
ax.set_ylabel('Age', size=14)
ax.set_xlabel('Birth year', size=14)
ax.set_title('Age distribution of the first events - fathers', size=20)
ax.set_xticks(np.arange(1860,2020,10))
# ax.set_xticklabels(bins,rotation=90)
plt.show()


# load all events

start = datetime.datetime.now()
b = df_info[['FINREGISTRYID']].merge(df_events[['FINREGISTRYID','EVENT_AGE']],'left',on='FINREGISTRYID')
# b = b.sort_values('AGE')
b['dup'] = b.duplicated(subset='FINREGISTRYID', keep='last')
age_df = b[b.dup == False]
end = datetime.datetime.now()
print(end - start)
# 0:08:59.482246

a = a.rename(columns={'ch_age':'ch_age_start','fa_age':'fa_age_start','mo_age':'mo_age_start'})

# add the end of the registry for each individual
a = a.merge(age_df[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='ch_id', right_on='FINREGISTRYID')
a = a.rename(columns={'EVENT_AGE':'ch_age_end'})
a = a.merge(age_df[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='mo_id', right_on='FINREGISTRYID')
a = a.rename(columns={'EVENT_AGE':'mo_age_end'})
a = a.merge(age_df[['FINREGISTRYID','EVENT_AGE']], how='left', left_on='fa_id', right_on='FINREGISTRYID')
a = a.rename(columns={'EVENT_AGE':'fa_age_end'})

# calculate the length of registry coverage for each individual
a['ch_age_delta'] = a.ch_age_end - a.ch_age_start
a['fa_age_delta'] = a.fa_age_end - a.fa_age_start
a['mo_age_delta'] = a.mo_age_end - a.mo_age_start

a = a[['ch_id', 'ch_birth', 'ch_year', 'ch_age_start', 'ch_age_end', 'ch_age_delta',
       'mo_id', 'mo_birth', 'mo_year', 'mo_age_start', 'mo_age_end', 'mo_age_delta',
       'fa_id', 'fa_birth', 'fa_year', 'fa_age_start', 'fa_age_end', 'fa_age_delta']]

a['ch_age_delta'] = a.ch_age_end - a.ch_age_start
a['fa_age_delta'] = a.fa_age_end - a.fa_age_start
a['mo_age_delta'] = a.mo_age_end - a.mo_age_start

df = a[(a.ch_year <= 1910) & (a.ch_year > 2010)]

#create matrix for heatmap
def get_matrix(df, who):
	'''
	df - pd.DataFrame
	who - string: 'ch', 'mo', 'fa'
	'''

    # select a subset we need from dataframe according to who
    sub_df = df[[who+'_year', who+'_age_delta']]
    sub_df = sub_df[sub_df.isna() == False] # remove rows without registry records
    sub_df[who+'_age_delta'] = sub_df[who+'_age_delta'].fillna(-0.5)
    
    sub_df[who+'_age_delta_bins'] = pd.cut(x=sub_df[who+'_age_delta'], bins=range(-1,99), labels=range(-1,98)) # max length in df is 97.xx -> 97
    matrix = pd.crosstab(sub_df[who+'_age_delta_bins'], sub_df[who+'_year'])
    
    # fill up missing rows and cols with 0
    missing_rows = set(range(-1,98)).difference(matrix.index)
    missing_cols = set(np.arange(1870.0,2010.0)).difference(matrix.columns)
    t = pd.DataFrame(np.zeros([len(missing_rows),len(matrix.columns)]), index=missing_rows, columns=matrix.columns)
    matrix = pd.concat([matrix,t], axis=0)
    matrix = pd.concat([matrix,pd.DataFrame(np.zeros([99,len(missing_cols)]), index=range(-1,98), columns=missing_cols)], axis=1)
    matrix = matrix.sort_index() # sort rows
    matrix = matrix.sort_index(axis=1) # sort cols
    
    return matrix

ch_mat = get_matrix(df, 'ch')
mo_mat = get_matrix(df, 'mo')
fa_mat = get_matrix(df, 'fa')

# create heatmaps for length of registry history by birth year
vals = [mo_mat, ch_mat, ch_mat, fa_mat]
titles = ['Mother', 'Child', 'Child', 'Father']

fig = plt.figure(figsize=(17,17))

grid = AxesGrid(fig, 111,
                nrows_ncols=(2, 2),
                axes_pad=0.3,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )

for val, title, ax in zip(vals,titles,grid):
    im = ax.imshow(val, extent=(1870.0, 2009.0, 97, -1))
    ax.set_title(title)
    ax.set_xlabel('Birth year')
    ax.set_ylabel('Length of records')

grid.cbar_axes[0].colorbar(im)

for cax in grid.cbar_axes:
    cax.toggle_label(True)
    
plt.show()


def get_length_img(df, who):
	'''
	df - pd.DataFrame
	who - string: 'ch', 'mo', 'fa'
	'''

	# create a sub dataframe for heatmap
	sub_df = df[[who+'_year', who+'_age_delta']]
	sub_df = sub_df[(sub_df[who+'_year'] >= 1910) & (sub_df[who+'_year'] < 2010)]
	sub_df[who+'_age_delta'] = sub_df[who+'_age_delta'].fillna(-0.5)
	sub_df[who+'_age_delta_bins'] = pd.cut(x=sub_df[who+'_age_delta'], bins=range(-1,99), labels=range(-1,98)) # max length in df is 97.xx
	matrix = pd.crosstab(sub_df[who+'_year'], sub_df[who+'_age_delta_bins'])

	# create a heatmap for length of registry history by birth year
	plt.figure(figsize=(10,10))
	extent = sub_df[who+'_year'].min(), sub_df[who+'_year'].max(), sub_df[who+'_age_delta_bins'].min(), sub_df[who+'_age_delta_bins'].max()
	plt.imshow(matrix, extent=extent)
	plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
	plt.xlabel('Birth year of '+who_dict[who], size=14)
	plt.ylabel('Length of registry coverage', size=14)
	plt.title('Length of registry history by birth year - '+who_dict[who], size=20)
	plt.colorbar(cax=plt.axes([0.85, 0.1, 0.02, 0.8]))
	plt.show()