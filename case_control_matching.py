"""
METHOD 1: exact matching
@Author: feiyiwang
"""

import re
import tqdm
import json
import itertools
from functools import reduce
import pandas as pd
from basic_tools import eps
from plot_tools import draw_distribution
from stat_tools import get_cohend, sum_cases

# define matching ratio and matching factors
OUTCOME = 'M13_RHEUMA'
MATCH_NUM = 3
MATCH_FACTORS = ['sex', 'ch_year_range', 'fa_year_range', 'mo_year_range', 'sib_number', 'province']


def case_control_matching(dataset, outcome, matching_factors, matching_number):
    """
    :param dataset: a DataFrame of study population
    :param outcome: a string of a given disease
    :param matching_factors: a list of matching factors
        e.g. ['sex','ch_year_range','fa_year_range','mo_year_range','sib_number', 'province']
    :param matching_number: an int indicating how many controls we need to select from the cohort
    :return: a DataFrame of cases; a DataFrame of controls
    """
    # split the data into cases and cohort for controls
    ep_index = str(eps.index(outcome))
    cases = dataset[~dataset['ch_age'+ep_index].isna()]  # 15102
    control_for_match = dataset[dataset['ch_age'+ep_index].isna()]
    print('To begin with, we split the study population into', len(cases), 'potential cases and',
          len(control_for_match), 'potential controls.')

    # list all possible permutations using the selected matching factors
    match_var = [dataset[col].unique().tolist() for col in matching_factors]
    match_permutation = list(itertools.product(*match_var))  # 3296000 -> 10560 -> 1734
    print(len(match_permutation), 'potential permutations were found.')

    # filter out all the impossible permutations
    match_permutation_keep = []
    for i in tqdm.tqdm(match_permutation):
        try:
            conditions = [(cases[matching_factors[j]] == i[j]) for j in range(len(matching_factors))]
            match = cases[reduce(lambda x, y: x & y, conditions)]
            if len(match) != 0:
                match_permutation_keep.append({'conditions': i, 'id_list': match.ID.tolist()})
        except Exception as e:
            print(e, i)
    print(len(match_permutation_keep), 'possible permutations were found.')

    # filter out all the permutations that cannot find enough controls
    match_failed, matched_ids = [], {}
    for i in match_permutation_keep:
        try:
            conditions = [(control_for_match[matching_factors[j]] == i['conditions'][j])
                          for j in range(len(matching_factors))]
            potential_control = control_for_match[reduce(lambda x, y: x & y, conditions)]
            if len(potential_control) < matching_number * len(i['id_list']):
                match_failed.append(i)
            else:
                potential_control = potential_control.sample(n=matching_number * len(i['id_list']))
                index_num = match_permutation_keep.index(i)
                for each in potential_control.ID.tolist()+i['id_list']:
                    matched_ids[each] = index_num
        except:
            match_failed.append(i)
    print(len(match_failed), 'permutations failed to match to enough controls.')

    # remove all the failed permutations from cases
    cases_to_remove = []
    for i in tqdm.tqdm(match_failed):
        cases_to_remove += i['id_list']
    print(len(cases_to_remove), 'permutations that could not be used to match enough controls were removed from cases.')

    # final cases and controls
    matched_data = dataset.merge(pd.DataFrame(list(matched_ids.items()), columns=['ID', 'subclass']), 'right', on='ID')
    print('Finally, we have', len(matched_data), 'rows in the matched data.')

    return matched_data


def remove_unnecessary_columns(dataset, outcome, threshold=20):
    n_cases_mo, ep_remove_mo = sum_cases(dataset, 'mother', threshold)
    n_cases_fa, ep_remove_fa = sum_cases(dataset, 'father', threshold)
    remained_eps = list(set(eps).difference(set(ep_remove_mo).intersection(set(ep_remove_fa))))
    remained_eps = [i for i in eps if i in remained_eps]
    dataset = dataset.drop(columns=['fa_ep' + str(eps.index(i)) for i in ep_remove_fa] +
                                   ['mo_ep' + str(eps.index(i)) for i in ep_remove_mo] +
                                   ['ch_ep' + str(eps.index(i)) for i in eps if i != outcome] +
                                   ['fa_age' + str(eps.index(i)) for i in ep_remove_fa] +
                                   ['mo_age' + str(eps.index(i)) for i in ep_remove_mo] +
                                   ['ch_age' + str(eps.index(i)) for i in eps if i != outcome])
    col_dict = {'ch_ep' + str(eps.index(outcome)): 'outcome', 'ch_age' + str(eps.index(outcome)): 'outcome_age'}
    for i in remained_eps:
        ep_index_old = str(eps.index(i))
        ep_index_new = str(remained_eps.index(i))
        if 'mo_ep' + ep_index_old in dataset.columns:
            col_dict['mo_ep' + ep_index_old] = 'mo_ep' + ep_index_new
        if 'mo_age' + ep_index_old in dataset.columns:
            col_dict['mo_age' + ep_index_old] = 'mo_age' + ep_index_new
        if 'fa_ep' + ep_index_old in dataset.columns:
            col_dict['fa_ep' + ep_index_old] = 'fa_ep' + ep_index_new
        if 'fa_age' + ep_index_old in dataset.columns:
            col_dict['fa_age' + ep_index_old] = 'fa_age' + ep_index_new
    dataset = dataset.rename(columns=col_dict)
    return dataset, remained_eps


def test_match_quality(dataset, outcome):
    ep_index = str(eps.index(outcome))
    for i in ['sex', 'ch_year', 'fa_year', 'mo_year', 'number_of_sib', 'province']:
        print(i+': ', get_cohend(dataset[dataset['ch_ep'+ep_index] == 1][i],
                                 dataset[dataset['ch_ep'+ep_index] == 0][i]))
    print('---------------------------------------')
    for i in ['ses', 'job', 'edulevel', 'edufield', 'number_of_children', 'in_social_assistance_registries',
              'in_vaccination_registry', 'in_infect_dis_registry', 'in_malformations_registry',
              'in_cancer_registry', 'ever_married', 'lang']:
        try:
            print(i + ': ', get_cohend(dataset[dataset['ch_ep' + ep_index] == 1][i],
                                       dataset[dataset['ch_ep' + ep_index] == 0][i]))
        except ValueError:
            print(i)


print('Start to load data...')
# load data
study_population = pd.read_csv('df.csv')
print('Start to match...')
data = case_control_matching(study_population, OUTCOME, MATCH_FACTORS, MATCH_NUM)

print('Test the quality of the data...')
# reasonable matching
# draw_distribution(data, 'edulevel', OUTCOME, 'Case-control')
# draw_distribution(study_population, 'edulevel', OUTCOME, 'Study population')
test_match_quality(data, OUTCOME)

print('Start to rename the columns and then save the data...')
# only diseases whose n_cases > 20
data, eps_remain = remove_unnecessary_columns(data, OUTCOME)
# add number of ADs individuals have to the data
df_ch = study_population[['ID']+[i for i in study_population.columns if re.findall(r'ch_ep\d+', i)
                                 and (eps[int(re.findall(r'_ep(\d+)', i)[0])] in eps_remain)]]
data_ch = data[['ID']]
data_ch = data_ch.merge(df_ch, 'left', on='ID')
data['num_of_eps_ch'] = data_ch.sum(axis=1)
# add number of ADs individuals' parents have to the data
eps_mo = [i for i in data.columns if re.findall(r'mo_ep\d+', i)]
eps_fa = [i for i in data.columns if re.findall(r'fa_ep\d+', i)]
data['num_of_eps_fa'] = data[eps_fa].sum(axis=1)
data['num_of_eps_mo'] = data[eps_mo].sum(axis=1)
# save the data
data.to_csv('data_'+OUTCOME+'.csv', index=None)
with open('eps_'+OUTCOME+'.json', 'w') as f:
    json.dump(eps_remain, f)
print('Done!')
