"""
METHOD 1: exact matching
@Author: feiyiwang
"""

import tqdm
import json
import itertools
from functools import reduce
import pandas as pd
from basic_tools import eps
from stat_tools import get_cohend, sum_cases

# define matching ratio and matching factors
OUTCOME = 'T1D_STRICT'
MATCH_NUM = 3
MATCH_FACTORS = ['sex', 'ch_year_range', 'fa_year_range', 'mo_year_range', 'sib_number', 'province']
MATCH_FACTORS_IN_MODEL = ['sex', 'ch_year', 'fa_year', 'mo_year', 'number_of_sib', 'province']


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
                match_permutation_keep.append({'conditions': i, 'n_cases': len(match)})
        except Exception as e:
            print(e, i)
    print(len(match_permutation_keep), 'possible permutations were found.')

    # filter out all the permutations that cannot find enough controls
    match_failed = []  # >>66% -> 66% -> 39%
    control_ids = []
    for i in tqdm.tqdm(match_permutation_keep):
        try:
            conditions = [(control_for_match[matching_factors[j]] == i['conditions'][j]) for j in range(len(matching_factors))]
            potential_control = control_for_match[reduce(lambda x, y: x & y, conditions)]
            if len(potential_control) < matching_number * i['n_cases']:
                match_failed.append(i)
            else:
                potential_control = potential_control.sample(n=matching_number * i['n_cases'])
                control_ids += potential_control.ID.tolist()
        except:
            match_failed.append(i)
    print(len(match_failed), 'permutations failed to match to enough controls.')

    # remove all the failed permutations from cases
    cases_to_remove = []
    for i in tqdm.tqdm(match_failed):
        try:
            conditions = [(cases[matching_factors[j]] == i['conditions'][j]) for j in range(len(matching_factors))]
            potential_not_cases = cases[reduce(lambda x, y: x & y, conditions)]
            for j in potential_not_cases.ID:
                cases_to_remove.append(j)
        except:
            print(i)
    print(len(cases_to_remove), 'permutations that could not be used to match enough controls were removed from cases.')
    print('This occupies', '{:.2%}'.format(len(cases_to_remove) / len(cases)), 'of the cases.')

    # final cases and controls
    cases = cases[~cases.ID.isin(cases_to_remove)]
    controls = control_for_match[control_for_match.ID.isin(control_ids)]
    print('Finally, we have', len(cases), 'cases and', len(controls), 'controls.')

    return cases, controls


def remove_unnecessary_columns(dataset, outcome, threshold=20):
    n_cases_mo, ep_remove_mo = sum_cases(dataset, 'mother', threshold)
    n_cases_fa, ep_remove_fa = sum_cases(dataset, 'father', threshold)
    eps_remain = list(set(eps).difference(set(ep_remove_mo).intersection(set(ep_remove_fa))))
    eps_remain = [i for i in eps if i in eps_remain]
    dataset = dataset.drop(columns=['fa_ep' + str(eps.index(i)) for i in ep_remove_fa] +
                                  ['mo_ep' + str(eps.index(i)) for i in ep_remove_mo] +
                                  ['ch_ep' + str(eps.index(i)) for i in eps if i != outcome] +
                                  ['fa_age' + str(eps.index(i)) for i in ep_remove_fa] +
                                  ['mo_age' + str(eps.index(i)) for i in ep_remove_mo] +
                                  ['ch_age' + str(eps.index(i)) for i in eps if i != outcome])
    col_dict = {'ch_ep' + str(eps.index(outcome)): 'outcome', 'ch_age' + str(eps.index(outcome)): 'outcome_age'}
    for i in eps_remain:
        ep_index_old = str(eps.index(i))
        ep_index_new = str(eps_remain.index(i))
        if 'mo_ep' + ep_index_old in dataset.columns:
            col_dict['mo_ep' + ep_index_old] = 'mo_ep' + ep_index_new
        if 'mo_age' + ep_index_old in dataset.columns:
            col_dict['mo_age' + ep_index_old] = 'mo_age' + ep_index_new
        if 'fa_ep' + ep_index_old in dataset.columns:
            col_dict['fa_ep' + ep_index_old] = 'fa_ep' + ep_index_new
        if 'fa_age' + ep_index_old in dataset.columns:
            col_dict['fa_age' + ep_index_old] = 'fa_age' + ep_index_new
    dataset = dataset.rename(columns=col_dict)
    return dataset, eps_remain


def test_match_quality(dataset, outcome, matching_factors_in_model):
    ep_index = str(eps.index(outcome))
    for i in matching_factors_in_model:
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
case, control = case_control_matching(study_population, OUTCOME, MATCH_FACTORS, MATCH_NUM)

print('Start to merge cases and controls...')
# merge the data
data = pd.concat([case, control], axis=0, ignore_index=True)
data = data.sample(frac=1)
# add group id to the data for conditional regression
data['subclass'] = (data.groupby(MATCH_FACTORS).cumcount() == 0).astype(int).cumsum()

print('Test the quality of the data...')
# reasonable matching
test_match_quality(data, OUTCOME, MATCH_FACTORS_IN_MODEL)

print('Start to rename the columns and then save the data...')
# only diseases whose n_cases > 20
data, eps_remain = remove_unnecessary_columns(data, OUTCOME)
# save the data
data.to_csv('data_'+OUTCOME+'.csv', index=None)
with open('eps_'+OUTCOME+'.json', 'w') as f:
    json.dump(eps_remain, f)
print('Done!')
