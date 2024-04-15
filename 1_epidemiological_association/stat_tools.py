import re
import pandas as pd
import numpy as np
from scipy import stats
from basic_tools import eps


# function to calculate Cohen's d for independent samples
def get_cohend(x1, x2):
    """
    :param x1: a Series of one samples
    :param x2: a Series of another samples
    :return: a float of standardized difference
    """
    # calculate the size of samples
    n1, n2 = len(x1), len(x2)
    # calculate the variance of the samples
    s1, s2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    # calculate the pooled standard deviation
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    # calculate the means of the samples
    u1, u2 = np.mean(x1), np.mean(x2)
    # calculate the effect size
    return (u1 - u2) / s


def sum_cases(dataset, who, threshold):
    """
    :param dataset: a DataFrame of study population
    :param who: a string of population group - child, father or mother
    :param threshold: an int as a threshold of minimum n_cases
    :return: a dictionary of n_cases; a list of disease names to remove
    """
    who_dict = {'child': 'ch', 'mother': 'mo', 'father': 'fa'}
    x = dataset[[i for i in dataset.columns if re.match(who_dict[who]+'_ep\d', i)]]
    n_cases = x.sum()
    n_cases_dict = dict(zip(eps, n_cases))
    ep_remove = list(dict(filter(lambda item: item[1] < threshold, n_cases_dict.items())).keys())
    return n_cases_dict, ep_remove


def get_summary(model, X, y):
    """
    :param model: a well-trained model instance
    :param X: a DataFrame of inputs
    :param y: a Series of output
    :return: a DataFrame of key statistics
    """
    coef = np.append(model.intercept_, model.coef_)
    X_matrix = np.append(np.ones((len(X),1)), X, axis=1)
    y_hat = model.predict(X)
    degree_of_freedom = X_matrix.shape[0] - X_matrix.shape[1]
    MSE = (sum((y_hat - y) ** 2)) / degree_of_freedom
    var_coef = MSE * (np.linalg.inv(np.dot(X_matrix.T, X_matrix)).diagonal())
    std_err_coef = np.sqrt(var_coef)
    t_stat_coef = coef / std_err_coef
    p_values_coef = [2 * (1 - stats.t.cdf(np.abs(i), degree_of_freedom)) for i in t_stat_coef]
    t_half_alpha = stats.t.ppf(1 - 0.025, degree_of_freedom)
    ci1_coef = [beta - t_half_alpha * se_beta for beta, se_beta in zip(coef, std_err_coef)]
    ci2_coef = [beta + t_half_alpha * se_beta for beta, se_beta in zip(coef, std_err_coef)]
    summary_df = pd.DataFrame({
        'coef': np.round(coef, 4),
        'std_err': np.round(std_err_coef, 3),
        't_stat': np.round(t_stat_coef, 3),
        'p_value': np.round(p_values_coef, 3),
        'ci_1': np.round(ci1_coef, 3),
        'ci_2': np.round(ci2_coef, 3),
    }, index=['const'] + X.columns.tolist())
    return summary_df