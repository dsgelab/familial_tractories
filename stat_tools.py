import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor


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


# singular matrix
# need to use lasso or ridge to remove collinearity
# np.linalg.det(x_train)
def vif(X):
    """
    :param X: a DataFrame of inputs
    :return: a DataFrame of input names sorted by VIF
    """
    vif_info = pd.DataFrame({
        'vif': [variance_inflation_factor(X.values, i) for i in range(X.shape[1])],
        'col': X.columns
    })
    return vif_info.sort_values('vif', ascending=False)