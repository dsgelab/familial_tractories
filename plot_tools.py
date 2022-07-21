import re
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Create a color palette: https://www.w3schools.com/colors/colors_picker.asp
# https://matplotlib.org/2.0.2/examples/color/named_colors.html
palette = dict(zip(['Father', 'Mother', 'Child'], ['blue', 'tomato', 'black']))


def get_ep_name(col_name, eps):
    """
    :param col_name: a string of column name to find out disease name
    :param eps: a list of diseases
    :return: a string of disease name or col_name if col_name is irrelevant to disease
    """
    if re.match(r'^\w{2}_ep\d+$', col_name):
        ep = eps[int(re.findall(r'_ep(\d+)$', col_name)[0])]
        if col_name.startswith('ch_'):
            return 'children with '+ep
        elif col_name.startswith('fa_'):
            return 'father with '+ep
        elif col_name.startswith('mo_'):
            return 'mother with '+ep
        else:
            return col_name
    else:
        return col_name


def draw_distribution(dataset, main_col, sub_col, title=''):
    """
    :param dataset: a DataFrame of study population
    :param main_col: a string of column name for the variable to display on x-axis
    :param sub_col: a string of column name for the variable to be categorized in legend
    :param title: a string of suptitle in case there is an extra explanation of the plot
    :return: a distribution plot of main_col by sub_col
    """
    summary = pd.crosstab(dataset[main_col], dataset[sub_col])
    summary = (100. * summary / summary.sum()).round(1)
    summary.plot(kind='bar', figsize=(10, 6))
    full_title = 'Distribution of '+get_ep_name(main_col)+' by '+get_ep_name(sub_col)
    plt.title(full_title, size=20)
    plt.suptitle(title, size=16)
    plt.xlabel(get_ep_name(main_col), size=14)
    plt.ylabel('Percentage (%)', size=14)
    plt.legend(dataset[sub_col].unique(), prop={'size':14}, title=get_ep_name(sub_col), title_fontsize=14)
    plt.xticks(rotation=0, size=14)
    plt.yticks(rotation=0, size=14)
    plt.tight_layout()
    plt.show()


def draw_grouped_bar_plot(dataset, col_group, col_num, title, ylabel='Number of individuals'):
    """
    :param dataset: a DataFrame of summary statistics
    :param col_group: a string of column name for grouping
    :param col_num: a string of column name for numbers to display
    :param title: a string of title in plot
    :param ylabel: a string of title on y axis
    :return: a distribution plot of col_num by col_group
    """
    dataset = dataset[['endpoint',col_group,col_num]]
    dataset.pivot(index='endpoint', columns=[col_group], values=[col_num]).plot.bar(figsize=(20,6))
    plt.legend(title='')
    plt.ylabel(ylabel, size=14)
    plt.title(title, size=20)
    plt.show()


def plot_odds_ratio(dataset, outcome):
    """
    :param dataset: a DataFrame of summary statistics
    :param outcome: a string which indicates the outcome disease name
    :return: a odds ratio plot of all the diseases in the list
    """
    dataset = dataset.sort_values(by='endpoint')
    dataset.index = range(len(dataset))
    plt.figure(figsize=(15, 5))
    plt.box(False)
    plt.grid()
    eps_sig = []
    for i, row in dataset.iterrows():
        alpha = 1 if row.pval <= 0.05 / len(dataset) else 0.08
        if alpha == 1:
            eps_sig.append(row.endpoint)
        plt.plot((i, i), (row.hr_025, row.hr_975), color='green', alpha=alpha)
        plt.plot(i, (row.hr_025+row.hr_975)/2, 's', color='green', alpha=alpha)

    plt.xticks(range(len(dataset)), dataset.endpoint.tolist(), rotation=90)
    plt.ylabel('Odds ratio', size=12)
    plt.axhline(y=1.0, color='black', linestyle='--', linewidth=1)
    plt.grid()
    plt.title(outcome, size=20)
    plt.show()
    return eps_sig


# "endpoint","who","se","pval","hr","hr_025","hr_975","note"
def process_crossed_data(dataset, note):
    """
    :param dataset: a DataFrame of summary statistics
    :param note: 'boy', 'girl'
    :return: a DataFrame of processed summary statistics
    """
    df1 = dataset[(dataset.note == note) & (dataset.who == 'father')][['endpoint',"pval", "hr_025", "hr_975", "se"]]
    df1 = df1.rename(columns={"pval": 'p_fa', "hr_025": 'lower_fa', "hr_975": 'upper_fa', "se": 'se_fa'})
    df2 = dataset[(dataset.note == note) & (dataset.who == 'mother')][['endpoint',"pval", "hr_025", "hr_975", "se"]]
    df2 = df2.rename(columns={"pval": 'p_mo', "hr_025": 'lower_mo', "hr_975": 'upper_mo', "se": 'se_mo'})
    df1 = df1.merge(df2, 'outer', on='endpoint')

    df1['hr_fa'] = (df1.lower_fa + df1.upper_fa) / 2
    df1['hr_mo'] = (df1.lower_mo + df1.upper_mo) / 2
    df1['hr_test'] = (df1.hr_fa - df1.hr_mo) / np.sqrt(df1.se_fa ** 2 + df1.se_mo ** 2)  # t-test
    df1['hr_p'] = 2 * scipy.stats.norm.cdf(-np.abs(df1.hr_test))
    df1['hr_significant'] = [True if i < 0.05 / len(df1) else False for i in df1.hr_p]

    return df1


def plot_crossed_odds_ratio(data, note, color):
    """
    :param data: a DataFrame of summary statistics
    :param note: a string which indicates the note content, boy or girl
    :param color: a string which indicates the color to display
    :return: a odds ratio plot of the diseases by groups
    """
    dataset = process_crossed_data(data, note)
    plt.figure(figsize=(8, 8))
    plt.box(False)

    for _, row in dataset.iterrows():
        x, y = row.hr_fa, row.hr_mo
        if row.hr_significant:
            alpha = 1
            plt.annotate(row.endpoint, (x*1.05, y*1.05))
        else:
            alpha = .1
        plt.plot(x, y, 's', color=color, alpha=alpha)
        plt.plot((row.lower_fa, row.upper_fa), (y, y), color=color, alpha=alpha)
        plt.plot((x, x), (row.lower_mo, row.upper_mo), color=color, alpha=alpha)

    plt.title(note, size=16)
    plt.xlabel('Father', size=12)
    plt.ylabel('Mother', size=12)
    plt.xscale('log')
    plt.yscale('log')
    ticks = [.5, 1, 3, 12]
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.plot([ticks[0], ticks[-1]], [ticks[0], ticks[-1]], color='grey', linestyle='--', linewidth=1)
    ax = plt.gca()
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.get_yaxis().set_major_formatter(ScalarFormatter())
    plt.grid()
    plt.show()


def plot_comparisons(dataset):
    """
    :param dataset: a DataFrame of summary statistics
    :return: a odds ratio plot of all the diseases in the list
    """
    dataset = dataset.sort_values(by='endpoint')
    dataset.index = range(len(dataset))
    # set up the figure: frameon=False removes frames
    fig, (ax2, ax1) = plt.subplots(nrows=2, sharex=True, figsize=(15, 5), subplot_kw=dict(frameon=False))
    plt.subplots_adjust(hspace=.5)
    for i, row in dataset.iterrows():
        alpha_rg = 1 if row.p_rg <= 0.05 / len(dataset) else 0.08
        alpha_hr = 1 if row.p_hr <= 0.05 / len(dataset) else 0.08
        ax1.plot((i, i), (row.rg_025, row.rg_975), color='tomato', alpha=alpha_rg)
        ax1.plot(i, (row.rg_025+row.rg_975)/2, 's', color='tomato', alpha=alpha_rg)
        ax2.plot((i, i), (row.hr_025, row.hr_975), color='green', alpha=alpha_hr)
        ax2.plot(i, (row.hr_025+row.hr_975)/2, 's', color='green', alpha=alpha_hr)
    ax1.axhline(y=0.0, color='black', linestyle='--', linewidth=1)
    ax2.axhline(y=1.0, color='black', linestyle='--', linewidth=1)
    plt.xticks(range(len(dataset)), dataset.endpoint.tolist(), rotation=90)
    ax1.set_ylabel('Genetic correlation', size=12)
    ax2.set_ylabel('Registry-based odds ratio', size=12)
    plt.show()
