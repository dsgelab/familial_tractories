import re
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib as mpl

# Create a color palette: https://www.w3schools.com/colors/colors_picker.asp
# https://matplotlib.org/2.0.2/examples/color/named_colors.html
palette = dict(zip(['Father', 'Mother'], ['blue', 'tomato']))


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


def plot_odds_ratio(results, eps, outcome, group_delta=.1, bar_cap=.1):
    """
    :param results: a DataFrame of summary statistics
    :param eps: a list of diseases
    :param outcome: a string which indicates the outcome disease name
    :param group_delta: a float which indicates the distance between father's OR and mother's OR given the same disease
    :param bar_cap: a float which indicates the length of error bar cap
    :return: a odds ratio plot of all the diseases in the list
    """
    df_mother = results[results.who == 'mother']
    df_father = results[results.who == 'father']
    plt.figure(figsize=(10, len(eps) / 3))
    plt.grid()

    for lower, upper, ep, pval in zip(df_mother.hr_025, df_mother.hr_975, df_mother.endpoint, df_mother.pval):
        i = eps.index(ep)
        alpha = 1 if pval <= 0.05 / len(eps) else 0.05 if pval > 0.05 else 0.4
        plt.plot((lower + upper) / 2, i - group_delta, '.', color=palette['Mother'], alpha=alpha)
        plt.plot((lower, upper), (i - group_delta, i - group_delta), color=palette['Mother'], alpha=alpha)
        plt.plot([lower, lower], [i - group_delta - bar_cap, i - group_delta + bar_cap],
                 color=palette['Mother'], alpha=alpha)
        plt.plot([upper, upper], [i - group_delta - bar_cap, i - group_delta + bar_cap],
                 color=palette['Mother'], alpha=alpha)

    for lower, upper, ep, pval in zip(df_father.hr_025, df_father.hr_975, df_father.endpoint, df_father.pval):
        i = eps.index(ep)
        alpha = 1 if pval <= 0.05 / len(eps) else 0.05 if pval > 0.05 else 0.4
        plt.plot((lower + upper) / 2, i + group_delta, '.', color=palette['Father'], alpha=alpha)
        plt.plot((lower, upper), (i + group_delta, i + group_delta), color=palette['Father'], alpha=alpha)
        plt.plot([lower, lower], [i + group_delta - bar_cap, i + group_delta + bar_cap],
                 color=palette['Father'], alpha=alpha)
        plt.plot([upper, upper], [i + group_delta - bar_cap, i + group_delta + bar_cap],
                 color=palette['Father'], alpha=alpha)
    plt.yticks(range(len(eps)), eps)
    plt.xlabel('Odds ratio for ' + outcome + ' diagnosis', size=14)
    plt.axvline(x=1.0, color='grey', linestyle='--')
    # Create legend handles manually
    handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    # Create legend
    plt.legend(handles=handles)
    plt.grid()
    plt.show()


# "endpoint","who","se","pval","hr","hr_025","hr_975","note"
def process_crossed_data(data, note_tuple):
    """
    :param data: a DataFrame of summary statistics
    :param note_tuple: a tuple which indicates the two group names
    :return: two DataFrames of processed summary statistics
    """
    data_ = data[data.pval < 0.05/len(data.endpoint.unique())]
    data = data[data.endpoint.isin(data_.endpoint.unique().tolist())]
    res1, res2 = data[data.note == note_tuple[0]], data[data.note == note_tuple[1]]

    def process(who):
        """
        :param who: 'Mother', 'Father'
        :return: a DataFrame of processed summary statistics
        """
        df1 = res1[(res1.who == who) & (~res1.se.isna())][["endpoint", "pval", "hr_025", "hr_975", "se"]]
        df1 = df1.rename(columns={"pval": 'p1', "hr_025": 'lower1', "hr_975": 'upper1', "se": 'se1'})
        df2 = res2[(res2.who == who) & (~res2.se.isna())][["endpoint", "pval", "hr_025", "hr_975", "se"]]
        df2 = df2.rename(columns={"pval": 'p2', "hr_025": 'lower2', "hr_975": 'upper2', "se": 'se2'})
        df1 = df1.merge(df2, 'outer', on='endpoint')

        df1['hr1'] = (df1.lower1 + df1.upper1) / 2
        df1['hr2'] = (df1.lower2 + df1.upper2) / 2
        df1['hr_test'] = (df1.hr1 - df1.hr2) / np.sqrt(df1.se1 ** 2 + df1.se2 ** 2)  # t-test
        df1['hr_p'] = 2 * scipy.stats.norm.cdf(-np.abs(df1.hr_test))
        df1['hr_significant'] = [True if i < 0.05 / len(df1) else False for i in df1.hr_p]

        return df1

    return process('mother'), process('father')


def plot_crossed_odds_ratio(data, note_tuple, outcome):
    """
    :param data: a DataFrame of summary statistics
    :param note_tuple: a tuple which indicates the two group names
    :param outcome: a string which indicates the outcome disease name
    :return: a odds ratio plot of the diseases by groups
    """
    df1, df2 = process_crossed_data(data, note_tuple)
    plt.figure(figsize=(12, 12))

    def draw_cross(dataset, who):
        """
        :param dataset: a string of column name to find out disease name
        :param who: 'Mother', 'Father'
        :return: a string of disease name or col_name if col_name is irrelevant to disease
        """
        for _, row in dataset.iterrows():
            x, y = (row.lower1 + row.upper1) / 2, (row.lower2 + row.upper2) / 2
            if row.hr_significant:
                alpha = 1
                plt.plot(x, y, '.', color=palette[who], alpha=alpha)
                plt.plot((row.lower1, row.upper1), (y, y), color=palette[who], alpha=alpha)
                plt.plot((x, x), (row.lower2, row.upper2), color=palette[who], alpha=alpha)
                plt.annotate(row.endpoint, (x, y))
            else:
                alpha = .1
                plt.plot(x, y, '.', color=palette[who], alpha=alpha)

    draw_cross(df1, 'Mother')
    draw_cross(df2, 'Father')

    plt.title('Odds ratios for ' + outcome + ' diagnosis', size=20)
    plt.xlabel(note_tuple[0], size=14)
    plt.ylabel(note_tuple[1], size=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks([.35, 1, 5, 10, 20, 30])
    plt.yticks([.35, 1, 5, 10, 20, 30])
    plt.plot([.35, 30], [.35, 30], color='grey', linestyle='--')
    ax = plt.gca()
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.get_yaxis().set_major_formatter(ScalarFormatter())

    #     plt.axline([0, 0], slope=1, color='black', linestyle='--')
    # Create legend handles manually
    handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    # Create legend
    plt.legend(handles=handles)
    plt.grid()
    plt.show()
