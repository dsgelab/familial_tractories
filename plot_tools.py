import re
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib as mpl
from basic_tools import eps

# Create a color palette: https://www.w3schools.com/colors/colors_picker.asp
# https://matplotlib.org/2.0.2/examples/color/named_colors.html
palette = dict(zip(['Father', 'Mother'], ['navy', 'firebrick']))
palette_light = dict(zip(['Father', 'Mother'], ['royalblue', 'lightcoral']))
palette_bright = dict(zip(['Father', 'Mother'], ['lightblue', 'peachpuff']))


def get_ep_name(col_name):
    """
    :param col_name: a string of column name to find out disease name
    :return: a string of disease name or col_name if col_name is irrelevant to disease
    """
    if re.match('^\w{2}_ep\d+$',col_name):
        ep = eps[int(re.findall('_ep(\d+)$',col_name)[0])]
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
    summary = pd.crosstab(dataset[main_col],dataset[sub_col])
    summary = (100. * summary / summary.sum()).round(1)
    summary.plot(kind='bar',figsize=(10,6))
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


def plot_odds_ratio(results, eps, ep_remove_mo, ep_remove_fa, ep_index=0, group_delta=.3, bar_delta=.1):
    ep_remain_num = len(eps)*2 - len(ep_remove_mo) - len(ep_remove_fa)
    df_mother = results[results.who == 'mother']
    df_father = results[results.who == 'father']

    plt.figure(figsize=(10, len(eps) / 3))
    plt.grid()
    for lower, upper, ep, pval in zip(df_mother.hr_025, df_mother.hr_975, df_mother.endpoint, df_mother.pval):
        if ep not in ep_remove_mo:
            i = eps.index(ep)
            color = palette['Mother'] if pval <= 0.05/ep_remain_num else palette_bright['Mother'] if pval > 0.05 else palette_light['Mother']
            plt.plot((lower + upper)/2, i - group_delta, 'o', color=color)
            plt.plot((lower, upper), (i - group_delta, i - group_delta), color=color)
            plt.plot([lower, lower], [i - group_delta - bar_delta, i - group_delta + bar_delta], color=color)
            plt.plot([upper, upper], [i - group_delta - bar_delta, i - group_delta + bar_delta], color=color)

    for lower, upper, ep, pval in zip(df_father.hr_025, df_father.hr_975, df_father.endpoint, df_father.pval):
        if ep not in ep_remove_fa:
            i = eps.index(ep)
            color = palette['Father'] if pval <= 0.05/ep_remain_num else palette_bright['Father'] if pval > 0.05 else palette_light['Father']
            plt.plot((lower + upper)/2, i, 'o', color=color)
            plt.plot((lower, upper), (i, i), color=color)
            plt.plot([lower, lower], [i - bar_delta, i + bar_delta], color=color)
            plt.plot([upper, upper], [i - bar_delta, i + bar_delta], color=color)
    plt.yticks(range(len(eps)), eps)
    plt.xlabel('Odds ratio for ' + eps[ep_index] + ' diagnosis', size=14)
    plt.axvline(x=1.0, color='black', linestyle='--')
    # Create legend handles manually
    handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    # Create legend
    plt.legend(handles=handles)
    plt.grid()
    plt.show()


# "endpoint","who","se","pval","hr","hr_025","hr_975","note"
def process_crossed_data(data, note_set):

    res1, res2 = data[data.note == note_set[0]], data[data.note == note_set[1]]

    def process(who):
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


def plot_crossed_odds_ratio(data, note_set, name1, name2, ep_index=0):
    df1, df2 = process_crossed_data(data, note_set)
    plt.figure(figsize=(12, 12))

    def draw_cross(data, who):
        """
        :param data: a string of column name to find out disease name
        :param who: 'Mother', 'Father'
        :return: a string of disease name or col_name if col_name is irrelevant to disease
        """
        for _, row in data.iterrows():
            x, y = (row.lower1 + row.upper1) / 2, (row.lower2 + row.upper2) / 2
            if row.hr_significant == True:
                color = palette[who]
                plt.plot(x, y, '.', color=color)
                plt.plot((row.lower1, row.upper1), (y, y), color=color)
                plt.plot((x, x), (row.lower2, row.upper2), color=color)
                plt.annotate(row.endpoint, (x, y))
            else:
                color = palette_bright[who]
                plt.plot(x, y, '.', color=color)

    draw_cross(df1, 'Mother')
    draw_cross(df2, 'Father')

    plt.title('Odds ratios for ' + eps[ep_index] + ' diagnosis', size=20)
    plt.xlabel(name1, size=14)
    plt.ylabel(name2, size=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks([.2,0.5, 1, 5, 10,20, 30])
    plt.yticks([.2,0.5, 1, 5, 10,20, 30])
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