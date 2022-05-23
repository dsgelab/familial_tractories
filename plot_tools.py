import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from basic_tools import eps


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
    df_mother = results[results.who == 'mother']
    df_father = results[results.who == 'father']
    # Create a color palette: https://www.w3schools.com/colors/colors_picker.asp
    palette = dict(zip(['Father', 'Mother', 'Child'], ['#35a6d4', '#ffa0aa', '#73b400']))
    plt.figure(figsize=(10, len(eps) / 3))
    plt.grid()
    for lower, upper, ep in zip(df_mother.hr_025, df_mother.hr_975, df_mother.endpoint):
        if ep not in ep_remove_mo:
            i = eps.index(ep)
            plt.plot((lower + upper)/2, i - group_delta, 'o', color=palette['Mother'])
            plt.plot((lower, upper), (i - group_delta, i - group_delta), color=palette['Mother'])
            plt.plot([lower, lower], [i - group_delta - bar_delta, i - group_delta + bar_delta], color=palette['Mother'])
            plt.plot([upper, upper], [i - group_delta - bar_delta, i - group_delta + bar_delta], color=palette['Mother'])

    for lower, upper, ep in zip(df_father.hr_025, df_father.hr_975, df_father.endpoint):
        if ep not in ep_remove_fa:
            i = eps.index(ep)
            plt.plot((lower + upper)/2, i, 'o', color=palette['Father'])
            plt.plot((lower, upper), (i, i), color=palette['Father'])
            plt.plot([lower, lower], [i - bar_delta, i + bar_delta], color=palette['Father'])
            plt.plot([upper, upper], [i - bar_delta, i + bar_delta], color=palette['Father'])
    plt.yticks(range(len(eps)), eps)
    plt.xlabel('Odds ratio for ' + eps[ep_index] + ' diagnosis', size=14)
    plt.axvline(x=1.0, color='black', linestyle='--')
    # Create legend handles manually
    handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    # Create legend
    plt.legend(handles=handles)
    plt.grid()
    plt.show()