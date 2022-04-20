import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def plot_odds_ratio(table, demos, eps, ep_index, ep_remove, group_delta=.1, bar_delta=.1):
    # Create a color palette: https://www.w3schools.com/colors/colors_picker.asp
    palette = dict(zip(['Father', 'Mother', 'Child'], ['#35a6d4', '#ffa0aa', '#73b400']))
    plt.figure(figsize=(10, len(eps)/3))
    plt.grid()
    mo_df = table.loc[table.index.str.startswith('mo_'), ['ci_1', 'ci_2']]
    fa_df = table.loc[table.index.str.startswith('fa_'), ['ci_1', 'ci_2']]
    ch_df = table.loc[table.index.str.startswith('ch_'), ['ci_1', 'ci_2']]
    ch_demo_df = table.loc[
        ~(table.index.str.startswith('mo_') | table.index.str.startswith('fa_') | table.index.str.startswith('ch_') | (table.index == 'const')),
        ['ci_1', 'ci_2']]
    parent_list = [i for i in range(len(eps)) if i not in ep_remove]
    for lower, upper, i in zip(np.exp(mo_df.ci_1), np.exp(mo_df.ci_2), parent_list):
        plt.plot((lower + upper) / 2, i - group_delta, 'o', color=palette['Mother'])
        plt.plot((lower, upper), (i - group_delta, i - group_delta), color=palette['Mother'])
        plt.plot([lower, lower], [i - group_delta - bar_delta, i - group_delta + bar_delta], color=palette['Mother'])
        plt.plot([upper, upper], [i - group_delta - bar_delta, i - group_delta + bar_delta], color=palette['Mother'])
    for lower, upper, i in zip(np.exp(fa_df.ci_1), np.exp(fa_df.ci_2), parent_list):
        plt.plot((lower + upper) / 2, i, 'o', color=palette['Father'])
        plt.plot((lower, upper), (i, i), color=palette['Father'])
        plt.plot([lower, lower], [i - bar_delta, i + bar_delta], color=palette['Father'])
        plt.plot([upper, upper], [i - bar_delta, i + bar_delta], color=palette['Father'])
    ch_list = parent_list
    ch_list.remove(ep_index)
    for lower, upper, i in zip(np.exp(ch_df.ci_1), np.exp(ch_df.ci_2), ch_list):
        plt.plot((lower + upper) / 2, i + group_delta, 'o', color=palette['Child'])
        plt.plot((lower, upper), (i + group_delta, i + group_delta), color=palette['Child'])
        plt.plot([lower, lower], [i + group_delta - bar_delta, i + group_delta + bar_delta], color=palette['Child'])
        plt.plot([upper, upper], [i + group_delta - bar_delta, i + group_delta + bar_delta], color=palette['Child'])
    for lower, upper, i in zip(np.exp(ch_demo_df.ci_1), np.exp(ch_demo_df.ci_2), range(len(eps), len(ch_demo_df) + len(eps))):
        plt.plot((lower + upper) / 2, i, 'o', color=palette['Child'])
        plt.plot((lower, upper), (i, i), color=palette['Child'])
        plt.plot([lower, lower], [i - bar_delta, i + bar_delta], color=palette['Child'])
        plt.plot([upper, upper], [i - bar_delta, i + bar_delta], color=palette['Child'])
    plt.yticks(range(len(eps + demos)), eps + demos)
    plt.xlabel('Odds ratio for ' + eps[ep_index] + ' diagnosis', size=14)
    plt.axvline(x=1.0, color='black', linestyle='--')
    # Create legend handles manually
    handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    # Create legend
    plt.legend(handles=handles)
    plt.grid()
    plt.show()
