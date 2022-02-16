import matplotlib.pyplot as plt
import matplotlib as mpl


def plot_ci(table, eps, demos, group_delta=.1, bar_delta=.1):
    # Create a color palette
    palette = dict(zip(['Father', 'Mother', 'Child'], ['#4fcae8', '#ffa0aa', '#73b400']))
    plt.figure(figsize=(10,6))
    plt.grid()
    mo_df = table.loc[table.index.str.startswith('mo_'),['ci_1','ci_2']]
    fa_df = table.loc[table.index.str.startswith('fa_'),['ci_1','ci_2']]
    ch_df = table.loc[~(table.index.str.startswith('mo_')|table.index.str.startswith('fa_')|(table.index=='const')),
                   ['ci_1','ci_2']]
    for lower,upper,i in zip(mo_df.ci_1, mo_df.ci_2, range(len(mo_df))):
        plt.plot((lower+upper)/2,i-group_delta,'o',color=palette['Mother'])
        plt.plot((lower,upper),(i-group_delta,i-group_delta),color=palette['Mother'])
        plt.plot([lower,lower], [i-group_delta-bar_delta,i-group_delta+bar_delta], color=palette['Mother'])
        plt.plot([upper,upper], [i-group_delta-bar_delta,i-group_delta+bar_delta], color=palette['Mother'])
    for lower,upper,i in zip(fa_df.ci_1, fa_df.ci_2, range(len(fa_df))):
        plt.plot((lower+upper)/2,i+group_delta,'o',color=palette['Father'])
        plt.plot((lower,upper),(i+group_delta,i+group_delta),color=palette['Father'])
        plt.plot([lower,lower], [i+group_delta-bar_delta,i+group_delta+bar_delta], color=palette['Father'])
        plt.plot([upper,upper], [i+group_delta-bar_delta,i+group_delta+bar_delta], color=palette['Father'])
    for lower,upper,i in zip(ch_df.ci_1, ch_df.ci_2, range(len(ch_df))):
        plt.plot((lower+upper)/2,i+len(eps),'o',color=palette['Child'])
        plt.plot((lower,upper),(i+len(eps),i+len(eps)),color=palette['Child'])
        plt.plot([lower,lower], [i+len(eps)-bar_delta,i+len(eps)+bar_delta], color=palette['Child'])
        plt.plot([upper,upper], [i+len(eps)-bar_delta,i+len(eps)+bar_delta], color=palette['Child'])
    plt.yticks(range(len(eps+demos)), eps+demos)
    plt.xlabel('Odds ratio for '+eps[0]+' diagnosis', size=14)
    plt.axvline(x=0.0, color='black', linestyle='--')
    # Create legend handles manually
    handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    # Create legend
    plt.legend(handles=handles)
    plt.show()