import re
import tqdm
import pandas as pd
import numpy as np
import statsmodels.api as sm
from plot_tools import plot_odds_ratio

# https://github.com/FINNGEN/CS-PRS-pipeline
T1D = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_phecode-250.1-both_sexes.tsv.sscore'
M13_RHEUMA = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_RA_GWASmeta_European_v2.txt.sscore'
G6_MS = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_discovery_metav3.0.meta.sscore'
L12_PSORIASIS = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_phecode-696.4-both_sexes.tsv.sscore'
K11_IBD = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_EUR.IBD.gwas_info03_filtered.assoc.sscore'
E4_HYTHY_AI_STRICT = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_UKBB_HYPO.txt.sscore'
T2D = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_Mahajan.NatGenet2018b.T2D.European.txt.sscore'
# GEST_DIABETES
K11_CROHN = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_EUR.CD.gwas_info03_filtered.assoc.sscore'
J10_ASTHMA = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_gabriel_asthma_meta-analysis_36studies_format_repository_NEJM.txt.sscore'

eps_dict = {
    'T1D': T1D, 'M13_RHEUMA': M13_RHEUMA, 'G6_MS': G6_MS, 'E4_HYTHY_AI_STRICT': E4_HYTHY_AI_STRICT,
    'K11_IBD': K11_IBD, 'L12_PSORIASIS': L12_PSORIASIS, 'T2D': T2D, 'K11_CROHN': K11_CROHN, 'J10_ASTHMA': J10_ASTHMA
}

eps_with_prs = ['T1D_STRICT', 'M13_RHEUMA', 'G6_MS', 'E4_HYTHY_AI_STRICT',
                'K11_IBD', 'L12_PSORIASIS', 'T2D', 'K11_CROHN', 'J10_ASTHMA']

eps = ['T1D_STRICT', 'M13_RHEUMA', 'M13_SJOGREN', 'M13_SYSTSLCE', 'M13_DERMATOPOLY', 'M13_WEGENER', 'D3_ALLERGPURPURA',
       'M13_MCTD', 'SLE_FG', 'I9_RHEUFEV', 'G6_MS', 'G6_DISSOTH', 'AUTOIMMUNE_HYPERTHYROIDISM', 'E4_THYROIDITAUTOIM',
       'E4_HYTHY_AI_STRICT', 'E4_GRAVES_OPHT_STRICT', 'E4_ADDISON', 'D3_AIHA_OTHER', 'D3_ITP', 'D3_ANAEMIA_B12_DEF',
       'K11_COELIAC', 'K11_IBD', 'G6_MYASTHENIA', 'G6_OTHDEMYEL', 'G6_GUILBAR', 'H7_IRIDOCYC_ANTER', 'CHIRBIL_PRIM',
       'L12_PSORIASIS', 'L12_VITILIGO', 'L12_ALOPECAREATA', 'L12_PEMPHIGOID', 'L12_DERMATHERP', 'N14_IGA_NEPHROPATHY',
       'T2D', 'GEST_DIABETES', 'K11_CROHN', 'J10_ASTHMA']
# eps <- fromJSON(file=paste0(OUTCOME,"/eps_",OUTCOME,".json")) # as same as eps above if OUTCOME = 'T1D_STRICT'
who_dict = {'child': 'ch', 'mother': 'mo', 'father': 'fa'}

fam_path = '/finngen/library-red/finngen_R9/kinship_1.0/data/finngen_R9_pedigree.fam'
event_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_endpoint_longitudinal_1.0.txt.gz'

# load family pedigree
fam = pd.read_csv(fam_path, sep='\t', header=None)
fam.columns = ['family_id', 'child_id', 'father_id', 'mother_id', 'sex', 'phenotype']
fam = fam[(fam.father_id.str.startswith('FG')) & (fam.mother_id.str.startswith('FG'))]
# load confounding vars
confounding = pd.read_csv('confounding.csv')
confounding = confounding[confounding.PC1 != 0]
confounding = fam[['child_id']].merge(confounding.rename(columns={'finngen_id': 'child_id'}), 'left', on='child_id')
confounding = confounding[~confounding.age.isna()]
confounding = confounding.sort_values(by='child_id')
# clean fam data
fam = confounding[['child_id']].merge(fam, 'left', on='child_id')
fam = fam.sort_values(by='child_id')
confounding.index = range(len(confounding))
fam.index = range(len(fam))

events = pd.read_csv(event_path, sep='\t')
events_sub = events[events.ENDPOINT == eps[0]]
fam['outcome_bool'] = np.select([
    (fam['child_id'].isin(events_sub.FINNGENID)), (~fam['child_id'].isin(events_sub.FINNGENID))
], [1, 0])


for ep in tqdm.tqdm(eps):
    if ep in eps_with_prs:
        ss_ep = pd.read_csv(eps_dict[ep], sep='\t')
        num = eps.index(ep)
        fam = fam.merge(ss_ep[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'mother_id',
                        'SCORE1_AVG': 'mo_prs' + str(num)}), 'left', 'mother_id')
        fam = fam.merge(ss_ep[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'father_id',
                        'SCORE1_AVG': 'fa_prs' + str(num)}), 'left', 'father_id')
        if ep == eps[0]:
            fam = fam.merge(ss_ep[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'child_id',
                            'SCORE1_AVG': 'outcome_prs'}), 'left', 'child_id')

for i in fam.columns:
    if re.match(r'^\w+_prs\d{0,3}$', i):
        fam[i+'_scaled'] = (fam[i] - np.mean(fam[i]))/np.sqrt(np.var(fam[i]))


def regression(index_num, who, model_name):
    col_name = who_dict[who] + '_prs' + str(index_num) + '_scaled'
    x = pd.concat([fam[[col_name]], confounding.iloc[:, 1:]], axis=1)
    x = sm.add_constant(x)
    if model_name == 'logit':
        y = fam['outcome_bool']
        model = sm.Logit(y, x).fit(disp=0)
        res_df = model.summary2().tables[1]
        pval = res_df.loc[col_name, 'P>|z|']
    else:
        y = fam['outcome_prs_scaled']
        model = sm.OLS(y, x).fit(disp=0)
        res_df = model.summary2().tables[1]
        pval = res_df.loc[col_name, 'P>|t|']
    hr_025 = np.exp(res_df.loc[col_name, '[0.025'])
    hr_975 = np.exp(res_df.loc[col_name, '0.975]'])
    se = res_df.loc[col_name, 'Std.Err.']
    coef = res_df.loc[col_name, 'Coef.']
    return [coef, who, se, pval, hr_025, hr_975]


def model_loop(note, res_df, model_name):
    for ep in tqdm.tqdm(eps):
        if ep in eps_with_prs:
            index_num = eps.index(ep)
            res_fa = regression(index_num, 'father', model_name)
            res_mo = regression(index_num, 'mother', model_name)
            res_df = res_df.append(pd.Series([ep] + res_fa + [note], index=res_df.columns), ignore_index=True)
            res_df = res_df.append(pd.Series([ep] + res_mo + [note], index=res_df.columns), ignore_index=True)
    return res_df


res = pd.DataFrame(columns=['endpoint', 'coef', 'who', 'se', 'pval', 'hr_025', 'hr_975', 'note'])
# model_name can be regression or logit
res = model_loop('all', res, 'regression')

plot_odds_ratio(res, eps, eps[0])

