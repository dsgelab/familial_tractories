import re
import pandas as pd
import numpy as np
import statsmodels.api as sm

# https://github.com/FINNGEN/CS-PRS-pipeline
T1D = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_phecode-250.1-both_sexes.tsv.sscore'
M13_RHEUMA = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_RA_GWASmeta_European_v2.txt.sscore'
G6_MS = 'discovery_metav3.0.meta.gz'
L12_PSORIASIS = 'phecode-696.4-both_sexes.tsv.gz'
K11_IBD = 'EUR.IBD.gwas_info03_filtered.assoc.gz'
HYPOTHYROIDISM = 'UKBB_HYPO.txt.gz'
T2D = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_Mahajan.NatGenet2018b.T2D.European.txt.sscore'
# GEST_DIABETES
K11_CROHN = 'EUR.CD.gwas_info03_filtered.assoc.gz'
J10_ASTHMA = '/finngen/library-red/finngen_R9/prs_1.0/data/finngen_R9_gabriel_asthma_meta-analysis_36studies_format_repository_NEJM.txt.gz'

eps = ['T1D_STRICT', 'M13_RHEUMA', 'M13_SJOGREN', 'M13_SYSTSLCE', 'M13_DERMATOPOLY', 'M13_WEGENER', 'D3_ALLERGPURPURA',
       'M13_MCTD', 'SLE_FG', 'I9_RHEUFEV', 'G6_MS', 'G6_DISSOTH', 'AUTOIMMUNE_HYPERTHYROIDISM', 'E4_THYROIDITAUTOIM',
       'E4_HYTHY_AI_STRICT', 'E4_GRAVES_OPHT_STRICT', 'E4_ADDISON', 'D3_AIHA_OTHER', 'D3_ITP', 'D3_ANAEMIA_B12_DEF',
       'K11_COELIAC', 'K11_IBD', 'G6_MYASTHENIA', 'G6_OTHDEMYEL', 'G6_GUILBAR', 'H7_IRIDOCYC_ANTER',  'CHIRBIL_PRIM',
       'L12_PSORIASIS', 'L12_VITILIGO', 'L12_ALOPECAREATA', 'L12_PEMPHIGOID', 'L12_DERMATHERP', 'N14_IGA_NEPHROPATHY',
       'T2D', 'GEST_DIABETES', 'K11_CROHN', 'J10_ASTHMA']

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

ss_T1D = pd.read_csv(T1D, sep='\t')
fam = fam.merge(ss_T1D[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'child_id', 'SCORE1_AVG': 'outcome_prs'}),
                'left', 'child_id')
fam = fam.merge(ss_T1D[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'mother_id', 'SCORE1_AVG': 'mo_prs0'}),
                'left', 'mother_id')
fam = fam.merge(ss_T1D[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'father_id', 'SCORE1_AVG': 'fa_prs0'}),
                'left', 'father_id')

ss_M13_RHEUMA = pd.read_csv(M13_RHEUMA, sep='\t')
fam = fam.merge(ss_M13_RHEUMA[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'mother_id', 'SCORE1_AVG': 'mo_prs1'}),
                'left', 'mother_id')
fam = fam.merge(ss_M13_RHEUMA[['IID', 'SCORE1_AVG']].rename(columns={'IID': 'father_id', 'SCORE1_AVG': 'fa_prs1'}),
                'left', 'father_id')

# for i in tqdm.tqdm(range(len(eps))):
#
#     fam['ch_ep' + str(i)] = np.select([
#         (fam['ID'].isin(events_sub.ID)), (~df['ID'].isin(events_sub.ID))
#     ], [1, 0])
#     df = df.merge(df_events_sub[['ID', 'AGE']].rename(columns={'AGE': 'ch_age' + str(i)}), 'left', on='ID')

# regression
y = fam['outcome_prs']
x = pd.concat([fam[['fa_prs0', 'mo_prs0']], confounding.iloc[:, 1:]], axis=1)
x = sm.add_constant(x)
model = sm.OLS(y, x).fit(disp=0)
res = model.summary2().tables[1]
np.exp(res.loc['father_prs', '[0.025'])
np.exp(res.loc['father_prs', '0.975]'])
np.exp(res.loc['mother_prs', '[0.025'])
np.exp(res.loc['mother_prs', '0.975]'])

# logit
y = fam['outcome_bool']
model = sm.Logit(y.astype(int), x).fit(disp=0)
res = model.summary2().tables[1]
# This method does not work at all

# logit
y = fam['outcome_bool']
x = pd.concat([fam[['fa_ep0', 'mo_ep0']], confounding.iloc[:, 1:]], axis=1)
model = sm.Logit(y.astype(int), x).fit(disp=0)
res = model.summary2().tables[1]
# This method works not well either
