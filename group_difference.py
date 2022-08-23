import json
import tqdm
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

OUTCOME = 'T1D_STRICT'

snp_strings = '''11:2160994:A:T, 1:113834946:A:G, 11:2140628:A:G, 11:2173684:A:G, 6:34649486:T:C, 
6:35281810:C:T, 6:34969284:C:T, 2:203880280:T:C, 6:35386132:G:T, 6:36095300:G:A, 21:42421796:G:A, 
12:56007301:G:A, 14:100840110:T:C, 7:5397122:C:T, 6:36154795:C:T, 19:48703417:G:A, 21:44294411:C:T, 
14:100835645:A:C, 10:6068912:C:A, 22:29789624:C:T, 6:34968003:A:C, 10:88291560:A:G, 10:6052734:C:T, 
12:8901640:G:T, 6:21919156:G:A, 16:11072145:G:A, 6:34658871:G:A, 16:75200974:C:A, 12:111446804:T:C, 
6:90267049:G:A, 16:28494339:G:C, 3:32294869:A:G, 6:23202837:T:C, 6:22219423:A:G, 9:114331148:C:T, 
12:55976939:C:T, 6:34612214:GACCAGCCTATGGCAAC:G, 6:35225659:G:A'''
snp_strings = ', '+snp_strings
snp_strings = snp_strings.replace(', ', ', chr')
snp_strings = snp_strings.replace(':', '_')
snp_strings = snp_strings[2:]

# obtain the individual data directly in notebook
# !plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 --snps $snp_strings --recode --out selected_snp
ped_head = {
    0: 'family_id',
    1: 'finngen_id',
    2: 'father_id',
    3: 'mother_id',
    4: 'sex',
    5: 'phenotype'
}
# In .ped file, the first 6 cols are listed as ped_head
# and then each two cols in the rest are a pair of alleles tied to a specific SNP
ped = pd.read_csv('selected_snp.ped', sep=' ', header=None).rename(columns=ped_head)
# get the right order of the selected SNPs
snp_ref = pd.read_csv('selected_snp.map', sep='\t', header=None)[1].tolist()
# add dosages to the dataframe
for i in tqdm.tqdm(range(len(snp_ref))):
    allele1 = np.select([(ped[6 + 2 * i] == snp_ref[i].split('_')[-1]),
                         (ped[6 + 2 * i] != snp_ref[i].split('_')[-1])], [1, 0])
    allele2 = np.select([(ped[6 + 2 * i + 1] == snp_ref[i].split('_')[-1]),
                         (ped[6 + 2 * i + 1] != snp_ref[i].split('_')[-1])], [1, 0])
    dosage = allele1 + allele2
    ped[snp_ref[i]] = dosage

fam_path = '/finngen/library-red/finngen_R9/kinship_1.0/data/finngen_R9_pedigree.fam'
fam = pd.read_csv(fam_path, sep='\t', header=None).rename(columns=ped_head)
fam = fam[(fam.fater_id.str.startswith('FG')) & (fam.mother_id.str.startswith('FG'))] # 10581

# load hla dataset
hla_df = pd.read.csv('hla_df.csv')
hla_df = hla_df.merge(ped[['finngen_id']+snp_ref], 'left', on='finngen_id')
hla_df_fam = hla_df.merge(fam[['finngen_id', 'father_id', 'mother_id']], 'inner', on='finngen_id')

event_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_endpoint_longitudinal_1.0.txt.gz'
events = pd.read_csv(event_path, sep='\t')
eps = json.load(open('eps_'+OUTCOME+'.json', 'r'))
for ep in tqdm.tqdm(eps):
    events_sub = events[events.ENDPOINT == ep]
    hla_df_fam[ep+'_fa'] = np.select([
        (hla_df_fam.father_id.isin(events_sub.finngen_id)), (~hla_df_fam.father_id.isin(events_sub.finngen_id))
    ], [1, 0])
    hla_df_fam[ep + '_mo'] = np.select([
        (hla_df_fam.mother_id.isin(events_sub.finngen_id)), (~hla_df_fam.mother_id.isin(events_sub.finngen_id))
    ], [1, 0])

eps_ads = ['T1D_STRICT', 'M13_RHEUMA', 'SLE_FG', 'AUTOIMMUNE_HYPERTHYROIDISM', 'E4_HYTHY_AI_STRICT',
           'D3_ANAEMIA_B12_DEF', 'K11_COELIAC', 'L12_VITILIGO', 'L12_ALOPECAREATA', 'L12_DERMATHERP']

ch_cols, fa_cols, mo_cols = [i for i in eps_ads[1:]], [i+'_fa' for i in eps_ads], [i+'_mo' for i in eps_ads]
ch_cols_, fa_cols_, mo_cols_ = [i for i in eps[1:]], [i+'_fa' for i in eps], [i+'_mo' for i in eps]
hla_df_fam['ch_ads'] = hla_df_fam[ch_cols].sum(axis=1)
hla_df_fam['fa_ads'] = hla_df_fam[fa_cols].sum(axis=1)
hla_df_fam['mo_ads'] = hla_df_fam[mo_cols].sum(axis=1)
hla_df_fam['pa_ads'] = hla_df_fam[fa_cols + mo_cols].sum(axis=1)
hla_df_fam['pa_eps'] = hla_df_fam[fa_cols_ + mo_cols_].sum(axis=1)
hla_df_fam['trio_ads'] = hla_df_fam[fa_cols + mo_cols + ch_cols].sum(axis=1)

hla_df_fam_0 = hla_df_fam[(hla_df_fam.fa_ads == 0) & (hla_df_fam.mo_ads == 0) & (hla_df_fam.T1D_STRICT == 1)] # 738
hla_df_fam_1 = hla_df_fam[((hla_df_fam.fa_ads > 0) | (hla_df_fam.mo_ads > 0)) & (hla_df_fam.T1D_STRICT == 1)] # 461
hla_df_fam_2 = hla_df_fam[(hla_df_fam.fa_ads > 0) & (hla_df_fam.mo_ads > 0) & (hla_df_fam.T1D_STRICT == 1)] # 50

selected_genes = ['DQB1*03:02', 'DQB1*02:01', 'DRB5*01:01', 'DQB1*04:02', 'DQB1*06:04', 'DRB1*04:01', 'DRB3*01:01',
                  'DQB1*03:01', 'B*39:01', 'A*02:01', 'A*24:02', 'A*03:01', 'DQA1*01:04', 'DRB4*01:03N', 'DPB1*03:01',
                  'DRB1*09:01', 'DRB1*10:01', 'DRB1*13:01']

# chi-squared test with similar proportions
# to check the SNP AF:
# https://gnomad.broadinstitute.org/variant/1-113834946-A-G?dataset=gnomad_r3
# Comparison 1: children whose parents have selected ADs vs children whose parents have no AD
# if child inherits healthy blocks of genes from parents in group 1?
for i in snp_ref+selected_genes:
    s0 = hla_df_fam_0[i].round().value_counts()
    s1 = hla_df_fam_1[i].round().value_counts()
    s2 = hla_df_fam_2[i].round().value_counts()
    table = [
        [dict(s0).get(0.0, 0), dict(s0).get(1.0, 0), dict(s0).get(2.0, 0)],
        [dict(s1).get(0.0, 0)+dict(s2).get(0.0, 0), dict(s1).get(1.0, 0)+dict(s2).get(1.0, 0),
         dict(s1).get(2.0, 0)+dict(s2).get(2.0, 0)]
    ]
    if table[0][2] == 0 and table[1][2] == 0:
        table = [table[0][:2], table[1][:2]]
    if table[0][1] == 0 and table[1][1] == 0:
        pass
    else:
        try:
            res = chi2_contingency(table) # stat, p, dof, expected
            if res[1] < 0.05:
                print(i, ' : ', res[1])
                print(table)
        except ValueError:
            print(i, ' : ')
            print(table)

# Comparison 2: children who have selected ADs vs children without any ADs
hla_df['ch_ads'] = hla_df[ch_cols].sum(axis=1)
hla_df['ch_eps'] = hla_df[ch_cols_].sum(axis=1)
hla_df_0 = hla_df[(hla_df.ch_eps == 0) & (hla_df.T1D_STRICT == 1)] # ch_eps:2099 ch_ads:2605
hla_df_1 = hla_df[(hla_df.ch_ads != 0) & (hla_df.T1D_STRICT == 1)] # ch_ads:795
for i in snp_ref+selected_genes:
    s0 = hla_df_0[i].round().value_counts()
    s1 = hla_df_1[i].round().value_counts()
    table = [
        [dict(s0).get(0.0, 0), dict(s0).get(1.0, 0), dict(s0).get(2.0, 0)],
        [dict(s1).get(0.0, 0), dict(s1).get(1.0, 0), dict(s1).get(2.0, 0)]
    ]
    if table[0][2] == 0 and table[1][2] == 0:
        table = [table[0][:2], table[1][:2]]
    try:
        res = chi2_contingency(table) # stat, p, dof, expected
        if res[1] < 0.05:
            print(i, ' : ', res[1])
            print(table)
    except ValueError:
        print(i, ' : ')
        print(table)

