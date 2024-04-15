import pandas as pd
import numpy as np
import statsmodels.api as sm

OUTCOME = 'T1D_EARLY'
ORIGEN_FACTORS = ['sex', 'ch_year', 'province']
SELECTED_AIDS = ['E4_HYTHY_AI_STRICT', 'K11_COELIAC', 'M13_RHEUMA', 'D3_SARCOIDOSIS', 'D3_ANAEMIA_B12_DEF',
                 'E4_GRAVES_STRICT', 'M13_MCTD', 'M13_SJOGREN', 'L12_VITILIGO', 'CHIRBIL_PRIM', 'D3_AIHA_OTHER',
                 'G6_MYASTHENIA', 'E4_ADDISON', 'L12_PSORIASIS', 'SLE_FG', 'L12_ALOPECAREATA','T1D_STRICT']

END = pd.to_datetime('2020-01-01 00:00:00')
START = pd.to_datetime('1960-01-01 00:00:00')

df_info = pd.read_csv('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv')
df_info.date_of_birth = pd.to_datetime(df_info.date_of_birth)
df = df_info[['FINREGISTRYID','sex','date_of_birth','death_date','id_mother','id_father']]#.date_of_birth
df = df.rename(columns={'FINREGISTRYID': 'ID'})

df_geo = pd.read_csv('/data/projects/project_akarvane/geo/living_province.csv')
df = df.merge(df_geo[['ID', 'maakunta']], 'left', on='ID')

first_event_path = '/data/processed_data/endpointer/R8/densified_first_events_DF8_all_endpoints_2021-09-04.txt'
# first_event_path = '/data/processed_data/endpointer/R10/densified_first_events_DF10_no_omits_2022-09-20.txt'
df_events = pd.read_csv(first_event_path)
df_events = df_events.rename(columns={'FINNGENID': 'ID'})

for ep in ['DEATH',OUTCOME]:
    sub = df_events[df_events.ENDPOINT == ep][['ID','AGE']]
    sub.columns = ['ID',ep+'_onset']
    df = df.merge(sub, 'left')
    df[ep] = np.select([df.ID.isin(sub.ID),(~df.ID.isin(sub.ID))],[1,0])

for ep in SELECTED_AIDS:
    sub = df_events[df_events.ENDPOINT == ep]
    df[ep] = np.select([df.ID.isin(sub.ID),(~df.ID.isin(sub.ID))],[1,0])
    df['fa_'+ep] = np.select([df.id_father.isin(sub.ID),(~df.id_father.isin(sub.ID))],[1,0])
    df['mo_'+ep] = np.select([df.id_mother.isin(sub.ID),(~df.id_mother.isin(sub.ID))],[1,0])

# remove those once stayed/or already live abroad
family = df[~df.ID.isin(df[df.maakunta == 21.0].ID)]

family = family[family.duplicated(subset=['ID']) == False]
family = family[(family.date_of_birth < END) & (df.date_of_birth >= START)]
family = family[(~family.id_father.isna())&(~family.id_mother.isna())]
family['ch_year'] = family['date_of_birth'].dt.year

family['age_2019end'] = pd.to_datetime(END) - family.date_of_birth
family['age_2019end'] = round(family['age_2019end']/np.timedelta64(1,'Y'),2)
family['period1'] = family.age_2019end
family['period2'] = family.T1D_EARLY_onset
family['period3'] = family.DEATH_onset
family['time_to_event'] = family[['period1','period2','period3']].min(axis=1)
family['time_to_event'] = family['time_to_event'].mask(family['time_to_event'] > 20, 20)

family['mo_other_AIDs'] = family[['mo_'+ep for ep in SELECTED_AIDS]].sum(axis=1)
family['fa_other_AIDs'] = family[['fa_'+ep for ep in SELECTED_AIDS]].sum(axis=1)

family['group_t1d'] = np.select([
    (family.mo_T1D_STRICT!=1)&(family.fa_T1D_STRICT!=1),
    (family.mo_T1D_STRICT==1)&(family.fa_T1D_STRICT!=1),
    (family.mo_T1D_STRICT!=1)&(family.fa_T1D_STRICT==1),
    (family.mo_T1D_STRICT==1)&(family.fa_T1D_STRICT==1)
],[0,1,2,3])

family['group_aid'] = np.select([
    (family.mo_other_AIDs!=1)&(family.fa_other_AIDs!=1),
    (family.mo_other_AIDs==1)&(family.fa_other_AIDs!=1),
    (family.mo_other_AIDs!=1)&(family.fa_other_AIDs==1),
    (family.mo_other_AIDs==1)&(family.fa_other_AIDs==1)
],[0,1,2,3])

family.to_csv('family_20240202.csv',index=None)