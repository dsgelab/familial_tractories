import logging
import pandas as pd

# select the endpoints you plan to look into
# a list of ADs: https://risteys.finngen.fi/phenocode/AUTOIMMUNE
eps = ['T1D_STRICT', 'M13_RHEUMA', 'M13_RELAPSPOLYCHONDR', 'M13_SJOGREN', 'M13_SYSTSLCE', 'M13_DERMATOPOLY', # E4_DM1
       'M13_WEGENER', 'M13_MICROPOLYANG', 'M13_CHURGSTRAUSS', 'D3_ALLERGPURPURA', 'M13_BEHCET', 'M13_MCTD',
       'M13_HYPERANG', 'SLE_FG',  #'M13_SLE',
       'I9_RHEUFEV', 'G6_MS', 'G6_ADEM', 'G6_DISSOTH', 'G6_NARCOCATA', 'AUTOIMMUNE_HYPERTHYROIDISM',
       'E4_THYROIDITAUTOIM', 'E4_AUTOPOLYFAI', 'E4_HYTHY_AI_STRICT', 'E4_GRAVES_OPHT_STRICT', 'E4_ADDISON',
       'AUTOHEP', 'D3_AIHA_DRUG', 'D3_AIHA_OTHER', 'D3_ITP', 'D3_ANAEMIA_B12_DEF', 'K11_COELIAC', 'K11_IBD',
       'G6_MYASTHENIA', 'G6_OTHDEMYEL', 'G6_MYOMUSCINOTH', 'G6_GUILBAR', 'H7_IRIDOCYC_ANTER',  'CHIRBIL_PRIM',
       'L12_PSORIASIS', 'L12_VITILIGO', 'L12_ALOPECAREATA', 'L12_PEMPHIGOID', 'L12_DERMATHERP',
       'N14_HENOCHSCHONLEIN_NEPHRITIS', 'N14_IGA_NEPHROPATHY', 'T2D', 'GEST_DIABETES', 'K11_CROHN', 'J10_ASTHMA']
# eps = ['T1D_STRICT', 'E4_THYROIDITAUTOIM', 'K11_COELIAC', 'D3_ANAEMIA_B12_DEF', 'M13_RHEUMA',
#        'L12_VITILIGO', 'GRAVES_OPHT', 'K11_CROHN']

who_dict = {'child': 'ch', 'mother': 'mo', 'father': 'fa', 'parent': 'pa'}


def load_data(event_path, info_path, geo_path, pedigree_path):
    logging.info('Data loading ...')
    # Get first events
    df_events = pd.read_csv(event_path)
    df_events = df_events.rename(columns={'FINNGENID': 'ID'})
    # Get demographic information
    df_info = pd.read_csv(info_path)
    df_info['ch_year'] = df_info['date_of_birth'].str.split('-').str[0].astype(float)
    df_info = df_info.rename(columns={'FINREGISTRYID': 'ID'})
    # Get geographic information
    df_geo = pd.read_csv(geo_path).rename(columns={'FINREGISTRYID': 'ID'})
    df_info = df_info.merge(df_geo[['ID', 'maakunta']], 'left', on='ID')
    # Load pedigree
    pedigree = pd.read_csv(pedigree_path, sep='\t')
    # merge pedigree info into df_info
    logging.info('Data is loaded.')
    return df_events, df_info, pedigree
