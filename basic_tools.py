import logging
import pandas as pd

who_dict = {'child': 'ch', 'mother': 'mo', 'father': 'fa', 'parent': 'pa'}

# select the endpoints you plan to look into
# a list of ADs: https://risteys.finngen.fi/phenocode/AUTOIMMUNE
eps = ['T1D_STRICT', 'M13_RHEUMA', 'M13_RELAPSPOLYCHONDR', 'M13_SJOGREN', 'M13_SYSTSLCE', 'M13_DERMATOPOLY',  # E4_DM1
       'M13_WEGENER', 'M13_MICROPOLYANG', 'M13_CHURGSTRAUSS', 'D3_ALLERGPURPURA', 'M13_BEHCET', 'M13_MCTD',
       'M13_HYPERANG', 'SLE_FG', 'I9_RHEUFEV', 'G6_MS', 'G6_ADEM', 'G6_DISSOTH', 'G6_NARCOCATA',
       'AUTOIMMUNE_HYPERTHYROIDISM', 'E4_THYROIDITAUTOIM', 'E4_AUTOPOLYFAI', 'E4_HYTHY_AI_STRICT',
       'E4_GRAVES_OPHT_STRICT', 'E4_ADDISON', 'AUTOHEP', 'D3_AIHA_DRUG', 'D3_AIHA_OTHER', 'D3_ITP',
       'D3_ANAEMIA_B12_DEF', 'K11_COELIAC', 'K11_IBD', 'G6_MYASTHENIA', 'G6_OTHDEMYEL', 'G6_MYOMUSCINOTH',
       'G6_GUILBAR', 'H7_IRIDOCYC_ANTER', 'CHIRBIL_PRIM', 'L12_PSORIASIS', 'L12_VITILIGO', 'L12_ALOPECAREATA',
       'L12_PEMPHIGOID', 'L12_DERMATHERP', 'N14_HENOCHSCHONLEIN_NEPHRITIS', 'N14_IGA_NEPHROPATHY',
       'T2D', 'GEST_DIABETES', 'K11_CROHN', 'J10_ASTHMACOPDKELA']

# convert eps to interpretable names
eps_dict = {'T1D_STRICT': 'Type 1 diabetes',
            'M13_RHEUMA': 'Rheumatoid arthritis',
            'M13_SJOGREN': 'Sjögren syndrome',
            'M13_SYSTSLCE': 'Systemic sclerosis',
            'M13_WEGENER': 'Wegener granulomatosis',
            'D3_ALLERGPURPURA': 'Allergic purpura',
            'M13_MCTD': 'Mixed connective tissue disease',
            'SLE_FG': 'Systemic lupus erythematosus',
            'G6_MS': 'Multiple Sclerosis',
            'G6_DISSOTH': 'Autoimmune acute disseminated demyelination',
            'AUTOIMMUNE_HYPERTHYROIDISM': 'Autoimmune hyperthyroidism',
            'E4_THYROIDITAUTOIM': 'Autoimmune thyroiditis',
            'E4_HYTHY_AI_STRICT': 'Autoimmune hypothyroidism',
            'M13_DERMATOPOLY': 'Dermatopolymyositis',
            'E4_GRAVES_OPHT_STRICT': 'Graves opthalmopathy',
            'L12_PEMPHIGOID': 'Pemphigoid',
            'I9_RHEUFEV': 'Rheumatic fever incl heart disease',
            'E4_ADDISON': 'Adrenocortical insufficiency',
            'D3_AIHA_OTHER': 'Autoimmune haemolytic anaemias',
            'D3_ITP': 'Idiopathic thrombocytopenic purpura',
            'D3_ANAEMIA_B12_DEF': 'Vitamin B12 deficiency anaemia',
            'K11_COELIAC': 'Coeliac disease',
            'K11_IBD': 'Inflammatory bowel disease',
            'G6_MYASTHENIA': 'Myasthenia gravis',
            'G6_OTHDEMYEL': 'Autoimmune demyelinating diseases',
            'G6_GUILBAR': 'Guillain-Barre syndrome',
            'H7_IRIDOCYC_ANTER': 'Anterior Iridocyclitis',
            'CHIRBIL_PRIM': 'Primary biliary cholangitis',
            'L12_PSORIASIS': 'Psoriasis',
            'L12_VITILIGO': 'Vitiligo',
            'L12_ALOPECAREATA': 'Alopecia areata',
            'L12_DERMATHERP': 'Dermatitis herpetiformis',
            'N14_IGA_NEPHROPATHY': 'IgA nephropathy',
            'T2D': 'Type 2 diabetes',
            'GEST_DIABETES': 'Gestational diabetes',
            'K11_CROHN': 'Crohn disease',
            'J10_ASTHMA': 'Asthma',
            'J10_ASTHMACOPDKELA': 'Asthma'}

eps_selected = ['AUTOIMMUNE_HYPERTHYROIDISM', 'CHIRBIL_PRIM', 'D3_AIHA_OTHER', 'D3_ALLERGPURPURA', 'D3_ANAEMIA_B12_DEF',
                'D3_ITP', 'E4_ADDISON', 'E4_HYTHY_AI_STRICT', 'E4_THYROIDITAUTOIM', 'G6_DISSOTH', 'G6_GUILBAR', 'G6_MS',
                'G6_MYASTHENIA', 'G6_MYOMUSCINOTH', 'G6_OTHDEMYEL', 'GEST_DIABETES', 'H7_IRIDOCYC_ANTER',
                'J10_ASTHMACOPDKELA', 'K11_COELIAC', 'K11_CROHN', 'K11_IBD', 'L12_ALOPECAREATA', 'L12_DERMATHERP',
                'L12_PSORIASIS', 'L12_VITILIGO', 'M13_MCTD', 'M13_RHEUMA', 'M13_SJOGREN', 'M13_SYSTSLCE', 'M13_WEGENER',
                'N14_IGA_NEPHROPATHY', 'SLE_FG', 'T1D_STRICT', 'T2D']


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
