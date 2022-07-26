import tqdm
import warnings
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# import imputed HLA data
hla_data = pd.read_csv('R9_imputed_HLAs_v2.csv')
columns = hla_data.ID.tolist()
hla_data = hla_data.transpose()
hla_data.columns = columns
hla_data = hla_data.iloc[9:, :]

# convert to dosage probability
hla_df = hla_data[['A*01:01']]
for i in tqdm.tqdm(columns):
    snp = hla_data[i].str.extract(':(.+),(.+),(.+)')
    snp['dose_str'] = hla_data[i].extract('(.+):')
    snp['dose_int'] = np.select([(snp.dose_str == './.'), (snp.dose_str == '0/0'), (snp.dose_str == '0/1'),
                                 (snp.dose_str == '1/1')], [-1.0, 0.0, 1.0, 2.0])
    a = snp[0].astype(float)*0 + snp[1].astype(float)*1 + snp[2].astype(float)*2
    a1, a2 = a[~a.isna()], a[a.isna()]
    a2 = snp[snp.index.isin(a2.keys())].dose_int
    hla_df[i] = a1.append(a2)

# add outcome col
events = pd.read_csv('/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_endpoint_longitudinal_1.0.txt.gz',
                    sep='\t')
events = events[events.ENDPOINT == 'T1D_STRICT']
hla_df['outcome'] = np.select([(hla_df.index.isin(events.FINNGENID)), (~hla_df.index.isin(events.FINNGENID))], [1, 0])

# import confounding cols
confounding = pd.read_csv('confounding.csv')
confounding = confounding[confounding.PC1 != 0]
confounding = confounding[~confounding.age.isna()]
hla_df['finngen_id'] = hla_df.index
hla_df = hla_df.merge(confounding, 'inner', on='finngen_id')

# build a for-loop for logistic regression
res_df = pd.DataFrame(columns=['snp', 'coef', 'se', 'pval', 'hr_025', 'hr_975'])
for snp in tqdm.tqdm(columns):
    x = hla_df[[snp]+confounding.columns[1:]]
    x = sm.add_constant(x)
    y = hla_df[['outcome']]
    model = sm.Logit(y, x).fit(disp=0)
    res = model.summary2().talbes[1]
    hr_025 = np.exp(res.loc[snp, '[0.025'])
    hr_975 = np.exp(res.loc[snp, '0.975]'])
    pval = res.loc[snp, 'P>|z|']
    se = res.loc[snp, 'Std.Err.']
    coef = res.loc[snp, 'Coef.']
    row = [snp, coef, se, pval, hr_025, hr_975]
    res_df = res_df.append(pd.Series(row, index=res_df.columns), ignore_index=True)

res_df['or'] = np.exp(res_df.coef)
res_df.to_csv('hla_res.csv', index=None)

