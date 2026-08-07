"""Quick check of temporal dynamics using mortality_30d/90d/1y columns."""
import pandas as pd
import numpy as np
import statsmodels.api as sm

df = pd.read_csv('aki_paradox_cohort_with_sofa.csv')
aki = df[df['max_aki_stage'] > 0].copy()
s3 = aki[aki['max_aki_stage'] == 3]

for col, label in [('hospital_expire_flag', 'Hospital'), ('mortality_30d', '30-day'),
                    ('mortality_90d', '90-day'), ('mortality_1y', '1-year')]:
    if col in s3.columns:
        pure = s3[s3['ckd'] == 0]
        ckd = s3[s3['ckd'] == 1]
        pm = pure[col].mean() * 100
        cm = ckd[col].mean() * 100
        X = sm.add_constant(s3['ckd'])
        m = sm.Logit(s3[col], X).fit(disp=0)
        or_val = np.exp(m.params['ckd'])
        ci = np.exp(m.conf_int().loc['ckd'])
        p = m.pvalues['ckd']
        print(f"  {label}: Pure={pm:.1f}%, CKD={cm:.1f}%, OR={or_val:.3f} ({ci[0]:.3f}-{ci[1]:.3f}), P={p:.3f}")
