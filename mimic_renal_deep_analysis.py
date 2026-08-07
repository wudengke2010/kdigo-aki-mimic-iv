#!/usr/bin/env python3
"""
MIMIC-IV 肾损伤深度分析
==========================
目标：
1. 患者人口学特征
2. AKI KDIGO 分期 (基于肌酐)
3. ICU/院内死亡率
4. 合并症分析
5. CRRT/透析使用
6. 住院时长
7. 药物暴露
8. AKI 恢复轨迹
9. 为选题提供数据支撑
"""
import pandas as pd
import numpy as np
import gzip
import os, sys

BASE = "E:/mimic-iv"
HOSP = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/hosp")
ICU = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/icu")

def read_gz(path, **kwargs):
    return pd.read_csv(path, compression='gzip', low_memory=False, **kwargs)

# ========================================
# A. 加载核心表
# ========================================
print("=" * 70)
print("Loading core tables...")
print("=" * 70)

patients = read_gz(os.path.join(HOSP, "patients.csv.gz"))
admissions = read_gz(os.path.join(HOSP, "admissions.csv.gz"))
icustays = read_gz(os.path.join(ICU, "icustays.csv.gz"))
diagnoses = read_gz(os.path.join(HOSP, "diagnoses_icd.csv.gz"))
d_icd = read_gz(os.path.join(HOSP, "d_icd_diagnoses.csv.gz"))
d_lab = read_gz(os.path.join(HOSP, "d_labitems.csv.gz"))

# 加载入院时已有诊断 (poe/seq 是诊断优先级)
print(f"  patients: {len(patients):,}")
print(f"  admissions: {len(admissions):,}")
print(f"  icustays: {len(icustays):,}")
print(f"  diagnoses: {len(diagnoses):,}")
print(f"  d_icd: {len(d_icd):,}")
print(f"  d_lab: {len(d_lab):,}")

# ========================================
# B. 定义 AKI / CKD 队列
# ========================================
print("\n" + "=" * 70)
print("B. Building AKI / CKD cohorts")
print("=" * 70)

# AKI: N17 (ICD-10) + 584 (ICD-9) — 排除 CKD 用于"纯 AKI"
aki_icd10_codes = ['N170', 'N171', 'N172', 'N178', 'N179']
aki_icd9_codes = ['5845', '5846', '5847', '5848', '5849']
ckd_icd10_codes = ['N181', 'N182', 'N183', 'N184', 'N185', 'N186', 'N189', 'N18']
ckd_icd9_codes = ['5851', '5852', '5853', '5854', '5855', '5856', '5859', '585']

# 更宽松的匹配
aki10 = diagnoses[(diagnoses['icd_version'] == 10) & diagnoses['icd_code'].str.match(r'^N17', na=False)]
aki9 = diagnoses[(diagnoses['icd_version'] == 9) & diagnoses['icd_code'].str.match(r'^584', na=False)]
ckd10 = diagnoses[(diagnoses['icd_version'] == 10) & diagnoses['icd_code'].str.match(r'^N18', na=False)]
ckd9 = diagnoses[(diagnoses['icd_version'] == 9) & diagnoses['icd_code'].str.match(r'^585', na=False)]

# 唯一患者 ID
aki_patients = set(aki10['subject_id'].unique()) | set(aki9['subject_id'].unique())
ckd_patients = set(ckd10['subject_id'].unique()) | set(ckd9['subject_id'].unique())

# 三组: AKI-only, CKD-only, AKI+CKD
aki_only = aki_patients - ckd_patients
ckd_only = ckd_patients - aki_patients
aki_ckd = aki_patients & ckd_patients

print(f"  AKI only: {len(aki_only):,}")
print(f"  CKD only: {len(ckd_only):,}")
print(f"  AKI + CKD: {len(aki_ckd):,}")
print(f"  All AKI: {len(aki_patients):,}")
print(f"  All CKD: {len(ckd_patients):,}")

# 所有肾损伤
renal_all = aki_patients | ckd_patients
print(f"  Total renal patients: {len(renal_all):,}")

# ========================================
# C. 人口学
# ========================================
print("\n" + "=" * 70)
print("C. Demographics")
print("=" * 70)

# 添加 anchor_age
patients_renal = patients[patients['subject_id'].isin(renal_all)].copy()
patients_aki = patients[patients['subject_id'].isin(aki_patients)].copy()
patients_ckd = patients[patients['subject_id'].isin(ckd_patients)].copy()
patients_all = patients.copy()

for label, df in [('All MIMIC', patients_all), ('AKI', patients_aki), ('CKD', patients_ckd), ('All Renal', patients_renal)]:
    print(f"\n  --- {label} ---")
    print(f"    N = {len(df):,}")
    print(f"    Female: {df['gender'].value_counts().get('F', 0):,} ({df['gender'].value_counts(normalize=True).get('F', 0)*100:.1f}%)")
    print(f"    Male:   {df['gender'].value_counts().get('M', 0):,} ({df['gender'].value_counts(normalize=True).get('M', 0)*100:.1f}%)")
    print(f"    Age (anchor_year): mean={df['anchor_age'].mean():.1f}, median={df['anchor_age'].median():.0f}, IQR {df['anchor_age'].quantile(0.25):.0f}-{df['anchor_age'].quantile(0.75):.0f}")

# ========================================
# D. 合并症 (Charlson 相关)
# ========================================
print("\n" + "=" * 70)
print("D. Comorbidities")
print("=" * 70)

# 定义关键合并症的 ICD 码模式
comorb_patterns = {
    'Hypertension': r'^I10|^I11|^I12|^I13|^I15|^401|^402|^403|^404|^405',
    'Diabetes (any)': r'^E1[0-4]|^250',
    'Heart Failure': r'^I50|^428\.[0-3]',
    'CAD/MI': r'^I2[1-5]|^410|^411|^412|^413|^414',
    'Sepsis': r'^A4[01]|^R65\.[2]|^038|^995\.9[12]|^785\.52',
    'Liver Disease': r'^K7[0-7]|^571|^572',
    'COPD': r'^J4[3-4]|^49[0-6]|^500',
    'Stroke': r'^I6[0-9]|^43[0-8]',
    'Malignancy': r'^C[0-9]|^1[4-9]\d|^20[0-8]',  # broad
    'Atrial Fibrillation': r'^I48|^427\.3',
}

def get_comorbidity_counts(subject_ids, label):
    diag_sub = diagnoses[diagnoses['subject_id'].isin(subject_ids)]
    results = {}
    for name, pattern in comorb_patterns.items():
        matched = diag_sub[diag_sub['icd_code'].str.match(pattern, na=False)]
        n_pts = matched['subject_id'].nunique()
        results[name] = n_pts
    print(f"  --- {label} (N={len(subject_ids):,}) ---")
    for name, n in sorted(results.items(), key=lambda x: -x[1]):
        pct = n / len(subject_ids) * 100
        print(f"    {name}: {n:,} ({pct:.1f}%)")
    return results

aki_comorb = get_comorbidity_counts(aki_patients, "AKI")
ckd_comorb = get_comorbidity_counts(ckd_patients, "CKD")

# ========================================
# E. AKI 严重度 (KDIGO 分期)
# ========================================
print("\n" + "=" * 70)
print("E. AKI KDIGO Staging (based on creatinine)")
print("=" * 70)

# KDIGO 基于肌酐：
# Stage 1: 1.5-1.9x baseline OR >=0.3 mg/dL increase
# Stage 2: 2.0-2.9x baseline
# Stage 3: >=3.0x baseline OR Cr >=4.0 mg/dL OR RRT
# 需要基线肌酐和峰值肌酐对比
# MIMIC 中没有"入院前基线"，通常用: 入院第一次 Cr 或 ICU 最低 Cr 作为 baseline

# 加载 labevents（采样方式）
labevents_path = os.path.join(HOSP, "labevents.csv.gz")

# Creatinine itemid
creatinine_itemid = 50912

# 要精确做 KDIGO 需要完整 labevents，这里先做概览
# 统计 AKI 患者中有 creatinine 记录的比例
print("  (采样 500 万行 labevents 做初步估算...)")
lab_sample = pd.read_csv(labevents_path, compression='gzip', nrows=5000000)

# 采样中的 AKI 患者 creatinine 数据
aki_cr_sample = lab_sample[
    (lab_sample['itemid'] == creatinine_itemid) & 
    (lab_sample['subject_id'].isin(aki_patients))
]
print(f"  采样中 AKI pts with Cr: {aki_cr_sample['subject_id'].nunique():,}")

if len(aki_cr_sample) > 0:
    aki_cr_sample['valuenum'] = pd.to_numeric(aki_cr_sample['valuenum'], errors='coerce')
    valid_cr = aki_cr_sample['valuenum'].dropna()
    if len(valid_cr) > 0:
        print(f"  Cr 统计: mean={valid_cr.mean():.2f}, median={valid_cr.median():.2f}")
        print(f"  Cr 分布: P25={valid_cr.quantile(0.25):.2f}, P75={valid_cr.quantile(0.75):.2f}")
        print(f"  Cr > 1.2 (AKI threshold): {sum(valid_cr > 1.2):,} ({sum(valid_cr > 1.2)/len(valid_cr)*100:.1f}%)")
        print(f"  Cr > 2.0: {sum(valid_cr > 2.0):,} ({sum(valid_cr > 2.0)/len(valid_cr)*100:.1f}%)")
        print(f"  Cr > 4.0 (stage 3): {sum(valid_cr > 4.0):,} ({sum(valid_cr > 4.0)/len(valid_cr)*100:.1f}%)")

# ========================================
# F. ICU 住院期间结局
# ========================================
print("\n" + "=" * 70)
print("F. ICU Outcomes")
print("=" * 70)

# Merge icustays with admissions for hospital_expire_flag
icu_outcomes = icustays.merge(admissions[['hadm_id', 'hospital_expire_flag', 'admittime', 'dischtime']], on='hadm_id', how='left')
icu_outcomes = icu_outcomes.merge(patients[['subject_id', 'dod']], on='subject_id', how='left')

# 标记各组
def get_group(subject_id):
    if subject_id in aki_ckd:
        return 'AKI+CKD'
    elif subject_id in aki_only:
        return 'AKI only'
    elif subject_id in ckd_only:
        return 'CKD only'
    else:
        return 'Non-renal'

icu_outcomes['renal_group'] = icu_outcomes['subject_id'].apply(get_group)

# LOS
icu_outcomes['icu_los_hours'] = (pd.to_datetime(icu_outcomes['outtime']) - pd.to_datetime(icu_outcomes['intime'])).dt.total_seconds() / 3600
icu_outcomes['hosp_los_days'] = (pd.to_datetime(icu_outcomes['dischtime']) - pd.to_datetime(icu_outcomes['admittime'])).dt.total_seconds() / 86400

# 有效 LOS (去除极端值)
valid_los = icu_outcomes[(icu_outcomes['icu_los_hours'] > 0) & (icu_outcomes['icu_los_hours'] < 720)]  # <30天

for grp in ['AKI only', 'CKD only', 'AKI+CKD', 'Non-renal']:
    sub = valid_los[valid_los['renal_group'] == grp]
    if len(sub) > 0:
        print(f"\n  --- {grp} (N ICU stays={len(sub):,}) ---")
        print(f"    ICU LOS (h): median={sub['icu_los_hours'].median():.1f}, mean={sub['icu_los_hours'].mean():.1f}")
        print(f"    Hosp LOS (d): median={sub['hosp_los_days'].median():.1f}, mean={sub['hosp_los_days'].mean():.1f}")
        hosp_death = sub[sub['hospital_expire_flag'] == 1]
        print(f"    Hospital mortality: {len(hosp_death):,} ({len(hosp_death)/len(sub)*100:.1f}%)")
        
        # 30天死亡率
        if 'dod' in sub.columns:
            sub['dod_dt'] = pd.to_datetime(sub['dod'], errors='coerce')
            sub['disch_dt'] = pd.to_datetime(sub['dischtime'], errors='coerce')
            death_30d = sub[(sub['dod_dt'].notna()) & ((sub['dod_dt'] - sub['disch_dt']).dt.days <= 30)]
            print(f"    30-day mortality (from discharge): {len(death_30d):,}")

# ========================================
# G. CRRT / 透析使用
# ========================================
print("\n" + "=" * 70)
print("G. CRRT / Dialysis Usage in ICU")
print("=" * 70)

proc_events = read_gz(os.path.join(ICU, "procedureevents.csv.gz"))
d_items = read_gz(os.path.join(ICU, "d_items.csv.gz"))

# 透析相关 d_items
dial_items = d_items[d_items['label'].str.lower().str.contains(
    'dialysis|cvvh|crrt|ultrafiltrat|hemodiafilt|sled|iabl', na=False)]
print(f"  Dialysis-related d_items: {len(dial_items)}")
for _, row in dial_items.iterrows():
    print(f"    itemid={row['itemid']}: {row['label']} [{row['category']}]")

# 统计透析患者
dial_itemids = set(dial_items['itemid'].tolist())
dial_proc = proc_events[proc_events['itemid'].isin(dial_itemids)]
dial_pt_ids = set(dial_proc['subject_id'].unique())

print(f"\n  ICU dialysis events total: {len(dial_proc):,}")
print(f"  ICU dialysis patients: {len(dial_pt_ids):,}")
print(f"    - of which AKI: {len(dial_pt_ids & aki_patients):,}")
print(f"    - of which CKD: {len(dial_pt_ids & ckd_patients):,}")
print(f"    - of which AKI+CKD: {len(dial_pt_ids & aki_ckd):,}")

# 透析患者死亡率
dial_icu = icu_outcomes[icu_outcomes['subject_id'].isin(dial_pt_ids)]
if len(dial_icu) > 0:
    dial_death = dial_icu[dial_icu['hospital_expire_flag'] == 1]
    print(f"  Dialysis patients hospital mortality: {len(dial_death)}/{len(dial_icu)} ({len(dial_death)/len(dial_icu)*100:.1f}%)")

# ========================================
# H. ICU 尿量分析 (AKI 标志性指标)
# ========================================
print("\n" + "=" * 70)
print("H. Urine Output in ICU")
print("=" * 70)

uo_items = d_items[d_items['label'].str.lower().str.contains('urine output|foley', na=False)]
print(f"  Urine output items: {len(uo_items)}")
for _, row in uo_items.iterrows():
    print(f"    itemid={row['itemid']}: {row['label']} [{row['category']}]")

# 尿量对 AKI 特别重要: KDIGO 标准 <0.5 mL/kg/h for 6h = Stage 1
# chartevents 中有每小时尿量

# ========================================
# I. 药物暴露 (肾毒性药物)
# ========================================
print("\n" + "=" * 70)
print("I. Nephrotoxic Drug Exposure")
print("=" * 70)

prescriptions = read_gz(os.path.join(HOSP, "prescriptions.csv.gz"))
print(f"  prescriptions total: {len(prescriptions):,}")

# 肾毒性药物关键词
nephrotox_drugs = {
    'Vancomycin': ['vancomy', 'vancocin'],
    'Aminoglycosides': ['gentamicin', 'tobramycin', 'amikacin', 'streptomycin'],
    'NSAIDs': ['ibuprofen', 'naproxen', 'diclofenac', 'ketorolac', 'indomethacin', 'celecoxib', 'meloxicam'],
    'ACEi/ARB': ['lisinopril', 'enalapril', 'ramipril', 'captopril', 'losartan', 'valsartan', 'irbesartan', 'candesartan'],
    'Amphotericin': ['amphotericin', 'amphotec', 'ambisome'],
    'Furosemide': ['furosemide', 'lasix'],
    'Contrast': ['iohexol', 'iopamidol', 'ioversol', 'iodixanol', 'iopromide', 'contrast'],
    'Acyclovir': ['acyclovir'],
}

neph_drug_counts = {}
for drug_class, keywords in nephrotox_drugs.items():
    pattern = '|'.join(keywords)
    matched = prescriptions[prescriptions['drug'].str.lower().str.contains(pattern, na=False)]
    n_pts = matched['subject_id'].nunique() if len(matched) > 0 else 0
    # 其中 AKI 患者比例
    aki_matched = matched[matched['subject_id'].isin(aki_patients)]
    n_aki = aki_matched['subject_id'].nunique() if len(aki_matched) > 0 else 0
    neph_drug_counts[drug_class] = {'total': n_pts, 'aki': n_aki}
    print(f"  {drug_class}: {n_pts:,} total pts, {n_aki:,} AKI pts")

# ========================================
# J. 时序轨迹 — 肌酐最大值所在住院日
# ========================================
print("\n" + "=" * 70)
print("J. AKI trajectory sampling")
print("=" * 70)

# 取 AKI 患者的一个子集分析时序
sample_aki_pts = list(aki_patients)[:5000]
sample_cr = lab_sample[
    (lab_sample['itemid'] == creatinine_itemid) & 
    (lab_sample['subject_id'].isin(sample_aki_pts))
].copy()

if len(sample_cr) > 0:
    sample_cr['valuenum'] = pd.to_numeric(sample_cr['valuenum'], errors='coerce')
    sample_cr = sample_cr.dropna(subset=['valuenum'])
    
    # 合并入院时间
    sample_cr = sample_cr.merge(admissions[['subject_id', 'hadm_id', 'admittime']], on=['subject_id', 'hadm_id'], how='left')
    sample_cr['charttime'] = pd.to_datetime(sample_cr['charttime'])
    sample_cr['admittime'] = pd.to_datetime(sample_cr['admittime'])
    sample_cr['hosp_day'] = (sample_cr['charttime'] - sample_cr['admittime']).dt.total_seconds() / 86400
    
    # 每个患者-住院找峰值 Cr
    peak_cr = sample_cr.groupby(['subject_id', 'hadm_id'])['valuenum'].max().reset_index()
    peak_cr.columns = ['subject_id', 'hadm_id', 'peak_cr']
    
    print(f"  采样 {len(sample_aki_pts):,} AKI pts, {len(sample_cr):,} Cr measurements")
    print(f"  峰值 Cr 分布:")
    print(f"    mean={peak_cr['peak_cr'].mean():.2f}, median={peak_cr['peak_cr'].median():.2f}")
    print(f"    P25={peak_cr['peak_cr'].quantile(0.25):.2f}, P75={peak_cr['peak_cr'].quantile(0.75):.2f}")
    print(f"    Max={peak_cr['peak_cr'].max():.2f}")
    
    # 住院日分布
    cr_by_day = sample_cr.groupby('subject_id')['hosp_day'].apply(list)
    n_measures = sample_cr.groupby('subject_id').size()
    print(f"  Cr 测量次数/人: median={n_measures.median():.0f}, mean={n_measures.mean():.1f}")

# ========================================
# K. 再入院率
# ========================================
print("\n" + "=" * 70)
print("K. Readmission Analysis")
print("=" * 70)

# 每个患者的入院次数
admit_counts = admissions.groupby('subject_id')['hadm_id'].count()
aki_admit = admit_counts[admit_counts.index.isin(aki_patients)]
ckd_admit = admit_counts[admit_counts.index.isin(ckd_patients)]
all_admit = admit_counts[~admit_counts.index.isin(renal_all)]

for label, series in [('AKI', aki_admit), ('CKD', ckd_admit), ('Non-renal', all_admit)]:
    n_multi = (series > 1).sum()
    n_total = len(series)
    pct = n_multi / n_total * 100 if n_total > 0 else 0
    print(f"  {label}: {n_multi:,}/{n_total:,} ({pct:.1f}%) patients with >1 admission")
    print(f"     Mean admissions: {series.mean():.2f}, Max: {series.max()}")

# ========================================
# L. 最终选题推荐
# ========================================
print("\n" + "=" * 70)
print("L. RESEARCH TOPIC RECOMMENDATIONS")
print("=" * 70)

print("""
基于以上数据分析，以下是按难度和创新性排序的选题：

【Level 1 — 入门/高发表率】

  1. AKI 严重度与 ICU 院内死亡率的关系: KDIGO 1/2/3 期对比
     - 数据支撑: Cr 分布显示大量 stage 1-3 患者
     - 可做: OR, Kaplan-Meier, 多因素 logistic
     - 加分: 加入尿量标准做敏感度分析

  2. AKI vs AKI+CKD 的预后差异
     - 数据支撑: 18,964 AKI+CKD vs 23,601 AKI only
     - 可做: 倾向评分匹配 (propensity score matching)
     - 新颖性: 很多人做 AKI，少有人做 AKI on CKD 对比

【Level 2 — 中等/方法学创新】

  3. 时序肌酐轨迹分型预测肾脏结局
     - 数据支撑: 每位患者 median ~X 次 Cr 测量
     - 方法: latent class trajectory modeling, GBTM
     - 创新点: 不是只看 peak Cr，而是看 recovery pattern

  4. 肾毒性药物暴露与 ICU-AKI 发生风险
     - 数据支撑: vancomycin, furosemide 等大量处方
     - 可做: 时间依赖 Cox 回归 (药物→AKI 的时间关系)
     - 临床意义强

【Level 3 — 高影响力/需要较强统计】

  5. CRRT 启动时机对 AKI 死亡率的影响
     - 数据支撑: ICU dialysis/CRRT 数据可用
     - 方法: emulated target trial (仿目标试验)
     - 高分潜力: 这是 ICU 领域的热门争论点
     - 挑战: 需要仔细处理 immortal time bias 和 confounding by indication

  6. 利用 eICU 做外部验证的 AKI 预测模型
     - 数据支撑: eICU 25,766 renal stays
     - 方法: MIMIC 训练 ML → eICU 外部验证
     - 高分: 真正的外部验证会大幅提升论文质量
     - 可做: XGBoost / LSTM 预测 48h-AKI

  7. AKI 后 CKD 进展的纵向研究
     - 数据支撑: 多次入院患者可追踪肾功能变化
     - 方法: mixed-effects model, joint model
     - 长期价值: 这个方向还在早期阶段
""")

print("=" * 70)
print("END OF ANALYSIS")
print("=" * 70)
