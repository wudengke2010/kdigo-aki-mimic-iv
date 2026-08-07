#!/usr/bin/env python3
"""
MIMIC-IV 脊髓 vs 肾损伤 论文可行性分析
"""
import pandas as pd
import gzip
import os

BASE = "E:/mimic-iv"

# =========================================
# 1. 基础表行数统计
# =========================================
print("=" * 70)
print("1. 基础表行数统计")
print("=" * 70)

core_tables = {}
for root, dirs, files in os.walk(os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1")):
    for f in files:
        if f.endswith('.csv.gz'):
            path = os.path.join(root, f)
            rel = os.path.relpath(path, os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1"))
            core_tables[rel] = path

# Key tables to load
HOSP = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/hosp")
ICU = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/icu")

def read_csv_gz(path, **kwargs):
    return pd.read_csv(path, compression='gzip', **kwargs)

# =========================================
# 2. 患者和入院基础信息
# =========================================
print("=" * 70)
print("2. 患者和入院基础数据")
print("=" * 70)

patients = read_csv_gz(os.path.join(HOSP, "patients.csv.gz"))
admissions = read_csv_gz(os.path.join(HOSP, "admissions.csv.gz"))
icustays = read_csv_gz(os.path.join(ICU, "icustays.csv.gz"))

print(f"  patients: {len(patients):,} rows")
print(f"  admissions: {len(admissions):,} rows, unique patients: {admissions['subject_id'].nunique():,}")
print(f"  icustays: {len(icustays):,} rows, unique patients: {icustays['subject_id'].nunique():,}")

# =========================================
# 3. 脊髓相关 ICD 诊断
# =========================================
print("\n" + "=" * 70)
print("3. 脊髓相关 ICD 诊断分析")
print("=" * 70)

diagnoses = read_csv_gz(os.path.join(HOSP, "diagnoses_icd.csv.gz"))
d_icd = read_csv_gz(os.path.join(HOSP, "d_icd_diagnoses.csv.gz"))

print(f"  diagnoses_icd 总行数: {len(diagnoses):,}")
print(f"  ICD 版本分布: {diagnoses['icd_version'].value_counts().to_dict()}")

# 脊髓 ICD-10
spinal_icd10_patterns = [
    (r'^S14', 'Cervical spinal cord injury'),
    (r'^S24', 'Thoracic spinal cord injury'),
    (r'^S34', 'Lumbar/sacral spinal cord injury'),
    (r'^G82', 'Paraplegia/Tetraplegia'),
    (r'^G95', 'Other spinal cord diseases'),
    (r'^M50\.[0-9]', 'Cervical disc with myelopathy'),
    (r'^M51\.[0-9]', 'Thoracic/lumbar disc with myelopathy'),
    (r'^T09\.3', 'Spinal cord injury unspecified'),
    (r'^T91\.[13]', 'Sequelae spinal cord injury'),
    (r'^Q05', 'Spina bifida'),
    (r'^Q06', 'Congenital spinal cord anomalies'),
    (r'^G83', 'Other paralytic syndromes'),
    (r'^M48\.0[0-6]', 'Spinal stenosis'),
    (r'^M47\.[1-9]', 'Spondylosis with myelopathy'),
]

spinal_icd9_patterns = [
    (r'^806', 'Fracture vertebral column + SCI'),
    (r'^952', 'SCI without fracture'),
    (r'^344\.[0-1]', 'Quadriplegia'),
    (r'^344\.1', 'Paraplegia'),
    (r'^336', 'Other spinal cord disease'),
    (r'^721\.1', 'Cervical spondylosis with myelopathy'),
    (r'^722\.7', 'Disc disorder with myelopathy'),
    (r'^724\.0[0-9]', 'Spinal stenosis'),
    (r'^907\.2', 'Late effect SCI'),
]

def count_icd_matches(diag, d_icd, patterns, version):
    sub = diag[diag['icd_version'] == version]
    results = {}
    all_indices = []
    for pattern, label in patterns:
        mask = sub['icd_code'].str.match(pattern, na=False)
        matched = sub[mask]
        if len(matched) > 0:
            results[label] = {
                'patients': matched['subject_id'].nunique(),
                'admissions': matched['hadm_id'].nunique(),
                'rows': len(matched),
                'codes': matched['icd_code'].nunique()
            }
            all_indices.append(matched.index)
    # 去重汇总
    if all_indices:
        all_matched = pd.concat([sub.loc[idx] for idx in all_indices]).drop_duplicates()
    else:
        all_matched = pd.DataFrame()
    return results, all_matched

print("\n  --- ICD-10 脊髓 ---")
r10, m10 = count_icd_matches(diagnoses, d_icd, spinal_icd10_patterns, 10)
for label, info in r10.items():
    print(f"  {label}: {info['patients']:,} pts, {info['admissions']:,} admits, {info['rows']:,} rows, {info['codes']} codes")

print("\n  --- ICD-9 脊髓 ---")
r9, m9 = count_icd_matches(diagnoses, d_icd, spinal_icd9_patterns, 9)
for label, info in r9.items():
    print(f"  {label}: {info['patients']:,} pts, {info['admissions']:,} admits, {info['rows']:,} rows, {info['codes']} codes")

spinal_all = pd.concat([m10, m9]).drop_duplicates('subject_id') if len(m10) > 0 or len(m9) > 0 else pd.DataFrame()
spinal_ids = set(spinal_all['subject_id'].unique())
print(f"\n  >>> 脊髓合并患者数: {len(spinal_ids):,}")
print(f"  >>> 脊髓合并入院数: {spinal_all['hadm_id'].nunique():,}")

# ICU 交叉
spinal_icu = icustays[icustays['subject_id'].isin(spinal_ids)]
print(f"  >>> 其中入ICU: {spinal_icu['subject_id'].nunique():,} pts, {len(spinal_icu):,} stays")

# =========================================
# 4. 肾损伤相关 ICD
# =========================================
print("\n" + "=" * 70)
print("4. 肾损伤 ICD 诊断分析")
print("=" * 70)

renal_icd10_patterns = [
    (r'^N17', 'Acute kidney failure (AKI)'),
    (r'^N18', 'Chronic kidney disease (CKD)'),
    (r'^N19', 'Unspecified kidney failure'),
    (r'^Z49', 'Encounter for dialysis'),
    (r'^Z99\.2', 'Dependence on renal dialysis'),
    (r'^Z94\.0', 'Kidney transplant status'),
    (r'^I12', 'Hypertensive renal disease'),
    (r'^I13', 'Hypertensive heart+renal disease'),
    (r'^E10\.2', 'T1DM kidney complications'),
    (r'^E11\.2', 'T2DM kidney complications'),
    (r'^N99\.0', 'Postprocedural renal failure'),
]

renal_icd9_patterns = [
    (r'^584', 'Acute kidney failure'),
    (r'^585', 'Chronic kidney disease'),
    (r'^586', 'Renal failure unspecified'),
    (r'^V56', 'Encounter for dialysis'),
    (r'^V45\.1', 'Renal dialysis status'),
    (r'^V42\.0', 'Kidney transplant status'),
    (r'^403', 'Hypertensive renal disease'),
    (r'^404', 'Hypertensive heart+renal'),
]

print("\n  --- ICD-10 肾损伤 ---")
r10k, m10k = count_icd_matches(diagnoses, d_icd, renal_icd10_patterns, 10)
for label, info in r10k.items():
    print(f"  {label}: {info['patients']:,} pts, {info['admissions']:,} admits, {info['rows']:,} rows")

print("\n  --- ICD-9 肾损伤 ---")
r9k, m9k = count_icd_matches(diagnoses, d_icd, renal_icd9_patterns, 9)
for label, info in r9k.items():
    print(f"  {label}: {info['patients']:,} pts, {info['admissions']:,} admits, {info['rows']:,} rows")

renal_all = pd.concat([m10k, m9k]).drop_duplicates()
renal_ids = set(renal_all['subject_id'].unique())
print(f"\n  >>> 肾损伤合并患者数: {len(renal_ids):,}")
print(f"  >>> 肾损伤合并入院数: {renal_all['hadm_id'].nunique():,}")

renal_icu = icustays[icustays['subject_id'].isin(renal_ids)]
print(f"  >>> 其中入ICU: {renal_icu['subject_id'].nunique():,} pts, {len(renal_icu):,} stays")

# AKI / CKD 细分
aki10 = diagnoses[(diagnoses['icd_version']==10) & diagnoses['icd_code'].str.match(r'^N17', na=False)]
aki9 = diagnoses[(diagnoses['icd_version']==9) & diagnoses['icd_code'].str.match(r'^584', na=False)]
ckd10 = diagnoses[(diagnoses['icd_version']==10) & diagnoses['icd_code'].str.match(r'^N18', na=False)]
ckd9 = diagnoses[(diagnoses['icd_version']==9) & diagnoses['icd_code'].str.match(r'^585', na=False)]
dialysis10 = diagnoses[(diagnoses['icd_version']==10) & diagnoses['icd_code'].str.match(r'^Z49|^Z99\.2', na=False)]
aki_all = pd.concat([aki10, aki9]).drop_duplicates('subject_id')
ckd_all = pd.concat([ckd10, ckd9]).drop_duplicates('subject_id')

print(f"\n  --- 细分类 ---")
print(f"  AKI (N17/584): {len(aki_all):,} patients")
print(f"  CKD (N18/585): {len(ckd_all):,} patients")
print(f"  Dialysis (Z49/Z99.2): {dialysis10['subject_id'].nunique():,} patients")

# AKI on CKD (同时有N17+N18)
aki_ids = set(aki_all['subject_id'])
ckd_ids = set(ckd_all['subject_id'])
aki_on_ckd = aki_ids & ckd_ids
print(f"  AKI on CKD (同时有 AKI+CKD): {len(aki_on_ckd):,} patients")

# =========================================
# 5. 肾功能实验室检查
# =========================================
print("\n" + "=" * 70)
print("5. 肾功能实验室检查")
print("=" * 70)

d_lab = read_csv_gz(os.path.join(HOSP, "d_labitems.csv.gz"))
print(f"  d_labitems: {len(d_lab):,} items")

renal_lab_kw = ['creatinine', 'bun', 'urea nitrogen', 'gfr', 'egfr', 'cystatin',
                'urine output', 'dialysis', 'renal', 'kidney']
for kw in renal_lab_kw:
    matches = d_lab[d_lab['label'].str.lower().str.contains(kw, na=False)]
    if len(matches) > 0:
        for _, row in matches.head(3).iterrows():
            print(f"  [{kw}] itemid={row['itemid']}: {row['label']}")

# 估算 labevents 大小
labevents_path = os.path.join(HOSP, "labevents.csv.gz")
labevents_size = os.path.getsize(labevents_path)
print(f"\n  labevents.csv.gz 文件大小: {labevents_size / (1024**3):.1f} GB")

# 快速采样估算行数
import subprocess
result = subprocess.run(['wsl', 'zcat', labevents_path.replace('E:', '/mnt/e').replace('\\', '/'), '|', 'wc', '-l'],
                       capture_output=True, text=True, shell=True, timeout=120)
labevents_rows = None
try:
    labevents_rows = int(result.stdout.strip())
    print(f"  labevents 估计总行数: {labevents_rows:,}")
except:
    print(f"  labevents 行数计数失败，将使用采样")

# =========================================
# 6. 脊髓手术 (ICD Procedures)
# =========================================
print("\n" + "=" * 70)
print("6. 脊髓相关手术操作")
print("=" * 70)

procedures = read_csv_gz(os.path.join(HOSP, "procedures_icd.csv.gz"))
d_proc = read_csv_gz(os.path.join(HOSP, "d_icd_procedures.csv.gz"))
print(f"  procedures_icd: {len(procedures):,} rows")
print(f"  d_icd_procedures: {len(d_proc):,} codes")

# 关键词搜索
proc_kw = ['spine', 'spinal', 'vertebr', 'lamin', 'discect', 'fusion',
           'cord', 'myelo', 'canal', 'decompress', 'cervical', 'thoracic', 'lumbar']
proc_patient_set = set()
for kw in proc_kw:
    matches = d_proc[d_proc['long_title'].str.lower().str.contains(kw, na=False)]
    if len(matches) > 0:
        codes = matches['icd_code'].tolist()
        sub = procedures[procedures['icd_code'].isin(codes)]
        proc_patient_set.update(sub['subject_id'].unique())

print(f"  >>> 脊髓相关手术总患者: {len(proc_patient_set):,}")

# =========================================
# 7. ICU内肾损伤严重度 (chartevents 中的 urine output, dialysis)
# =========================================
print("\n" + "=" * 70)
print("7. ICU 数据 - 尿量和透析")
print("=" * 70)

d_items = read_csv_gz(os.path.join(ICU, "d_items.csv.gz"))
# 搜索尿量和透析相关 item
uo_kw = ['urine output', 'urine', 'foley', 'uop']
dial_kw = ['dialysis', 'ultrafiltrat', 'crrt', 'iabl', 'sled', 'cvvh', 'hemodialy']

for kw in uo_kw + dial_kw:
    matches = d_items[d_items['label'].str.lower().str.contains(kw, na=False)]
    if len(matches) > 0:
        for _, row in matches.head(3).iterrows():
            print(f"  [{kw}] itemid={row['itemid']}: {row['label']} ({row['category']})")

# ICU procedure events (CRRT, dialysis)
proc_events = read_csv_gz(os.path.join(ICU, "procedureevents.csv.gz"))
print(f"\n  procedureevents (ICU): {len(proc_events):,} rows")

# =========================================
# 8. ED 数据中脊髓/肾损伤
# =========================================
print("\n" + "=" * 70)
print("8. 急诊数据 (ED)")
print("=" * 70)

ED = os.path.join(BASE, "ed/physionet.org/files/mimic-iv-ed/2.2/ed")
edstays = read_csv_gz(os.path.join(ED, "edstays.csv.gz"))
ed_dx = read_csv_gz(os.path.join(ED, "diagnosis.csv.gz"))
print(f"  edstays: {len(edstays):,}")
print(f"  ed_diagnosis: {len(ed_dx):,}")

# 脊髓 ED
ed_spinal = ed_dx[ed_dx['icd_code'].str.match('|'.join([p[0] for p in spinal_icd10_patterns + spinal_icd9_patterns]), na=False)]
print(f"  ED脊髓诊断: {len(ed_spinal):,} rows, {ed_spinal['subject_id'].nunique():,} pts")

# 肾 ED
ed_renal = ed_dx[ed_dx['icd_code'].str.match('|'.join([p[0] for p in renal_icd10_patterns + renal_icd9_patterns]), na=False)]
print(f"  ED肾损伤诊断: {len(ed_renal):,} rows, {ed_renal['subject_id'].nunique():,} pts")

# =========================================
# 9. eICU 多中心数据
# =========================================
print("\n" + "=" * 70)
print("9. eICU 多中心数据")
print("=" * 70)

EICU = os.path.join(BASE, "eicu/physionet.org/files/eicu-crd/2.0")
eicu_patient = read_csv_gz(os.path.join(EICU, "patient.csv.gz"))
eicu_dx = read_csv_gz(os.path.join(EICU, "diagnosis.csv.gz"))
eicu_lab = read_csv_gz(os.path.join(EICU, "lab.csv.gz"))
eicu_intake = read_csv_gz(os.path.join(EICU, "intakeOutput.csv.gz"))

print(f"  eICU patients: {len(eicu_patient):,}")
print(f"  eICU diagnosis: {len(eicu_dx):,}")
print(f"  eICU lab: {len(eicu_lab):,}")
print(f"  eICU intakeOutput: {len(eicu_intake):,}")

# eICU 脊髓诊断
eicu_spinal = eicu_dx[eicu_dx['diagnosisstring'].str.lower().str.contains(
    'spinal cord|paraplegia|quadriplegia|tetraplegia|myelopathy|spine injury|SCI', na=False)]
print(f"\n  eICU 脊髓诊断: {len(eicu_spinal):,} rows, {eicu_spinal['patientunitstayid'].nunique():,} stays")

# eICU 肾损伤诊断
eicu_renal = eicu_dx[eicu_dx['diagnosisstring'].str.lower().str.contains(
    'acute kidney|acute renal|aki|ckd|chronic kidney|renal failure|kidney failure|dialysis|crrt|nephro', na=False)]
print(f"  eICU 肾损伤诊断: {len(eicu_renal):,} rows, {eicu_renal['patientunitstayid'].nunique():,} stays")

# =========================================
# 10. Notes 文本
# =========================================
print("\n" + "=" * 70)
print("10. Notes 文本数据")
print("=" * 70)

NOTE = os.path.join(BASE, "note/physionet.org/files/mimic-iv-note/2.2/note")
try:
    disc = read_csv_gz(os.path.join(NOTE, "discharge.csv.gz"), nrows=5)
    rad = read_csv_gz(os.path.join(NOTE, "radiology.csv.gz"), nrows=5)
    print(f"  出院小结: 有 (columns: {list(disc.columns)})")
    print(f"  放射报告: 有 (columns: {list(rad.columns)})")
except Exception as e:
    print(f"  Notes 读取: {e}")

# =========================================
# 最终总结
# =========================================
print("\n" + "=" * 70)
print("=== 最终分析总结 ===")
print("=" * 70)

print(f"""
【数据概览】
  MIMIC-IV v3.1 核心  : hosp(15表) + icu(8表), ~10 GB
  MIMIC-IV-Note v2.2  : 出院小结 + 放射报告, ~2 GB
  MIMIC-IV-ED v2.2    : 急诊, ~117 MB
  eICU-CRD v2.0       : 多中心ICU, ~5.2 GB
  MIMIC-IV-ECG v1.0   : 心电波形 (下载中)

【脊髓相关论文可行性】
  - ICD诊断脊髓患者  : {len(spinal_ids):,} 人
  - 其中入ICU         : {spinal_icu['subject_id'].nunique():,} 人  
  - 脊髓手术患者      : {len(proc_patient_set):,} 人
  - ED脊髓就诊        : {ed_spinal['subject_id'].nunique() if len(ed_spinal) > 0 else 0:,} 人
  - eICU脊髓诊断      : {eicu_spinal['patientunitstayid'].nunique() if len(eicu_spinal) > 0 else 0:,} stays
  - Notes文本支持     : 出院小结 + 放射报告

【肾损伤论文可行性】
  - 肾损伤ICD患者     : {len(renal_ids):,} 人
  - AKI 患者           : {len(aki_ids):,} 人
  - CKD 患者           : {len(ckd_ids):,} 人
  - AKI on CKD         : {len(aki_on_ckd):,} 人
  - 透析依赖           : {dialysis10['subject_id'].nunique():,} 人
  - 其中入ICU          : {renal_icu['subject_id'].nunique():,} 人
  - 丰富实验室         : creatinine, BUN, eGFR 等
  - ICU尿量数据        : 有 (chartevents)
  - ICU透析/CRRT       : 有 (procedureevents)
  - eICU肾损伤         : {eicu_renal['patientunitstayid'].nunique() if len(eicu_renal) > 0 else 0:,} stays
""")

# 保存关键统计为 CSV
stats = {
    'category': ['Total Patients', 'Total Admissions', 'Total ICU Stays',
                 'Spinal Cord Patients', 'Spinal Cord ICU', 'Spinal Surgery Patients',
                 'Renal Patients', 'AKI Patients', 'CKD Patients', 'AKI on CKD',
                 'Dialysis Patients', 'Renal ICU Patients'],
    'count': [len(patients), len(admissions), len(icustays),
              len(spinal_ids), spinal_icu['subject_id'].nunique(), len(proc_patient_set),
              len(renal_ids), len(aki_ids), len(ckd_ids), len(aki_on_ckd),
              dialysis10['subject_id'].nunique(), renal_icu['subject_id'].nunique()]
}
pd.DataFrame(stats).to_csv('mimic_stats.csv', index=False)
print("统计结果已保存到 mimic_stats.csv")
