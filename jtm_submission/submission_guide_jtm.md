# Journal of Translational Medicine 投稿操作指引

> 更新日期：2026-09-23 | 投稿包位置：`jtm_submission/` | 稿件：manuscript_jtm.docx

---

## 一、投稿前须知（已核实）

| 项目 | 详情 |
|------|------|
| 投稿入口 | https://translational-medicine.biomedcentral.com/ → 首页右侧 **"Submit manuscript"** |
| 投稿系统 | Editorial Manager（需注册账号，通讯作者邮箱 wudk2010@csu.edu.cn） |
| Article type | **Research article** |
| 期刊 Section | 提交时选择 **"Pain, Critical Care and Anesthesia"**（本文最佳匹配）；备选 "Data-driven Clinical Decision Processes" |
| 审稿模式 | 单盲（single-anonymous），审稿人知道作者身份 |
| 速度 | 官方中位 **9 天**首决；网友经验约 1 周送审、86 天 submission-to-accept |
| **APC（录用后付）** | **£3,050 / $4,090 / €3,350**（2026-09 官网现价，注意比早期估计 $2,890 已上调） |
| APC 减免 | 低收入国家自动减免；其他情况需**在投稿时**提出 case-by-case 申请（录用后再申请无效） |

## 二、上传文件清单

| 顺序 | 文件 | 系统中的 designation |
|:---:|------|------|
| 1 | `manuscript_jtm.docx`（主稿，含 Title page + 摘要 + 正文 + Declarations + References + Tables + Figure legends） | Manuscript |
| 2 | `cover_letter_jtm.md` → 上传时直接粘贴进 Cover Letter 输入框（或转成 PDF 上传） | Cover Letter |
| 3 | `strobe_checklist_jtm.docx` | STROBE Checklist（若有该 designation；否则作 Additional file / Report） |
| 4–10 | `figures/Figure1_flowchart.png` … `FigureS3_phenotype_forest.png`（7 张） | Figure（逐张上传） |
| 11 | 补充材料（可选打包）：Tables S1–S4 + 补充报告 → `v2_outputs/sensitivity_report.md`、`table_S1.md` | Additional files / Supplementary Material |

## 三、Editorial Manager 分步操作

1. **登录/注册**：用 wudk2010@csu.edu.cn 注册（若已有 BMC/Springer Nature 账号可直接登录；注意 EM 账号与 BMC 网站账号独立）。
2. **Submit a new manuscript** → 选择 Article Type: **Research Article**。
3. **填写基本信息**：
   - Title：照抄稿件标题（含副标题）
   - Abstract：从 DOCX 复制（349 词，含 Background/Methods/Results/Conclusions 小标题）
   - Keywords：6 个（Acute kidney injury; KDIGO; Mortality; Transient AKI; Persistent AKI; Critical care）
   - Section：选 **Pain, Critical Care and Anesthesia**
4. **作者信息**：逐位添加。**投稿后作者姓名/顺序/单位不可再改**，务必核对：
   - Jiqiang Liu — 单位 1,2
   - Dengke Wu（Corresponding，勾选）— 单位 1,2，wudk2010@csu.edu.cn
5. **上传文件**（按上表顺序），每个文件选择正确 designation。
6. **Cover Letter** 输入框粘贴 `cover_letter_jtm.md` 全文（保留分段加粗格式）。
7. **回答系统问题**（常见必答项）：
   - 数据可用性：选 "Data included in a publicly available repository / all data are available" → 填 PhysioNet 两个链接 + GitHub 仓库链接
   - Ethics：研究使用已去标识化公开数据，IRB 批准且豁免知情同意
   - Competing interests：None
   - Funding：照抄稿件 Funding 段三个基金号
   - Reporting guidelines：STROBE（+ TRIPOD for external validation）
8. **建议审稿人（可选但有帮助）**：Cover Letter 已留有空间；如填写需真实可验证（机构邮箱/ORCID），可建议 2–3 名 AKI 预后研究方向的国内/国际学者。
9. **Build PDF & Approve**：系统合成 PDF → 逐页检查（表格、图片、上标）→ Approve → 完成投稿。投稿号（如 JTM-D-26-xxxxx）会即时生成。

## 四、合规自查（已全部满足 ✅）

- ✅ 摘要 349 词（≤350），四段式结构
- ✅ 关键词 6 个（3–10 之间）
- ✅ Title page：完整作者名 + 单位 + 通讯作者
- ✅ Declarations 全部 8+1 分节（Ethics/Consent/Availability/Competing/Funding/Contributions/Acknowledgements/Authors' info/GenAI）
- ✅ LLM 使用声明已在 **Methods 末尾**（JTM 硬性要求）+ Declarations（双保险）
- ✅ 彩图 ≥300 dpi（7 张全部达标）
- ⚠️ 线条图（Figure 1 流程图）官方建议 ≥1000 dpi，目前 350 dpi；BMC 实际接收 PNG 通常不因线条图 dpi 桌拒，若系统技术检查退回，需用矢量 PDF 重导高分辨率版本（源文件 `v2_outputs/fig1_flowchart_v2.png` 由 matplotlib 生成，可一键改 dpi 重跑）
- ✅ 数据可用性：PhysioNet 凭证访问 + GitHub 全代码公开
- ✅ 无一稿多投（CKJ 2026-09-07 桌拒、Renal Failure 2026-09-21 桌拒，均已结案）

## 五、投稿后状态解读

| 状态 | 含义 | 预期 |
|------|------|------|
| New submission → Technical check | 格式/文件核查 | 1–3 天 |
| Editor assigned / With editor | 编辑评估是否送审 | 3–10 天（官方中位 9 天首决） |
| Under review | 外审中 | 4–8 周 |
| Decision (reject / revise / accept) | — | — |

**若桌拒**：立即启动路线下一站 **Annals of Intensive Care**（1区Top, IF 6.9）或 **J Crit Care**（IF 3.5，免 APC）。稿件为 BMC 格式，投 AIC（Springer）需微调，投 J Crit Care（Elsevier）需转格式——届时告知即可。

## 六、录用后付款

APC $4,090 在录用后才会收到 invoice；付款前可确认单位是否有 Springer Nature 开放获取协议（中南大学可能参与 transformative agreement，可覆盖部分费用）。查询：https://www.springernature.com/gp/open-research/institutional-agreements
