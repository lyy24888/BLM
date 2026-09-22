# 📚 GWAS 学习笔记：原理、R语言实战与数据库资源

本笔记整合了 GWAS（全基因组关联分析）的核心原理、常见结果解读、完整的 R 语言实践流程，以及常用的数据库资源和下游分析方向，适合作为入门与查阅的参考资料。

---

## 第一部分：GWAS 基本原理

### 1.1 核心概念
GWAS 本质上是**对全基因组上每一个 SNP 分别做一次回归分析**，检验该 SNP 的基因型与表型是否存在统计学关联。

### 1.2 统计模型
- **连续型表型（如身高、BMI、血压）**：使用线性回归
$$
  Y = \beta_0 + \beta_1 G + \beta_2 \text{年龄} + \beta_3 \text{性别} + \beta_4 PC1 + \cdots + \varepsilon
$$
  其中 \(G\) 为 SNP 的基因型剂量（0, 1, 2），\(\beta_1\) 表示每多一个效应等位基因，表型平均改变多少。

- **二分类表型（如病例/对照）**：使用逻辑回归
$$
\log\left(\frac{P(Y=1)}{1-P(Y=1)}\right) = \beta_0 + \beta_1 G + \text{协变量}
$$
  结果通常报告为 $OR = exp(β₁)$。

- **假设检验**：$(H_0: \beta_1 = 0)$（即 SNP 与表型无关联）。

### 1.3 关键分析结果

| 指标                  | 含义                    |
| --------------------- | ----------------------- |
| **效应等位基因 (EA)** | 增加/降低表型的等位基因 |
| **EAF**               | 效应等位基因频率        |
| **Beta / OR**         | 效应大小                |
| **SE**                | 标准误，越小越精确      |
| **P 值**              | 关联显著性              |
| **N**                 | 样本量                  |

- **全基因组显著阈值**：$(P < 5 \times 10^{-8})$（Bonferroni 校正：约 100 万个独立常见变异，0.05 / 1,000,000 = 5e-8）。
- **提示性显著**：$(P < 1 \times 10^{-5})$，需独立队列复制验证。

### 1.4 经典可视化

- **曼哈顿图 (Manhattan Plot)**：横轴为染色体位置，纵轴为 $(-log_{10}(P))$。高峰表示该区域存在关联信号。
- **QQ 图**：评估 P 值整体分布。尾部上翘提示存在真实关联信号；整体抬高提示群体分层或技术偏差。
- **Lambda GC**：接近 1 较好，明显大于 1 需警惕假阳性。

> ⚠️ **重要提醒**：GWAS 发现的是**统计关联，不是因果**。真正因果变异和机制还需 fine-mapping、功能实验和独立队列验证。

---

## 第二部分：R语言全流程实践脚本

以下脚本覆盖**模拟数据生成 → 读取 → 预处理 → 曼哈顿图/QQ图可视化 → 显著位点提取 → 效应量置信区间计算**全流程。

### 2.1 环境准备
```r
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(data.table)
library(qqman)
library(ggplot2)
set.seed(42)
```

### 2.2 模拟数据生成
```r
n_snps <- 5000
gwas <- data.table(
  SNP   = paste0("rs", 1:n_snps),
  CHR   = sample(1:22, n_snps, replace = TRUE),
  BP    = sample(10000:250000000, n_snps),
  EA    = sample(c("A", "T", "C", "G"), n_snps, replace = TRUE),
  EAF   = round(runif(n_snps, 0.05, 0.95), 2),
  BETA  = round(rnorm(n_snps, 0, 0.05), 4),
  SE    = round(runif(n_snps, 0.01, 0.1), 4)
)
# 人为植入强信号
signal_idx <- c(100, 1500, 3200, 4800)
gwas[signal_idx, BETA := c(0.35, -0.28, 0.42, 0.22)]
gwas[signal_idx, SE   := c(0.04, 0.05, 0.06, 0.03)]
# 计算 Z 和 P
gwas[, Z := BETA / SE]
gwas[, P := 2 * pnorm(-abs(Z))]
gwas[, LOG10_P := -log10(P)]
```

### 2.3 预处理
```r
# 检查缺失值
print(colSums(is.na(gwas)))
# 按 P 值排序
setorder(gwas, P)
# 染色体类型转换
gwas[, CHR := as.numeric(CHR)]
```

### 2.4 可视化
```r
# 曼哈顿图 (qqman)
png("manhattan_plot.png", width = 1400, height = 600, res = 120)
manhattan(gwas, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          col = c("grey50", "grey20"),
          main = "Manhattan Plot (Simulated GWAS)")
dev.off()

# QQ 图
png("qq_plot.png", width = 600, height = 600, res = 120)
qq(gwas$P, main = "Q-Q Plot")
dev.off()

# ggplot2 增强版曼哈顿图
gwas[, CHR_FAC := factor(CHR, levels = 1:22)]
p_manhattan <- ggplot(gwas, aes(x = BP, y = LOG10_P, color = CHR_FAC)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed") +
  scale_color_manual(values = rep(c("#1f77b4", "#ff7f0e"), 11)) +
  labs(x = "Chromosome", y = expression(-log[10](P)), title = "Manhattan Plot") +
  theme_minimal() + theme(legend.position = "none")
ggsave("manhattan_ggplot2.png", p_manhattan, width = 12, height = 5, dpi = 150)
```

### 2.5 下游分析
```r
# 提取显著位点
sig_df <- gwas[P < 1e-5]
print(sig_df[, .(SNP, CHR, BP, EA, EAF, BETA, SE, P)])

# 计算 95% 置信区间
sig_df[, CI_lower := BETA - 1.96 * SE]
sig_df[, CI_upper := BETA + 1.96 * SE]

# 森林图
sig_df[, SNP := factor(SNP, levels = SNP[order(BETA)])]
p_forest <- ggplot(sig_df, aes(x = SNP, y = BETA)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "SNP", y = "Effect Size (Beta) with 95% CI") +
  theme_minimal()
ggsave("forest_plot.png", p_forest, width = 8, height = 4, dpi = 150)
```

---

## 第三部分：GWAS 常用数据库及数据下载方式

### 3.1 综合 GWAS 汇总数据库
| 数据库             | 网址                                 | 数据规模           | 访问方式                    | 特点                     |
| ------------------ | ------------------------------------ | ------------------ | --------------------------- | ------------------------ |
| **GWAS Catalog**   | https://www.ebi.ac.uk/gwas/          | 全球GWAS研究汇编   | 网页 + FTP + API            | 最全面，提供关联TSV下载  |
| **OpenGWAS (IEU)** | https://gwas.mrcieu.ac.uk/           | 50,000+ 汇总数据集 | 网页 + API + `ieugwasr` R包 | 手动整理，支持MR、共定位 |
| **FinnGen**        | https://finngen.gitbook.io/          | 芬兰人群GWAS       | 在线申请 + 邮件获取         | 北欧人群，表型丰富       |
| **GWAS Atlas**     | https://bigd.big.ac.cn/gwas/         | 跨物种GWAS         | 网页下载                    | 涵盖植物、动物等多物种   |
| **Pan-UKBB**       | https://pan.ukbb.broadinstitute.org/ | UK Biobank全表型   | 网页下载                    | 多 ancestry 的 UKBB 数据 |

### 3.2 队列/联盟 GWAS 数据库
| 联盟/队列      | 目标疾病/表型         | 样本量   | 下载链接                                                     |
| -------------- | --------------------- | -------- | ------------------------------------------------------------ |
| **UK Biobank** | 多种疾病和表型        | 500,000+ | https://sites.google.com/broadinstitute.org/ukbbgwasesults/home |
| **PGC**        | 精神疾病（11种）      | 400,000+ | https://www.med.unc.edu/pgc/results-and-downloads            |
| **GIANT**      | BMI、身高、腰围       | 100,000+ | http://portals.broadinstitute.org/collaboration/giant/       |
| **ENIGMA**     | 脑影像、精神/神经疾病 | 50,000+  | http://enigma.ini.usc.edu/research/download-enigma-gwas-results/ |
| **CHARGE**     | 房颤、血压、骨密度等  | 38,000+  | https://www.chargeconsortium.com/                            |

### 3.3 数据下载方式详解
- **GWAS Catalog**：网页搜索导出单条结果；批量下载整个Catalog；FTP下载每个研究的完整汇总统计。
- **OpenGWAS**：使用 `ieugwasr` R包程序化下载。
  ```r
  library(ieugwasr)
  datasets <- gwasinfo("ieu-a-2")           # 搜索数据集
  gwas_data <- associations(variants = c("rs123", "rs456"), id = "ieu-a-2")  # 提取指定变异
  full_data <- tophits(id = "ieu-a-2", pval = 1e-5)  # 获取top hits
  ```
- **FinnGen**：在线表单申请，邮件获取下载指引。
- **通用工具**：`wget` / `curl` 批量下载FTP；`data.table::fread()` 直接读取 `.tsv.gz`。

---

## 第四部分：下游分析扩展方向

| 分析类型           | R包/工具           | 用途                                    | 输入数据要求       |
| ------------------ | ------------------ | --------------------------------------- | ------------------ |
| **LD Clumping**    | PLINK              | 提取独立显著位点，去除LD冗余            | 基因型 + 参考面板  |
| **Fine-mapping**   | `susieR`           | 缩小因果变异可信集                      | 汇总统计 + LD矩阵  |
| **共定位**         | `coloc`            | 检验GWAS信号与eQTL/pQTL是否共享因果变异 | 两个性状的汇总统计 |
| **孟德尔随机化**   | `TwoSampleMR`      | 推断暴露→结局因果效应                   | 暴露和结局汇总统计 |
| **多基因风险评分** | `bigsnpr` / PRSice | 构建个体遗传风险预测                    | 训练集 + 验证集    |
| **功能注释**       | ANNOVAR / VEP      | 注释SNP基因位置和功能                   | SNP列表 (CHR:BP)   |

### 重点方法说明
- **共定位分析 (coloc)**：使用 `coloc.abf()` 或 `coloc.susie()`，输出 PP.H0~PP.H4 五个后验概率，PP.H4 高表示共定位证据强。
- **Fine-mapping**：推荐 `susieR::susie_rss()`，输入Z分数和LD矩阵即可获得可信集。

---

## 第五部分：学习要点总结

1. **GWAS = 逐SNP回归**：连续表型用线性回归看 Beta，病例对照用逻辑回归看 OR。
2. **关键结果**：效应大小、P值、频率，以及曼哈顿图、QQ图、独立显著位点。
3. **显著性阈值**：全基因组 \(5 \times 10^{-8}\)，提示性 \(1 \times 10^{-5}\)。
4. **关联 ≠ 因果**：需 fine-mapping、功能实验和独立队列验证。
5. **实战流程**：数据读取 → 预处理 → 可视化 → 下游分析（显著位点提取、CI计算、功能注释）。
6. **数据库资源**：GWAS Catalog、OpenGWAS、FinnGen、PGC、UK Biobank 等，支持网页、FTP、API、R包多种下载方式。
7. **下游扩展**：LD Clumping、Fine-mapping、共定位、孟德尔随机化、多基因风险评分、功能注释。

> 💡 **建议**：先用手上的模拟数据跑通 R 脚本，熟悉数据结构和可视化，再尝试从 OpenGWAS 或 GWAS Catalog 下载真实汇总统计进行实战。