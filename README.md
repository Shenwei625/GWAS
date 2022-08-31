# GWAS
## 1 背景
### 1.1 GWAS 简介
GWAS全称 Genome-wide association study 全基因组关联分析，即研究表型（关注的性状）和基因型（变异）之间的**关联**，试图找到影响表型的差异的遗传因素

+ 思路

对多个个体再全基因组范围的遗传变异性进行检测，获得基因型，进而将基因型与可观测的性状（表型）进行群体水平的统计学分析，根据检验标准筛选出最有可能影响该形状的遗传变异

### 1.2 Bonferroni 校正
我们再进行假设检验的时候，通常会设置一个零假设，之后计算出一个 p 值，即数据分布符合原假设的概率，p 值越低，即代表拒绝原假设的概率越大。我们通常认为 p ＜ 0.05 是一个判断是否显著的阈值。

**但是**，由于GWAS是对基因组上的多个遗传因素进行检验，在同时对多组数据进行处理和比较的时候，很可能其中部分数据因为随机效应而超过阈值，造成假阳性结果。而检验的次数越多，出现假阳性的概率就越大，因此急需一种方法来对结果的阈值进行校正。

+ 原理

Bonferroni校正即为最严格的多重检验矫正方法。在同一数据集上同时检验n个相互独立的假设，那么用于每一假设的统计显著水平，应为仅检验一个假设时的显著水平的1/n。如以显著水平0.05检验同一数据集上两个独立的假设，此时用于检验该两个假设应使用更严格的0.025；对于10000个基因的检验，若将p设置为1e-6，进行10000次比较之后犯错误的概率是10-6*10000 = 0.01，严格地控制了假阳性的出现。

由于GWAS标记之间的连锁不平衡，可能会存在多个标记或者SNP之间相互连锁的情况，也就是说它们之间的分布并不是完全独立的，所以假设GWAS数据集的每个关联测试都是独立的是不正确的。因此，应用Bonferroni校正通常会为我们提供最保守的p值阈值，其中可能会出现假阴性的情况，我们往往需要根据实际曼哈顿图的情况对阈值进行一些调整。

### 1.3 GWAS曼哈顿图分析
![](./Fig/MHD.jpg)

图片来源于文章《1,135 Genomes Reveal the Global Pattern of Polymorphism in *Arabidopsis thaliana*》，分析了SNP 与不同温度下拟南芥（10 与 16℃）开花时间之间的关联，黑色和灰色的点表示 1001 数据库中 SNP 的数据，彩色的点表示 RegMap 数据库中的数据；虚线表示 Bonferroni 校正后的 p 值，点线表示 permutation 校正（相比于 Bonferroni 校正松）后的 p 值。

观察左上角的图中的黑色与灰色的点（1001 数据库中拟南芥 SNP 与拟南芥 10℃ 开花时间的关联），我们可以发现在 1 号染色体上于拟南芥 10℃ 开花时间相关显著的 SNP 存在于 *FT*基因上；2 号染色体上可能相关的 SNP 存在于*SVP*基因上；同理 5 号染色体上*DOG1*和*VIN3*上的 SNP 与开花时间存在显著关联。


## 2 流程
### 2.1 数据下载
+ 使用文章[《Genome evolution across 1,011 *Saccharomyces cerevisiae* isolates》](https://doi.org/10.1038/s41586-018-0030-5)中的数据 

```bash
mkdir data
cd data

# 基因型文件
wget http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz
gzip -d 1011Matrix.gvcf.gz


# 表型文件
wget http://1002genomes.u-strasbg.fr/files/phenoMatrix_35ConditionsNormalizedByYPD.tab.gz
gzip -d phenoMatrix_35ConditionsNormalizedByYPD.tab.gz
# 包含了 971 个分离株在 35 中胁迫条件下的生长率（30℃，YPD 培养基）


# 分离株信息文件
wget http://1002genomes.u-strasbg.fr/isolates/page8/files/1002genomes.txt
# 包含了1011 株分离株的信息（ID、株系名、生态来源、地理来源等）


```

### 数据准备和预处理

```bash
sed -i '1013,1042d' 1002genomes.txt # 删除多余信息
tsv-summarize -H --count -g 4 1002genomes.txt | 
    mlr --itsv --omd cat # 检查不同分离株的来源
```
| Ecological origins | count |
| --- | --- |
| Wine | 248 |
| Beer | 59 |
| Sake | 47 |
| Unknown | 28 |
| Distillery | 29 |
| Bakery | 37 |
| Human, clinical | 107 |
| Soil | 38 |
| Industrial | 30 |
| Fruit | 47 |
| Nature | 52 |
| Water | 19 |
| Dairy | 27 |
| Tree | 64 |
| Fermentation | 36 |
| Cider | 17 |
| Palm wine | 30 |
| Human | 31 |
| Insect | 20 |
| Flower | 14 |
| Probiotic | 2 |
| Bioethanol | 27 |
| Lab strain | 2 |




## 3 参考
[1. GWAS 分析](https://zhuanlan.zhihu.com/p/158869408)

[2. 什么是Bonferroni校正？| 群体遗传专题](https://zhuanlan.zhihu.com/p/440376273)

[3. Fast-LMM](https://github.com/fastlmm/FaST-LMM)

