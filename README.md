# GWAS
## 1 背景
### 1.1 GWAS 简介
GWAS全称 Genome-wide association study 全基因组关联分析，即研究表型（关注的性状）和基因型（变异）之间的**关联**，试图找到影响表型的差异的遗传因素

+ 思路

对多个个体全基因组范围的遗传变异性进行检测，获得基因型，进而将基因型与可观测的性状（表型）进行群体水平的统计学分析，根据检验标准筛选出最有可能影响该形状的遗传变异

### 1.2 Bonferroni 校正
我们在进行假设检验的时候，通常会设置一个零假设，之后计算出一个 p 值，即数据分布符合原假设的概率，p 值越低，即代表拒绝原假设的概率越大。我们通常认为 p ＜ 0.05 是一个判断是否显著的阈值。

**但是**，由于GWAS是对基因组上的多个遗传因素进行检验，在同时对多组数据进行处理和比较的时候，很可能其中部分数据因为随机效应而超过阈值，造成假阳性结果。而检验的次数越多，出现假阳性的概率就越大，因此急需一种方法来对结果的阈值进行校正。

+ 原理

Bonferroni校正即为最严格的多重检验矫正方法。在同一数据集上同时检验n个相互独立的假设，那么用于每一假设的统计显著水平，应为仅检验一个假设时的显著水平的1/n。如以显著水平0.05检验同一数据集上两个独立的假设，此时用于检验该两个假设应使用更严格的0.025；对于10000个SNP的检验，若将p设置为1e-6，进行10000次比较之后犯错误的概率是10-6*10000 = 0.01，严格地控制了假阳性的出现。

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

mkdir info pheno ref
# 基因型文件
wget http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz


# 表型文件
wget http://1002genomes.u-strasbg.fr/files/phenoMatrix_35ConditionsNormalizedByYPD.tab.gz -P ./pheno
gzip -d ./pheno/phenoMatrix_35ConditionsNormalizedByYPD.tab.gz
# 包含了 971 个分离株在 35 中胁迫条件下的生长率（30℃，YPD 培养基）


# 分离株信息文件
wget http://1002genomes.u-strasbg.fr/isolates/page8/files/1002genomes.txt -P info
# 包含了1011 株分离株的信息（ID、株系名、生态来源、地理来源等）


# 用于 GWAS 的矩阵
wget http://1002genomes.u-strasbg.fr/files/1011GWASMatrix.tar.gz -P ref
gzip -d ./ref/1011GWASMatrix.tar.gz
# 包含了 .bed、.bim 和 .fam 三个文件
```

### 2.2 GWAS 常用文件格式
![](./Fig/file.png)

常用 Plink 进行基因型与表型的关联分析，Plink 常用的文件格式有两套：map/ped 和 bim/fam/bed，两组文件均没有列名，且每一列表示的意思是一定的。

>.map文件主要是图谱文件信息，包括染色体编号、SNP名称、染色体的摩尔位置（可选项，可以用0）、SNP的物理位置
>
>.ped文件主要包括SNP的信息，包括 Family ID(没有可以用个体 ID 代替)、个体 ID、父本编号、母本编号、性别（未知用0）和表型数据
>
>.bim文件储存每个遗传变异的相关信息，每行代表一个遗传变异，共六列（染色体位置、遗传变异的编号、遗传变异在基因组上的摩尔位置、碱基对的坐标、等位基因1、等位基因2）
>
>.fam文件储存样本家系信息，共六列（家系编号、个体编号、父系编号（0表示缺失）、母系编号、性别编号、表型值（-9表示缺失））
>
>.bed文件储存基因型信息


### 2.3 数据准备和预处理
+ Plink 下载
```bash
brew install plink2
```

+ 简化数据方便计算
```bash
cd ref
sed -i '1013,1042d' 1002genomes.txt # 删除多余信息
tsv-summarize -H --count -g 4 1002genomes.txt | 
    mlr --itsv --omd cat # 检查不同分离株的来源

# 为了简化计算量，我们选取 Beer、Tree 和 Human, clinical 这三种来源的酵母数据进行操作
for i in Beer Tree 'Human, clinical'; do
    tsv-filter --str-eq 4:"$i" 1002genomes.txt |
    tsv-select -f 1,4 >> select_info.tsv
done 

wc -l select_info.tsv
# 230
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

+ 数据质量控制

质控一般包含两个方向：[一个是样本的质量控制（缺失率 < 5%;杂合性等）；另一个是 SNP 位点的质量控制（MAF > 5%; 哈温平衡检验等）](https://zhuanlan.zhihu.com/p/149947873?from_voters_page=true)

```bash
# 筛选出其中的 SNP 位点
brew install bcftools
bcftools index --threads 4 1011Matrix.gvcf.gz
bcftools view -v snps 1011Matrix.gvcf.gz > snp.gvcf 

# 筛选出我们需要的样本信息
bcftools view --threads 4 -S <(cat ./info/select_info.tsv | cut -f 1) snp.gvcf > tem&&
    mv tem snp.gvcf

# 筛选出 MAF > 5% 的点



```
> 1. 为什么对MAF进行过滤
> MAF:minor allele frequency,次等位基因频率；某个一个位点有AA或AT或TT，那么就可以计算A的基因频率和T的基因频率，qA + qT = 1，这里谁比较小，谁就是最小等位基因频率，qA = 0.3，qT = 0.7，那么这个位点的 MAF 为 0.3（如果一个位点有三个等位基因，那么频率排在中间的是 MAF；如果一个位点有四个等位基因，那么频率）。之所以用这个过滤标准，是因为 MAF 如果非常小，那么意味着大部分位点都是相同的基因型，这些位点贡献的信息非常少，放在计算中增加计算量，增加了假阳性的可能。
>
>
>

## 3 参考
[1. GWAS 分析](https://zhuanlan.zhihu.com/p/158869408)

[2. GWAS分析基本流程及分析思路](https://www.jianshu.com/p/f27c620d0bb2)

[3. 什么是Bonferroni校正？| 群体遗传专题](https://zhuanlan.zhihu.com/p/440376273)

[4. plink格式文件的介绍及相互转换](https://blog.csdn.net/qq_22253901/article/details/121608557)



