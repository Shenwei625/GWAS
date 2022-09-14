# 使用 R 进行 GWAS 分析
## 1 背景
### 1.1 泊松分布
泊松分布适合于描述单位时间内随机事件发生的次数，如某一服务设施在一定时间内到达的人数，电话交换机接到呼叫的次数，汽车站台的候客人数，机器出现的故障数等。（事件是独立发生的；事件发生的概率在给定的固定时间内不随时间变化）

概率函数：P(X = k) = (λ<sup>k</sup>/k!)e<sup>-λ</sup>

λ 表示一段时间内事件发生的平均次数；k 表示一段时间内事件发生的次数。

### 1.2 二项分布
二项分布即重复n次独立的伯努利试验。在每次试验中只有两种可能的结果，而且两种结果发生与否互相对立，并且相互独立，与其它各次试验结果无关，事件发生与否的概率在每一次独立试验中都保持不变，则这一系列试验总称为n重伯努利实验

+ 二项分布与泊松分布的关系

当二项分布的n很大时，泊松分布可以作为二项分布的近似，常当n≧20,p≦0.05时，就可以用泊松公式近似得计算

### 1.3 伽马分布
是统计学中的连续概率函数，要等到n个随机事件都发生，需要经历多久时间。


### 1.4 glm 函数结果解读
```
Call:  glm(formula = AT1G10920 ~ chr1_10000545, family = Gamma(link = "inverse"), 
    data = DATA, na.action = na.omit)
# 套用的公式

Coefficients:
  (Intercept)  chr1_10000545  
    0.1518219      0.0001042  
# Coefficients：回归系数
# intercept：常数项

Degrees of Freedom: 451 Total (i.e. Null);  450 Residual
  (因为不存在，20个观察量被删除了)
Null Deviance:	    1.116 
Residual Deviance: 1.115 	AIC: 267.3
# Degrees of Freedom：自由度，计算某一统计量时，取值不受限制的变量个数
# Null Deviance：无效偏差，是指仅包括截距项、不包括解释变量的模型和饱和模型比较得到的偏差统计量的值
# Residual Deviance：残差，是指既包括截距项，又包括解释变量的模型和饱和模型比较得到的偏差统计量的值
# AIC值：Akaike information criterion，是衡量统计模型拟合优良性的一种标准，AIC越小，模型越好
```


## 2 数据预处理和导入

+ 使用 210217.SuppDataSet2.DiploidUnitNumberCalls.tsv 表型文件 Rexp.tsv
```bash
# 将缺失值改为 NA
sed -i 's/\t\t/\tNA\t/g' 210217.SuppDataSet2.DiploidUnitNumberCalls.tsv
cat 210217.SuppDataSet2.DiploidUnitNumberCalls.tsv | perl -ne'
    s/\t\n/\tNA\n/g;
    print "$_";
' > tem&&
mv tem 210217.SuppDataSet2.DiploidUnitNumberCalls.tsv

# 修改表头
sed -i 's/^0/SAMPLE/' 210217.SuppDataSet2.DiploidUnitNumberCalls.tsv
sed -i 's/^\t/SAMPLE\t/' Rexp.tsv

tsv-join -H --filter-file Rexp.tsv --key-fields SAMPLE --append-fields AT1G10920 <(tsv-select -H --fields SAMPLE,chr1_10000545 210217.SuppDataSet2.DiploidUnitNumberCalls.tsv) > test.tsv
```

```R
# 用法：glm(Y ~ X1 + X2 + X3,family = poisson(link = "log"), data = data, na.action=na.omit)

FILE <- "test.tsv"
DATA <- read.table(FILE, header = TRUE, sep = "\t")


```

