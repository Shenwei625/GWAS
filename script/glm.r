#! usr/bin/R

args <- commandArgs(T)
DATA <- read.table(args[1], header = TRUE, sep = "\t")
model <- glm(formula = DATA[,3] ~ DATA[,2] + as.factor(Admixture_Group),data = DATA, family = quasipoisson(link = "log"), na.action = na.omit)
NAME <- paste(sep="", args[1], ".lst")
write.table(coef(summary(model))[,4],file = NAME, quote = FALSE, sep = "\t")

