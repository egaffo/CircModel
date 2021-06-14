meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/DM1/analyses/meta_DM1.csv")
meta.data = as.data.table(meta.data)
meta.data <- meta.data[, .(sample_id = sample,
                           condition = ifelse(disease_class=="myotonic dystrophy type 1", "DM1","Normal"))]
meta.data = meta.data[order(meta.data$sample_id),][seq(1,nrow(meta.data), by = 2),]     # for odd rows

coldata <- data.frame(group = meta.data$condition,
                     sample = meta.data$sample,
                     row.names = meta.data$sample)
coldata$group <- factor(coldata$group)
coldata$sample <- as.character(coldata$sample)
coldata

library(caret)
grp1_name = "Normal"
grp2_name = "DM1"
variable_name = "group"

printIdx <- function() {
  NormSub <- sample(coldata$sample, size = 10, replace = FALSE)
  tumorSub <- sample(setdiff(coldata$sample, NormSub), size = 10, replace = FALSE)
  idx <- c(NormSub, tumorSub)
  
  idx
  
}

set.seed(5)
randomSubsets <- t(replicate(50, printIdx()))

# set.seed(5)
# numSubsets <- t(replicate(30, numIdx()))
# all(dist(numSubsets) > 0)

write.table(randomSubsets,file="/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/robustness_glmm/DM1/random_shuffle50.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
