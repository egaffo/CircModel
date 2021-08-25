#------------------------------------------------------------------------------------------------
###### DM1 data set 

meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/DM1/meta_DM1.csv")
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
NormIdx <- which(coldata$group == "Normal")
TumorIdx <- which(coldata$group == "DM1")
# half_1_2 <- round(c(length(NormIdx)/2, length(TumorIdx)/2))
# names(half_1_2) <- c(grp1_name,grp2_name)
# If one of index lists is bigger than the other, we select only the first n indexes for both lists
# So index_1 and index_2 lists are of equal lengths
# min_length <- min(length(NormIdx),length(TumorIdx))
# index_1 <- NormIdx[1:min_length]
# index_2 <- TumorIdx[1:min_length]
# half_1_2 <- rep(min_length/2,2)
# names(half_1_2) <- c(grp1_name,grp2_name)
# coldata.filt <- coldata[c(index_1,index_2),]
# re-compute indexes for grp1 and grp2
index_1 <- which(coldata[,variable_name] == grp1_name)
index_2 <- which(coldata[,variable_name] == grp2_name)
printIdx <- function() {
  NormSub <- sample(index_1, size = 3, replace = FALSE)
  tumorSub <- sample(index_2, size = 3, replace = FALSE)
  idx <- c(NormSub, tumorSub)
  
  index_1bis = index_1[-idx]
  index_2bis = index_2[-idx]
  
  NormSub2 <- sample(index_1bis, size = 2, replace = FALSE)
  tumorSub2 <- sample(index_2bis, size = 6, replace = FALSE)
  idxbis = c(NormSub2, tumorSub2)
  
  c(coldata$sample[idx], coldata$sample[idxbis])
}

set.seed(5)
randomSubsets <- t(replicate(30, printIdx()))

# set.seed(5)
# numSubsets <- t(replicate(30, numIdx()))
# all(dist(numSubsets) > 0)

write.table(randomSubsets,file="/blackhole/alessia/CircModel/power/random_subsets_eval_veri.txt", 
            quote=FALSE,row.names=FALSE,col.names=FALSE)

#------------------------------------------------------------------------------------------------
###### IPF data set 

meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/IPF/analyses/meta_IPF.csv")
meta.data = meta.data[seq(1,nrow(meta.data), by = 2),]     # for odd rows

coldata <- DataFrame(condition = meta.data$condition,
                     group = ifelse(meta.data$condition=="normal", "normal", "IPF"),
                     sample = meta.data$sample,
                     row.names = meta.data$sample)
coldata$condition <- factor(coldata$condition)
coldata$group <- factor(coldata$group)
coldata$sample <- as.character(coldata$sample)
coldata

library(caret)
grp1_name = "normal"
grp2_name = "IPF"
variable_name = "group"
NormIdx <- which(coldata$group == "normal")
TumorIdx <- which(coldata$group == "IPF")
# half_1_2 <- round(c(length(NormIdx)/2, length(TumorIdx)/2))
# names(half_1_2) <- c(grp1_name,grp2_name)
# If one of index lists is bigger than the other, we select only the first n indexes for both lists
# So index_1 and index_2 lists are of equal lengths
# min_length <- min(length(NormIdx),length(TumorIdx))
# index_1 <- NormIdx[1:min_length]
# index_2 <- TumorIdx[1:min_length]
# half_1_2 <- rep(min_length/2,2)
# names(half_1_2) <- c(grp1_name,grp2_name)
# coldata.filt <- coldata[c(index_1,index_2),]
# re-compute indexes for grp1 and grp2
index_1 <- which(coldata[,variable_name] == grp1_name)
index_2 <- which(coldata[,variable_name] == grp2_name)
printIdx <- function() {
  NormSub <- sample(index_1, size = 3, replace = FALSE)
  tumorSub <- sample(index_2, size = 4, replace = FALSE)
  idx <- c(NormSub, tumorSub)
  
  index_1bis = index_1[-idx]
  index_2bis = index_2[-idx]
  
  NormSub2 <- sample(index_1bis, size = 4, replace = FALSE)
  tumorSub2 <- sample(index_2bis, size = 4, replace = FALSE)
  idxbis = c(NormSub2, tumorSub2)
  
  c(coldata$sample[idx], coldata$sample[idxbis])
}

set.seed(5)
randomSubsets <- t(replicate(30, printIdx()))

write.table(randomSubsets,file="/blackhole/alessia/CircModel/power/IPF_random_subsets_eval_veri.txt", 
            quote=FALSE,row.names=FALSE,col.names=FALSE)
