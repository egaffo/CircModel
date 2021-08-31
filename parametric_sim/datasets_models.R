library(edgeR)
library(DESeq2)

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

by.group=coldata$condition
refdesign= model.matrix(~condition, coldata)

## load data previously filtered
load("/blackhole/alessia/CircModel/data/IPFData_list.RData")

IPFData_models <- lapply(IPFData_list, function(count.matrix){
  # count.matrix = IPFData_list$ccp2
  by.sample <- DGEList(count.matrix)
  dds <- DESeqDataSetFromMatrix(countData = ceiling(count.matrix[,rownames(coldata)]),
                                colData = coldata,
                                design = ~ condition)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  sf <- sizeFactors(dds)
  # Estimate the NB dispersion using deseq2 method
  dds <- estimateDispersions(dds, fitType = "local",
                             minmu = 1e-6)
  # Estimate coefficient 
  dds <- nbinomWaldTest(dds)
  # Estimate the log-overall mean
  centered.off <- getOffset(by.sample)
  centered.off <- centered.off - mean(centered.off)
  logmeans <- mglmOneGroup(by.sample$counts, offset = centered.off, dispersion = dispersions(dds))
  dispersion = dispersions(dds)
  # Estimate of dispersion with ZINB model
  zinb.prop <- rep(-Inf, nrow(count.matrix))
  zinb.disp <- dispersions(dds)
  zinb.mean <- exp(logmeans)
  
  library(pscl)
  for (i in which(rowSums(by.sample$counts==0)>0)){
    try({
      # i = 1
      zfit <- zeroinfl(by.sample$counts[i,] ~ 1 | 1, 
                       dist = "negbin", 
                       offset = log(sf),
                       EM = F,
                       start = list(count = coef(dds)[i,1],
                                    zero = -3,
                                    theta = 1/dispersions(dds)[i]))
      
      zinb.mean[i] <- mean(exp(zfit$coefficients$count))
      zinb.prop[i] <- zfit$coefficients$zero
      zinb.disp[i] <- 1/zfit$theta
      
    })
  }
  
  zinb.prop <- exp(zinb.prop)/(1+exp(zinb.prop))
  
  return(list(sample.mu = logmeans, sample.disp = dispersion, 
              zi.mu = zinb.mean, zi.disp = zinb.disp, zi.prop = zinb.prop, 
              relative.size = by.sample$samples$lib.size))
})
save(IPFData_models,file = "/blackhole/alessia/CircModel/data/IPFZinb_nb_Fit_detmet_models.RData")


