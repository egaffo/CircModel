# ------------------------
## load libraries
# ------------------------

library("GenomicRanges")
library("Biobase")
library("BiocParallel")
pkgs <- c("edgeR", 
          "glmGamPoi",
          "DESeq2",
          "phyloseq", 
          "plyr", 
          "reshape2",
          "ROCR",
          # "samr",
          "apeglm",
          "ZIM",
          "zinbwave",
          "AUC",
          "genefilter",
          "MAST",
          # "scde", # It is important to use flexmix v2.3-13 and scde 1.99.1
          # "Seurat",
          "crayon",
          "aod",
          "arm",
          "fdrtool",
          "lars",
          "emdist",
          # "baySeq",
          "snow",
          # "ShrinkBayes",
          # "DEsingle",
          "VGAM",
          "Matrix",
          "maxLik",
          "MASS",
          # "samr",
          #"EBseq",
          "NewWave",
          "limma",
          "data.table",
          "dplyr")
for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }

library(devtools)
install_github("egaffo/CREART", ref = "dev", force = T)
library(CREART)

# ------------------------------------------------------------------------
## load data and meta data from previous data formatting in GOF.Rmd step
# ------------------------------------------------------------------------


## DM1 data set
randomSubsets <- read.table("/blackhole/alessia/CircModel/power/random_subsets_eval_veri.txt",strings=FALSE)

meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/DM1/meta_DM1.csv")
meta.data = as.data.table(meta.data)
meta.data <- meta.data[, .(sample_id = sample,
                           condition = ifelse(disease_class=="myotonic dystrophy type 1", "DM1","Normal"))]
meta.data = meta.data[order(meta.data$sample_id),][seq(1,nrow(meta.data), by = 2),]     # for odd rows

coldata <- DataFrame(group = meta.data$condition,
                     sample = meta.data$sample,
                     row.names = meta.data$sample)
coldata$group <- factor(coldata$group)
coldata$sample <- as.character(coldata$sample)
coldata

## IPF data set
randomSubsets <- read.table("/blackhole/alessia/CircModel/power/IPF_random_subsets_eval_veri.txt",strings=FALSE)

meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/IPF/analyses/meta_IPF.csv")
meta.data = meta.data[seq(1,nrow(meta.data), by = 2),]     # for odd rows

coldata <- DataFrame(condition = meta.data$condition,
                     group = ifelse(meta.data$condition=="normal", "normal", "IPF"),
                     sample_id = meta.data$sample,
                     row.names = meta.data$sample)
coldata$condition <- factor(coldata$condition)
coldata$group <- factor(coldata$group)
coldata$sample_id <- as.character(coldata$sample_id)
coldata

#-----------------------------------------
## summarized experiment for single matrix
## only ccp2 output is of interest
#-----------------------------------------

# se <- lapply(lapply(Data_list, function(x) x[,match(randomSubsets[1,],colnames(x))]), function(x) SummarizedExperiment(x, colData = coldata[rownames(coldata)%in%randomSubsets[1,],]))
# se <- lapply(DM1Data_list, function(x) SummarizedExperiment(x[,rownames(coldata)], 
#                                                             colData = coldata))
# eset <- lapply(se, function(x) ExpressionSet(assay(x),
#                                              AnnotatedDataFrame(as.data.frame(colData(x)))))
# 
# e <- lapply(eset, function (x) {
#   pData(x)$condition <- factor(pData(x)$group)
#   levels(pData(x)$condition) <- c("A", "B")
#   x.new = x[which(rowSums(exprs(x))!=0),]
#   return(x.new)})

ccp2 = read.table("/blackhole/alessia/CircModel/data/IPF_ccp2.csv", header = T)
# ccp2 = RCurl::scp(host = "threesum", path = "/home/enrico/analysis/zicirc/IPF/data/IPF_ccp2.csv", keypasswd = "alessiasucks666", 
#            user="alessia")
# load("/blackhole/alessia/CircModel/data/IPFData_list.RData")

# ccp2 = IPFData_list$ccp2
chr <- sub(":.*", "", ccp2$circ_id)
start <- as.numeric(gsub(".*:(.*)\\-.*", "\\1", ccp2$circ_id))
end <- sub(".*-", "\\1", ccp2$circ_id)
circ_id = paste0(chr, ":", start+1, "-", end)
ccp2$circ_id = circ_id
ccp2.df = ccp2[,-1]
rownames(ccp2.df) = ccp2$circ_id
se.ccp2 <- SummarizedExperiment(assays = list(counts = as.matrix(ccp2.df)),
                                colData = coldata)
eset.ccp2 <- ExpressionSet(assay(se.ccp2),
                           AnnotatedDataFrame(as.data.frame(colData(se.ccp2))))
pData(eset.ccp2)$condition <- factor(pData(eset.ccp2)$group)
levels(pData(eset.ccp2)$condition) <- c("A", "B")
e = eset.ccp2

#-----------------------------------------
## summarized experiment for GLMM matrix
#-----------------------------------------

## DM1
load("/blackhole/alessia/CircModel/data/DM1Data_list.RData")
glmm.db <- rbindlist(lapply(DM1Data_list[1:6], function(x) data.frame(x, circ_id = rownames(x))), 
                     idcol = "method", use.names = TRUE)
glmm.melt <- rbindlist(lapply(DM1Data_list[1:6], function(x) reshape2::melt(x)), 
                       idcol = "method", use.names = TRUE)

count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
name.sep <- "."
count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample")

colData.dt <- as.data.table(count.data.merge)[, .N, by = .(variable, method, 
                                                           group,
                                                           sample.name.ext)][, N := NULL][]

colData <- data.frame(colData.dt, 
                      row.names = "sample.name.ext")

glmm.wide = dcast(glmm.melt, method+Var2~Var1, value.var = "value", fill=0)
colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
head(glmm.wide[,c(1:8)])

## IPF
glmm.db <- rbindlist(lapply(IPFData_list[1:6], function(x) data.frame(x, circ_id = rownames(x))), 
                     idcol = "method", use.names = TRUE)
glmm.melt <- rbindlist(lapply(IPFData_list[1:6], function(x) reshape2::melt(x)), 
                       idcol = "method", use.names = TRUE)

count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
name.sep <- "."
count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample_id")

colData.dt <- as.data.table(count.data.merge)[, .N, by = .(variable, method, 
                                                           group,
                                                           sample.name.ext)][, N := NULL][]
colData <- data.frame(colData.dt, 
                      row.names = "sample.name.ext")

# glmm.wide = dcast(glmm.melt, method+Var2~Var1, value.var = "value", fill=0)
# colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
# head(glmm.wide[,c(1:8)])

nreps = 30

# count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
# name.sep <- "."
# count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
# count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample")

# count.matrix.glmm.dt <- dcast(data = as.data.table(count.data.merge), 
#                          formula = circ_id ~ sample.name.ext, 
#                          fill = 0, fun.aggregate = sum, 
#                          value.var = "value")
# count.matrix.filtered.glmm.dt <- as.data.table(count.matrix.glmm.dt)#[rowSums(as.data.table(count.matrix.glmm.dt)[,-"circ_id"]>=2)>1,]
# count.matrix.glmm <- as.matrix(count.matrix.filtered.glmm.dt, 
#                           rownames = "circ_id")[, rownames(colData)]

## DM1
count.matrix.glmm <- CREART::get_combined_matrix("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/DM1/", 
                                         select_methods = unique(colData$method))
## IPF
count.matrix.glmm <- CREART::get_combined_matrix("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/realdata/IPF/analyses/", 
                                                 select_methods = unique(colData$method))

glmm.long = melt(count.matrix.glmm)
glmm.long$method = sub(".*\\.","",glmm.long$X2)
glmm.long$sample = sub("\\..*", "", glmm.long$X2)
glmm.wide = dcast(glmm.long, method+sample~X1, value.var = "value", fill=0)
colnames(glmm.wide) = c("MethodID", "SampleID", colnames(glmm.wide)[-c(1,2)])
head(glmm.wide[,c(1:8)])
se.glmm <- SummarizedExperiment(assays = list(counts = count.matrix.glmm),
                           colData = colData)
eset.glmm <- ExpressionSet(assay(se.glmm),
                      AnnotatedDataFrame(as.data.frame(colData(se.glmm))))

#------------------------------------------------
## load function to perform DE analysis
#------------------------------------------------

source("/blackhole/alessia/CircModel/DEscripts.R")

#-----------------------------------------------
## useful function for zinb-wave weights 
#-----------------------------------------------

computeExactWeights <- function (model, x) 
{
  mu <- getMu(model)
  pi <- getPi(model)
  theta <- getTheta(model)
  theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
  nb_part <- dnbinom(t(x), size = theta, mu = mu)
  zinb_part <- pi * ( t(x) == 0 ) + (1 - pi) *  nb_part
  zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
  zinbwg <- t(zinbwg)
  zinbwg[x > 0] <- 1
  zinbwg[zinbwg < 1e-15] <- 1e-15
  rownames(zinbwg) = rownames(x)
  zinbwg
}

#---------------------------------------
## specify the algorithms to be compared
#---------------------------------------

algos <- list("DESeq2"= runDESeq2,
              "DESeq2-ZI"=runDESeq2.ZI,
              "DESEq2-glmGamPoi"= runDESeq2_gampoi,
              # "DESeq2-apeglm"=runDESeq2.apeglm, 
              "DESeq2-ZINBWave"= DESeq_zinbweights,
              "edgeR"=runEdgeR,
              "edgeR-robust"=runEdgeRRobust,
              "voom" = runVoom,
              "edgeR-ZINBWave"= edgeR_zinbweights,
              "circMeta" = runPois.ztest
              ) #"EBSeq"=runEBSeq)


namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

#---------------------------------
## simulation results
#---------------------------------

resTest <- list()
resHeldout <- list()
lfcTest <- list()
lfcHeldout <- list()
ratemodel <- list()

set.seed(12388)

# ---------------------------------------------------------------------------------------------------------------------
## run benchmark of NB and ZINB models across evaluation and verification datasets
# ---------------------------------------------------------------------------------------------------------------------

## multi-machines
# hosts <- c("anhal", "anhal")
# param <- SnowParam(workers = hosts, type = "SOCK")
# FUN <- function(i) system("hostname", intern=TRUE)
# bplapply(1:4, FUN, BPPARAM = param)

## single machine
# library(RNAseqData.HNRNPC.bam.chr14)
# fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
# library(GenomicAlignments) ## for GenomicRanges and readGAlignments()
# gr <- GRanges("chr14", IRanges((1000:3999)*5000, width=1000))
# param <- ScanBamParam(which=range(gr))
# FUN <- function(fl, param) {
#   gal <- readGAlignments(fl, param = param)
#   sum(countOverlaps(gr, gal))
# }
# bplapply(fls[1:3], FUN, BPPARAM = MulticoreParam(), param = param)

resMethods <- bplapply(1:nreps, function(i) {   
  
  # i = 1
  cat(i," ")
  x = e$ccp2
  testSet <- as.character(randomSubsets[i,c(1:6)])
  heldOutSet <- as.character(randomSubsets[i,-c(1:6)])
  eTest <- x[,testSet]
  keep = which(rowSums(exprs(eTest))==0)
  eTest <- eTest[-keep,]
  eTest <- eTest[which(rowSums(exprs(eTest)) > 2),]
  
  eHeldout <- x[,heldOutSet]
  keep = which(rowSums(exprs(eHeldout))==0)
  eHeldout <- eHeldout[-keep,]
  eHeldout <- eHeldout[which(rowSums(exprs(eHeldout)) > 2),]
  
  testSetGLMM <- coldata[coldata$sample%in%as.character(randomSubsets[i,1:6]),]
  heldOutSetGLMM <- coldata[coldata$sample%in%as.character(randomSubsets[i,-c(1:6)]),]
  
  zinbmodelTest <- zinbwave::zinbFit(Y = round(exprs(eTest)),
                           X = model.matrix(~ pData(eTest)$condition), K = 0,
                           epsilon = 1e10, commondispersion = FALSE, verbose = FALSE)
  weightsTest <- computeExactWeights(model = zinbmodelTest, x = exprs(eTest))
  
  zinbmodelHeldout <- zinbwave::zinbFit(Y = round(exprs(eHeldout)),
                              X = model.matrix(~ pData(eHeldout)$condition), K = 0,
                              epsilon = 1e10, commondispersion = TRUE, verbose = FALSE)
  weightsHeldout <- computeExactWeights(model = zinbmodelHeldout, x = exprs(eHeldout))
  
  resTest0 <- lapply(namesAlgos, function(n) algos[[n]](e=eTest, w=NULL))
  resHeldout0 <- lapply(namesAlgos, function(n) algos[[n]](e=eHeldout, w=NULL))
  resIdx <- rownames(eHeldout)[rownames(eHeldout)%in%rownames(eTest)]
  resTest <- as.data.frame(c(lapply(resTest0, function(z) z$padj[resIdx])))
  resHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$padj[resIdx])))
  lfcTest <- as.data.frame(c(lapply(resTest0, function(z) z$beta[resIdx])))
  lfcHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$beta[resIdx])))
  rownames(resTest) <- resIdx
  rownames(resHeldout) <- resIdx
  rownames(lfcTest) <- resIdx
  rownames(lfcHeldout) <- resIdx
  
  list(resTest=resTest,resHeldout=resHeldout,lfcTest=lfcTest,lfcHeldout=lfcHeldout)
}, BPPARAM =  MulticoreParam())

resTes <- lapply(resMethods, "[[", "resTest") #lapply(res, function(x) lapply(x, "[[", "resTest"))
resHeldout <- lapply(resMethods, "[[", "resHeldout") #lapply(res, function(x) lapply(x, "[[", "resHeldout"))
lfcTest <- lapply(resMethods, "[[", "lfcTest") #lapply(res, function(x) lapply(x, "[[", "lfcTest"))
lfcHeldout<- lapply(resMethods, "[[", "lfcHeldout") #lapply(res, function(x) lapply(x, "[[", "lfcHeldout"))

save(resTes,resHeldout,lfcTest,lfcHeldout,
     namesAlgos,
     file="/blackhole/alessia/CircModel/power/DM1_sensitivityPrecision_CCP2_NBZINBmodels.RData")

# ----------------------------------------------------------------------------------------------------------------------------------------
## run benchmark of NB vs GLMM models across evaluation and verification datasets (using unbalanced verification and evalutation datasets)
# ----------------------------------------------------------------------------------------------------------------------------------------

resTest <- list()
resHeldout <- list()
lfcTest <- list()
lfcHeldout <- list()
ratemodel <- list()
resTestGLMM_NB <- list()
resHeldoutGLMM_NB <- list()
lfcTestGLMM_NB <- list()
lfcHeldoutGLMM_NB <- list()
resTestGLMM_ZINB <- list()
resHeldoutGLMM_ZINB <- list()
lfcTestGLMM_ZINB <- list()
lfcHeldoutGLMM_ZINB <- list()

set.seed(12388)
#library("future.apply")

#plan(multisession, workers = 3)

# res <- future_lapply(X = e, future.seed = T, FUN =  function(x) { 
#   for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }
#   library(dplyr)
#   library(data.table)
#   library(plyr)

## define BIOCPARALLEL parameters
hosts <- c("grigri", "grigri", "anhal")
param <- SnowParam(workers = hosts, type = "SOCK")
param

## summary stats of simulated datasets
stats_sim <- function(matrix){
  
  # matrix : count.matrices for one verification and evaluation
    
    n.circ <- nrow(matrix)
    zero <- sum(matrix==0)/length(matrix)
    summarydata <- data.frame(zero,
                              n.circ)

  names(summarydata) <- c("perc.zeros", "n.circular")

  return(summarydata)
}

# ccp2 = as.data.frame(ccp2)
# ccp2$circ_id = rownames(ccp2)

resMethods <- bplapply(1:30, function(i) {   
    
    cat(i," ")
    # i = 1

    testSet <- as.character(randomSubsets[i,c(1:7)])
    heldOutSet <- as.character(randomSubsets[i,-c(1:7)])
    
    eTest <- e[,testSet]
    eTest_filt <- CREART::smallest_group_filter(x = as.data.table(ccp2[,c("circ_id", testSet)]), 
                                                cond = as.data.table(coldata[coldata$sample_id%in%testSet,]),
                                                rthr = 1)
    
    keep = eTest_filt$circ_id
    eTest <- eTest[keep,]
    
    summary_stat_test <- stats_sim(matrix = as.matrix(exprs(eTest)))
    
    eHeldout <- e[,heldOutSet]
    eHeldout_filt <- CREART::smallest_group_filter(x = as.data.table(ccp2[,c("circ_id", heldOutSet)]), 
                                                   cond = as.data.table(coldata[coldata$sample_id%in%heldOutSet,]),
                                                   rthr = 1)
    keep = eHeldout_filt$circ_id
    eHeldout <- eHeldout[keep,]
    
    summary_stat_heldout <- stats_sim(as.matrix(exprs(eHeldout)))
    
    testSetGLMM <- coldata[coldata$sample_id%in%as.character(randomSubsets[i,1:7]),]
    heldOutSetGLMM <- coldata[coldata$sample_id%in%as.character(randomSubsets[i,-c(1:7)]),]
    
    ## GLMM
    ## glmm-NB
    #test data
    pheno <- colData %>% dplyr::filter(variable%in%testSetGLMM$sample_id) %>% 
      dplyr::rename(SampleID = variable) %>% 
      dplyr::rename(MethodID = method) %>% 
      dplyr::rename(condition = group)

    circularcounts <- count.matrix.glmm[rownames(count.matrix.glmm)%in%rownames(eTest), rownames(pheno)]
    colnames(circularcounts) = sub("\\..*", "", colnames(circularcounts))
    testIDx = rownames(circularcounts)
    dge <- edgeR::DGEList(circularcounts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    offsets <- dge$samples$norm.factors # norm.factors per samples
    pheno$condition <- as.factor(pheno$condition)
    pheno$MethodID <- as.factor(pheno$MethodID)
    pheno$ln.lib.size <- offsets
    # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
    circularcounts <- t(circularcounts)
    allsamples_test <- cbind(pheno,circularcounts) %>% as.data.frame()

    
    # verification data
    pheno <- colData %>% dplyr::filter(variable%in%heldOutSetGLMM$sample_id) %>% 
      dplyr::rename(SampleID = variable) %>% 
      dplyr::rename(MethodID = method) %>% 
      dplyr::rename(condition = group)
    
    circularcounts <- count.matrix.glmm[rownames(count.matrix.glmm)%in%rownames(eHeldout), rownames(pheno)]
    colnames(circularcounts) = sub("\\..*", "", colnames(circularcounts))
    heldIDx = rownames(circularcounts)
    dge <- edgeR::DGEList(circularcounts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    offsets <- dge$samples$norm.factors # norm.factors per samples
    pheno$condition <- as.factor(pheno$condition)
    pheno$MethodID <- as.factor(pheno$MethodID)
    pheno$ln.lib.size <- offsets
    # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
    circularcounts <- t(circularcounts)
    allsamples_ver <- cbind(pheno,circularcounts) %>% as.data.frame()
    
    
    resIDx = intersect(testIDx, heldIDx)
    allsamples_test = allsamples_test[, colnames(allsamples_test)%in%c(colnames(allsamples_test)[1:4], resIDx)]
    allsamples_ver = allsamples_ver[, colnames(allsamples_ver)%in%c(colnames(allsamples_ver)[1:4], resIDx)]
    
    ## test GLMM-NB
    fit_GLMM_NB_test <- lapply(5:ncol(allsamples_test),
                               function(x){glmmTMB::glmmTMB(allsamples_test[,x] ~ condition + (1 | SampleID),
                                                            data=allsamples_test,
                                                            family=glmmTMB::nbinom2,
                                                            ziformula= ~0)})
    summaries <- lapply(fit_GLMM_NB_test, summary)
    pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionnormal",1]})))
    resTestGLMM_NB = data.frame(GLMM_NB = pvalues,
                                row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
    lfcTestGLMM_NB <- data.frame(GLMM_NB = lfc, row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
    
    ## verification
    fit_GLMM_NB_ver <- lapply(5:ncol(allsamples_ver),
                          function(x){glmmTMB::glmmTMB(allsamples_ver[,x] ~ condition + (1 | SampleID),
                                                       data=allsamples_ver,
                                                       family=glmmTMB::nbinom2,
                                                       ziformula= ~0)})
    summaries <- lapply(fit_GLMM_NB_ver, summary)
    pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionnormal",1]})))
    resHeldoutGLMM_NB = data.frame(GLMM_NB = pvalues,
                            row.names = colnames(allsamples_ver)[5:ncol(allsamples_ver)])
    lfcHeldoutGLMM_NB <- data.frame(GLMM_NB = lfc, row.names = colnames(allsamples_ver)[5:ncol(allsamples_ver)])

    
    resHeldoutGLMM_NB = data.frame(resHeldoutGLMM_NB[resIDx, ], row.names = resIDx)
    resTestGLMM_NB = data.frame(resTestGLMM_NB[resIDx, ], row.names = resIDx)
    lfcHeldoutGLMM_NB = data.frame(lfcHeldoutGLMM_NB[resIDx, ], row.names = resIDx)
    lfcTestGLMM_NB = data.frame(lfcTestGLMM_NB[resIDx, ], row.names = resIDx)
    
    ## glmm-ZINB
    #test GLMM-ZINB
    # fit_GLMM_ZINB_test <- lapply(5:ncol(allsamples_test),
                               # function(x){glmmTMB::glmmTMB(allsamples_test[,x] ~ condition + (1 | SampleID),
                               #                              data=allsamples_test,
                               #                              family=glmmTMB::nbinom2,
                               #                              control = glmmTMB::glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")),
                               #                              ziformula= ~1)})
    # summaries <- lapply(fit_GLMM_ZINB_test, summary)
    # pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    # pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    # lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionnormal",1]})))
    
    # resTestGLMM_ZINB = data.frame(GLMM_ZINB = pvalues,
    #                             row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
    # lfcTestGLMM_ZINB <- data.frame(GLMM_ZINB = lfc, row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
    
    #verification GLMM-ZINB
    # fit_GLMM_ZINB_ver <- lapply(5:ncol(allsamples_ver),
    #                           function(x){glmmTMB::glmmTMB(allsamples_ver[,x] ~ condition + (1 | SampleID),
    #                                                        data=allsamples_ver,
    #                                                        family=glmmTMB::nbinom2,
    #                                                        ziformula= ~1)})
    # summaries <- lapply(fit_GLMM_ZINB_ver, summary)
    # pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    # pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    # lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionnormal",1]})))
    
    # resHeldoutGLMM_ZINB = data.frame(GLMM_ZINB = pvalues,
    #                                row.names = colnames(allsamples_ver)[5:ncol(allsamples_ver)])
    # lfcHeldoutGLMM_ZINB <- data.frame(GLMM_ZINB = lfc, row.names = colnames(allsamples_ver)[5:ncol(allsamples_ver)])
    
    # resHeldoutGLMM_ZINB = data.frame(resHeldoutGLMM_ZINB[resIDx, ], row.names = resIDx)
    # resTestGLMM_ZINB = data.frame(resTestGLMM_ZINB[resIDx, ], row.names = resIDx)
    # lfcHeldoutGLMM_ZINB = data.frame(lfcHeldoutGLMM_ZINB[resIDx, ], row.names = resIDx)
    # lfcTestGLMM_ZINB = data.frame(lfcTestGLMM_ZINB[resIDx, ], row.names = resIDx)
    
    
    eTest = eTest[resIDx,]
    eHeldout = eHeldout[resIDx,]
    
    ## run other DE methods for comparison
    zinbmodelTest <- zinbFit(Y = round(exprs(eTest)),
                             X = model.matrix(~ pData(eTest)$condition), K = 0,  maxiter.optimize = 2,
                             epsilon = 1000, commondispersion = FALSE, verbose = T, BPPARAM=BiocParallel::SerialParam())
    # zinbmodelTest <- newWave(round(exprs(eTest)), 
    #                          X = model.matrix(~ pData(eTest)$condition), 
    #                          K = 0, epsilon = 1e10,)
    
    weightsTest <- computeExactWeights(model = zinbmodelTest, x = exprs(eTest))

    zinbmodelHeldout <- zinbFit(Y = round(exprs(eHeldout)),  maxiter.optimize = 2,
                                X = model.matrix(~ pData(eHeldout)$condition), K = 0,
                                epsilon = 1000, commondispersion = TRUE, verbose = T, BPPARAM=BiocParallel::SerialParam())
    # zinbmodelHeldout <- newWave(round(exprs(eHeldout)), 
    #                          X = model.matrix(~ pData(eTest)$condition), 
    #                          K = 0, epsilon = 1e10,)
    weightsHeldout <- computeExactWeights(model = zinbmodelHeldout, x = exprs(eHeldout))

    resTest0 <- lapply(namesAlgos, function(n) algos[[n]](e=eTest, w=weightsTest))
    resTest0$circMeta = runPois.ztest(e = eTest)
    resTest0$circMeta = runVoom(e = eTest)
    
    resHeldout0 <- lapply(namesAlgos, function(n) algos[[n]](e=eHeldout, w=weightsHeldout))
    resHeldout0$circMeta = runPois.ztest(e = eHeldout)
    resHeldout0$circMeta = runVoom(e = eHeldout)
    
    resTest <- as.data.frame(c(lapply(resTest0, function(z) z$padj[resIDx])))
    resHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$padj[resIDx])))
    lfcTest <- as.data.frame(c(lapply(resTest0, function(z) z$beta[resIDx])))
    lfcHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$beta[resIDx])))
    rownames(resTest) <- resIDx
    rownames(resHeldout) <- resIDx
    rownames(lfcTest) <- resIDx
    rownames(lfcHeldout) <- resIDx
    resTest$GLMM_NB = resTestGLMM_NB$resTestGLMM_NB.resIDx...[match(rownames(resTest), rownames(resTestGLMM_NB))]
    resHeldout$GLMM_NB = resHeldoutGLMM_NB$resHeldoutGLMM_NB.resIDx...[match(rownames(resHeldout), rownames(resHeldoutGLMM_NB))]
    lfcTest$GLMM_NB = lfcTestGLMM_NB$lfcTestGLMM_NB.resIDx...[match(rownames(lfcTest), rownames(lfcTestGLMM_NB))]
    lfcHeldout$GLMM_NB = lfcHeldoutGLMM_NB$lfcHeldoutGLMM_NB.resIDx...[match(rownames(lfcHeldout), rownames(lfcHeldoutGLMM_NB))]
    
    # resTest$GLMM_ZINB = resTestGLMM_ZINB$resTestGLMM_ZINB.resIDx...[match(rownames(resTest), rownames(resTestGLMM_ZINB))]
    # resHeldout$GLMM_ZINB = resHeldoutGLMM_ZINB$resHeldoutGLMM_ZINB.resIDx...[match(rownames(resHeldout), rownames(resHeldoutGLMM_ZINB))]
    # lfcTest$GLMM_ZINB = lfcTestGLMM_ZINB$lfcTestGLMM_ZINB.resIDx...[match(rownames(lfcTest), rownames(lfcTestGLMM))]
    # lfcHeldout$GLMM_ZINB = lfcHeldoutGLMM_ZINB$lfcHeldoutGLMM_ZINB.resIDx...[match(rownames(lfcHeldout), rownames(lfcHeldoutGLMM_ZINB))]
    
    list(resTest=resTest,resHeldout=resHeldout,lfcTest=lfcTest,lfcHeldout=lfcHeldout,
         weightTest = weightsTest, weightsHeldout = weightsHeldout,
         resTestGLMM=resTestGLMM_NB,resHeldoutGLMM=resHeldoutGLMM_NB,
         lfcTestGLMM=lfcTestGLMM_NB,lfcHeldoutGLMM=lfcHeldoutGLMM_NB,
         summary_stat_heldout=summary_stat_heldout,summary_stat_test=summary_stat_test
         # resTestGLMM_ZINB=resTestGLMM_ZINBNB,resHeldoutGLMM_ZINB=resHeldoutGLMM_ZINB,lfcTestGLMM_ZINB=lfcTestGLMM_ZINB,
         # lfcHeldoutGLMM_ZINB=lfcHeldoutGLMM_ZINB
         )
  }, BPPARAM =  MulticoreParam(workers = 5))
  # return(list(res=resMethods))
  # resMethods
# })

# resTesGLMMp <-lapply(resMethods, "[[", "resTestGLMM") #lapply(res, function(x) lapply(x, "[[", "resTestGLMM"))
# resHeldoutGLMMp <- lapply(resMethods, "[[", "resHeldoutGLMM") #lapply(res, function(x) lapply(x, "[[", "resHeldoutGLMM"))
# lfcTestGLMMp <- lapply(resMethods, "[[", "lfcTestGLMM") #lapply(res, function(x) lapply(x, "[[", "lfcTestGLMM"))
# lfcHeldoutGLMMp <- lapply(resMethods, "[[", "lfcHeldoutGLMM") #lapply(res, function(x) lapply(x, "[[", "lfcHeldoutGLMM"))
resTes <- lapply(resMethods, "[[", "resTest") #lapply(res, function(x) lapply(x, "[[", "resTest"))
resHeldout <- lapply(resMethods, "[[", "resHeldout") #lapply(res, function(x) lapply(x, "[[", "resHeldout"))
lfcTest <- lapply(resMethods, "[[", "lfcTest") #lapply(res, function(x) lapply(x, "[[", "lfcTest"))
lfcHeldout<- lapply(resMethods, "[[", "lfcHeldout") #lapply(res, function(x) lapply(x, "[[", "lfcHeldout"))
zinbweightTest <- lapply(resMethods, "[[", "weightTest")
zinbweightHeldout <- lapply(resMethods, "[[", "weightsHeldout")
summary_stat_test <- lapply(resMethods, "[[", "summary_stat_test")
summary_stat_heldout <- lapply(resMethods, "[[", "summary_stat_heldout")

save(resMethods, resTes, resHeldout,lfcTest,lfcHeldout,
     zinbweightTest, zinbweightHeldout,
     summary_stat_test, summary_stat_heldout,
     # resTesGLMMp,resHeldoutGLMMp,lfcTestGLMMp,lfcHeldoutGLMMp,
     namesAlgos,
     file="/blackhole/alessia/CircModel/power/IPF_sensitivityPrecision_CCP2_glmglmm_30rep.RData")

# -------------------
## type I error rate
# -------------------

randomSubsets <- read.table("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/robustness_glmm/DM1/random_shuffle50.txt",strings=FALSE)
resTest<- list()
resHeldout <- list()
lfcTest <- list()
lfcHeldout <- list()
ratemodel <- list()
nreps <- 50
set.seed(12388)
library("future.apply")

plan(multisession, workers = 2)
res <- future_lapply(X = e, future.seed = T, FUN =  function(x) { 
  for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }
  library(dplyr)
  library(data.table)
  library(plyr)
  resMethods <- bplapply(1:nreps, function(i) {   
    
    cat(i," ")
     i = 1
     x = e$ccp2
    NormSet <- as.character(randomSubsets[i,c(1:10)])
    TumorSet <- as.character(randomSubsets[i,-c(1:10)])
    data <- x[,c(NormSet, TumorSet)]
    keep = which(rowSums(exprs(data))==0)
    data <- data[-keep,]
    data <- data[which(rowSums(exprs(data)) >= 2),]
    resIdx = rownames(data)
    
    pData(data)$condition = ifelse(pData(data)$sample%in%NormSet, "B", "A")
    
    # zinbmodelTest <- zinbFit(Y = round(exprs(data)), 
    #                          X = model.matrix(~ pData(data)$condition), K = 0,
    #                          epsilon = 1e10, commondispersion = FALSE, verbose = FALSE)
    # weightsTest <- computeExactWeights(model = zinbmodelTest, x = exprs(data))
    
    resTest0 <- lapply(namesAlgos, function(n) algos[[n]](e=data, w=NULL))
    
    resTest <- as.data.frame(c(lapply(resTest0, function(z) z$padj)))
    lfcTest <- as.data.frame(c(lapply(resTest0, function(z) z$beta)))
    rownames(resTest) <- resIdx
    rownames(lfcTest) <- resIdx
    signif <- lapply(resTest0, function(x) ifelse(x$padj <= 0.1, 1, 0))
    ratemodel <- as.data.frame(c(lapply(signif, function(x) mean(x))))
    
    cat(i," ")
    glmm.wide = glmm.wide %>% dplyr::filter(SampleID%in%randomSubsets[i,]) %>% 
      mutate(condition = case_when(
        SampleID%in%randomSubsets[i, 1:10] ~ "Tumor",
        SampleID%in%randomSubsets[i, -c(1:10)] ~ "Normal") %>%
          as.factor() %>%
          structure(levels = c("Tumor","Normal")))
    pheno <- glmm.wide[,c("MethodID", "SampleID", "condition")]
    circularcounts <- as.matrix(t(glmm.wide[,-c("MethodID", "SampleID", "condition")]))
    circularcounts <- circularcounts[rownames(data), ]
    
    # circularcounts <- circularcounts[which(rowSums(circularcounts >= 4) >= 5), ]
    
    dge <- edgeR::DGEList(circularcounts)
    dge <- calcNormFactors(dge, method = "TMM")
    offsets <- dge$samples$norm.factors # norm.factors per samples
    pheno$condition <- as.factor(pheno$condition)
    pheno$MethodID <- as.factor(pheno$MethodID)
    pheno$ln.lib.size <- offsets
    # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
    circularcounts <- t(circularcounts[,pheno$SampleID])
    allsamples <- cbind(pheno,circularcounts) %>% as.data.frame()
    
    ## glmmTMB
    fitTMB <- lapply(5:ncol(allsamples),
                     function(x){glmmTMB::glmmTMB(allsamples[,x] ~ condition + (1 | SampleID),
                                                  data=allsamples,
                                                  family=nbinom2,
                                                  ziformula= ~0)})
    summaries <- lapply(fitTMB, summary)
    pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    pvalues[is.na(pvalues)] = 1
    padj = p.adjust(pvalues, method = "BH")
    
    rateTMB = list()
    rateTMB.adj = list()
    
    for (p in 1:length(pval)){
      # p=1
      signif.adj <- ifelse(padj <= pval[p], 1, 0)
      signif <- ifelse(pvalues <= pval[p], 1, 0)
      rateTMB[[p]] = mean(signif)
      rateTMB.adj[[p]] = mean(signif.adj)
      
    }
    rateGLMM = merge(rbindlist(lapply(rateTMB, function(x) data.frame(rateGLMM = x)), idcol = "pvalue"),
                 rbindlist(lapply(rateTMB.adj, function(x) data.frame(rateGLMM.adj = x)), idcol = "pvalue"), 
                 by = "pvalue")
    rateGLMM$pvalue = pval
    resTestGLMM <- as.data.frame(pvalues)
    rownames(resTestGLMM) = colnames(allsamples[,5:ncol(allsamples)])
    
    list(rate=ratemodel,resTest=resTest, resGLMM = resTestGLMM, rateGLMM = rateGLMM)
  })
  # return(list(res=resMethods))
  resMethods
})

resRate <- lapply(res, function(x) lapply(x, "[[", "rate"))
resTest <- lapply(res, function(x) lapply(x, "[[", "resTest"))
resRateGLMM <- lapply(res, function(x) lapply(x, "[[", "rateGLMM"))
resTestGLMM <- lapply(res, function(x) lapply(x, "[[", "resTestGLMM"))

save(resRate,resTest,namesAlgos,resRateGLMM,resTestGLMM,
     file="/blackhole/alessia/CircModel/data/DM1_typeIerrorGLM_GLMM_10rep.RData")

