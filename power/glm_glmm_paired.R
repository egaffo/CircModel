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
          "limma",
          "data.table",
          "dplyr")
for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }

# ------------------------------------------------------------------------
## load data and meta data from previous data formatting in GOF.Rmd step
# ------------------------------------------------------------------------

load("/blackhole/alessia/CircModel/data/DM1Data_list.RData")
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

#-----------------------------------------
## summarized experiment for single matrix
#-----------------------------------------

# se <- lapply(lapply(Data_list, function(x) x[,match(randomSubsets[1,],colnames(x))]), function(x) SummarizedExperiment(x, colData = coldata[rownames(coldata)%in%randomSubsets[1,],]))
se <- lapply(DM1Data_list, function(x) SummarizedExperiment(x[,rownames(coldata)], 
                                                            colData = coldata))
eset <- lapply(se, function(x) ExpressionSet(assay(x),
                                             AnnotatedDataFrame(as.data.frame(colData(x)))))

e <- lapply(eset, function (x) {
  pData(x)$condition <- factor(pData(x)$group)
  levels(pData(x)$condition) <- c("A", "B")
  x.new = x[which(rowSums(exprs(x))!=0),]
  return(x.new)})

#-----------------------------------------
## summarized experiment for GLMM matrix
#-----------------------------------------

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

nreps = 30

count.data.melt <- as.data.table(reshape2::melt(glmm.db, id.vars = c("method", "circ_id")))
name.sep <- "."
count.data.melt[, sample.name.ext := paste0(variable, name.sep, method)]
count.data.merge <- merge(count.data.melt, coldata, by.x = "variable", by.y = "sample")

count.matrix.glmm.dt <- dcast(data = as.data.table(count.data.merge), 
                         formula = circ_id ~ sample.name.ext, 
                         fill = 0, fun.aggregate = sum, 
                         value.var = "value")
count.matrix.filtered.glmm.dt <- as.data.table(count.matrix.glmm.dt)#[rowSums(as.data.table(count.matrix.glmm.dt)[,-"circ_id"]>=2)>1,]
count.matrix.glmm <- as.matrix(count.matrix.filtered.glmm.dt, 
                          rownames = "circ_id")[, rownames(colData)]

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
              "circMeta" = runPois.ztest) #"EBSeq"=runEBSeq)

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

resTest<- list()
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

resMethods <- bplapply(1:nreps, function(i) {   
    
    cat(i," ")
    # i = 1
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
    
    ## GLMM
    ## glmm-NB
    #test
    glmm.wide.test = glmm.wide %>% dplyr::filter(SampleID%in%testSetGLMM$sample) %>% 
      dplyr::mutate(condition = case_when(
        SampleID%in%testSetGLMM$sample[testSetGLMM$group=="DM1"] ~ "Tumor",
        SampleID%in%testSetGLMM$sample[testSetGLMM$group=="Normal"] ~ "Normal") %>%
          as.factor() %>%
          structure(levels = c("Tumor","Normal"))) %>% 
      relocate(condition)
    pheno <- glmm.wide.test[,c("MethodID", "SampleID", "condition")]
    circularcounts <- as.matrix(t(glmm.wide.test[,-c(1:3)]))
    circularcounts <- circularcounts[rownames(circularcounts)%in%rownames(eTest), ]
    testIDx = rownames(circularcounts)
    dge <- edgeR::DGEList(circularcounts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    offsets <- dge$samples$norm.factors # norm.factors per samples
    pheno$condition <- as.factor(pheno$condition)
    pheno$MethodID <- as.factor(pheno$MethodID)
    pheno$ln.lib.size <- offsets
    # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
    circularcounts <- t(circularcounts[,pheno$SampleID])
    allsamples_test <- cbind(pheno,circularcounts) %>% as.data.frame()
    fit_GLMM_NB_test <- lapply(5:ncol(allsamples_test),
                     function(x){glmmTMB::glmmTMB(allsamples_test[,x] ~ condition + (1 | SampleID),
                                                  data=allsamples_test,
                                                  family=glmmTMB::nbinom2,
                                                  ziformula= ~0,
                                                  control = glmmTMB::glmmTMBControl(parallel = 2))})
    summaries <- lapply(fit_GLMM_NB_test, summary)
    pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionNormal",1]})))
    resTestGLMM_NB = data.frame(GLMM_NB = pvalues,
                         row.names = colnames(allsamples)[5:ncol(allsamples_test)])
    lfcTestGLMM_NB <- data.frame(GLMM_NB = lfc, row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
    
    # verification
    glmm.wide.heldout = glmm.wide %>% dplyr::filter(SampleID%in%heldOutSetGLMM$sample) %>% 
      dplyr::mutate(condition = case_when(
        SampleID%in%heldOutSetGLMM$sample[heldOutSetGLMM$group=="DM1"] ~ "Tumor",
        SampleID%in%heldOutSetGLMM$sample[heldOutSetGLMM$group=="Normal"] ~ "Normal") %>%
          as.factor() %>%
          structure(levels = c("Tumor","Normal"))) %>% 
      relocate(condition)
    pheno <- glmm.wide.heldout[,c("MethodID", "SampleID", "condition")]
    circularcounts <- as.matrix(t(glmm.wide.heldout[,-c(1:3)]))
    circularcounts <- circularcounts[rownames(eHeldout), ]
    heldIDx = rownames(circularcounts)
    dge <- edgeR::DGEList(circularcounts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    pheno$condition <- as.factor(pheno$condition)
    pheno$MethodID <- as.factor(pheno$MethodID)
    # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
    circularcounts <- t(circularcounts[,pheno$SampleID])
    allsamples_ver <- cbind(pheno,circularcounts) %>% as.data.frame()
    fit_GLMM_NB_ver <- lapply(5:ncol(allsamples_ver),
                          function(x){glmmTMB::glmmTMB(allsamples_ver[,x] ~ condition + (1 | SampleID),
                                                       data=allsamples_ver,
                                                       family=glmmTMB::nbinom2,
                                                       ziformula= ~0,
                                                       control = glmmTMB::glmmTMBControl(parallel = 2))})
    summaries <- lapply(fit_GLMM_NB_ver, summary)
    pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionNormal",1]})))
    resHeldoutGLMM_NB = data.frame(GLMM_NB = pvalues,
                            row.names = colnames(allsamples)[5:ncol(allsamples)])
    lfcHeldoutGLMM_NB <- data.frame(GLMM_NB = lfc, row.names = colnames(allsamples)[5:ncol(allsamples)])

    resIDx = intersect(testIDx, heldIDx)
    
    resHeldoutGLMM_NB = data.frame(resHeldoutGLMM_NB[resIDx, ], row.names = resIDx)
    resTestGLMM_NB = data.frame(resTestGLMM_NB[resIDx, ], row.names = resIDx)
    lfcHeldoutGLMM_NB = data.frame(lfcHeldoutGLMM_NB[resIDx, ], row.names = resIDx)
    lfcTestGLMM_NB = data.frame(lfcTestGLMM_NB[resIDx, ], row.names = resIDx)
    
    ## glmm-ZINB
    #test
    fit_GLMM_ZINB_test <- lapply(5:ncol(allsamples_test),
                               function(x){glmmTMB::glmmTMB(allsamples[,x] ~ condition + (1 | SampleID),
                                                            data=allsamples_test,
                                                            family=glmmTMB::nbinom2,
                                                            ziformula= ~1,
                                                            control = glmmTMB::glmmTMBControl(parallel = 2))})
    summaries <- lapply(fit_GLMM_ZINB_test, summary)
    pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionNormal",1]})))
    
    resTestGLMM_ZINB = data.frame(GLMM_ZINB = pvalues,
                                row.names = colnames(allsamples)[5:ncol(allsamples_test)])
    lfcTestGLMM_ZINB <- data.frame(GLMM_ZINB = lfc, row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
    
    #verification
    fit_GLMM_ZINB_ver <- lapply(5:ncol(allsamples_ver),
                              function(x){glmmTMB::glmmTMB(allsamples_ver[,x] ~ condition + (1 | SampleID),
                                                           data=allsamples_ver,
                                                           family=glmmTMB::nbinom2,
                                                           ziformula= ~1,
                                                           control = glmmTMB::glmmTMBControl(parallel = 2))})
    summaries <- lapply(fit_GLMM_ZINB_ver, summary)
    pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
    pvalues[is.na(pvalues)] = 1
    # signif <- ifelse(pvalues < pval, 1, 0)
    # rateTMB <- mean(signif)
    lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditionNormal",1]})))
    
    resHeldoutGLMM_ZINB = data.frame(GLMM_ZINB = pvalues,
                                   row.names = colnames(allsamples)[5:ncol(allsamples)])
    lfcHeldoutGLMM_ZINB <- data.frame(GLMM_ZINB = lfc, row.names = colnames(allsamples)[5:ncol(allsamples)])
    
    resHeldoutGLMM_ZINB = data.frame(resHeldoutGLMM_ZINB[resIDx, ], row.names = resIDx)
    resTestGLMM_ZINB = data.frame(resTestGLMM_ZINB[resIDx, ], row.names = resIDx)
    lfcHeldoutGLMM_ZINB = data.frame(lfcHeldoutGLMM_ZINB[resIDx, ], row.names = resIDx)
    lfcTestGLMM_ZINB = data.frame(lfcTestGLMM_ZINB[resIDx, ], row.names = resIDx)
    
    
    ## run other DE methods for comparison
    zinbmodelTest <- zinbFit(Y = round(exprs(eTest)),
                             X = model.matrix(~ pData(eTest)$condition), K = 0,
                             epsilon = 1e10, commondispersion = FALSE, verbose = T, BPPARAM=BiocParallel::SerialParam())
    weightsTest <- computeExactWeights(model = zinbmodelTest, x = exprs(eTest))
    zinbmodelHeldout <- zinbFit(Y = round(exprs(eHeldout)),
                                X = model.matrix(~ pData(eHeldout)$condition), K = 0,
                                epsilon = 1e10, commondispersion = TRUE, verbose = T, BPPARAM=BiocParallel::SerialParam())
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
    resTest$GLMM_NB = resTestGLMM_NB$resTestGLMM_NB.resIDx...[match(rownames(resTest), rownames(resTestGLMM_NB))]
    resHeldout$GLMM_NB = resHeldoutGLMM_NB$resHeldoutGLMM_NB.resIDx...[match(rownames(resHeldout), rownames(resHeldoutGLMM_NB))]
    lfcTest$GLMM_NB = lfcTestGLMM_NB$lfcTestGLMM_NB.resIDx...[match(rownames(lfcTest), rownames(lfcTestGLMM))]
    lfcHeldout$GLMM_NB = lfcHeldoutGLMM_NB$lfcHeldoutGLMM_NB.resIDx...[match(rownames(lfcHeldout), rownames(lfcHeldoutGLMM_NB))]
    
    resTest$GLMM_ZINB = resTestGLMM_ZINB$resTestGLMM_ZINB.resIDx...[match(rownames(resTest), rownames(resTestGLMM_ZINB))]
    resHeldout$GLMM_ZINB = resHeldoutGLMM_ZINB$resHeldoutGLMM_ZINB.resIDx...[match(rownames(resHeldout), rownames(resHeldoutGLMM_ZINB))]
    lfcTest$GLMM_ZINB = lfcTestGLMM_ZINB$lfcTestGLMM_ZINB.resIDx...[match(rownames(lfcTest), rownames(lfcTestGLMM))]
    lfcHeldout$GLMM_ZINB = lfcHeldoutGLMM_ZINB$lfcHeldoutGLMM_ZINB.resIDx...[match(rownames(lfcHeldout), rownames(lfcHeldoutGLMM_ZINB))]
    
    list(resTest=resTest,resHeldout=resHeldout,lfcTest=lfcTest,lfcHeldout=lfcHeldout,
         resTestGLMM=resTestGLMM,resHeldoutGLMM=resHeldoutGLMM,lfcTestGLMM=lfcTestGLMM,lfcHeldoutGLMM=lfcHeldoutGLMM,
         resTestGLMM_ZINB=resTestGLMM_ZINBNB,resHeldoutGLMM_ZINB=resHeldoutGLMM_ZINB,lfcTestGLMM_ZINB=lfcTestGLMM_ZINB,
         lfcHeldoutGLMM_ZINB=lfcHeldoutGLMM_ZINB)
  }, BPPARAM =  MulticoreParam(workers = 3))
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

save(resTes,resHeldout,lfcTest,lfcHeldout,
     # resTesGLMMp,resHeldoutGLMMp,lfcTestGLMMp,lfcHeldoutGLMMp,
     namesAlgos,
     file="/blackhole/alessia/CircModel/power/DM1_sensitivityPrecision_CCP2_glmglmm.RData")

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

