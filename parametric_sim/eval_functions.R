# Many functions of this script are adapted from the work of:
# A broken promise: microbiome differential abundance methods do not control the false discovery rate.
# The original code is available at https://users.ugent.be/~shawinke/ABrokenPromise/index.html

# It is mandatory to use flexmix v2.3-13 and scde 1.99.1.
# https://cran.r-project.org/src/contrib/Archive/flexmix/
# https://github.com/hms-dbmi/scde/releases
# Seurat v2.3.4
# https://github.com/satijalab/seurat/releases

# install.packages("devtools")
# library(devtools)
# install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")
# download.file(url = "https://github.com/hms-dbmi/scde/archive/1.99.1.tar.gz", destfile = "1.99.1.tar.gz")
# install.packages("1.99.1.tar.gz", repos = NULL, type = 'source', dependencies = TRUE)


pkgs <- c("edgeR", 
          "limma", 
          "DESeq2",
          "phyloseq", 
          "plyr", 
          "reshape2",
          "ROCR",
          "samr",
          "apeglm",
          "ZIM",
          "zinbwave",
          "BiocParallel",
          "AUC",
          "genefilter",
          "MAST",
          "scde", # It is important to use flexmix v2.3-13 and scde 1.99.1
          "Seurat",
          "crayon",
          "MAST",
          "aod",
          "arm",
          "fdrtool",
          "lars",
          "emdist",
          "baySeq",
          "snow",
          "ShrinkBayes",
          "DEsingle",
          "VGAM",
          "Matrix",
          "maxLik",
          "MASS")
for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }


register(SerialParam())

### edgeR normalisations: TMM and RLE
normEdgeR <- function(physeq, method = c('TMM', 'RLE', 'upperquartile'))
{
  # require(edgeR)
  CircTab <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq))
  {
    otuTab <- t(otuTab)
  } else {}
  
  if (method == "upperquartile")
  {
    scaledCounts <- t(CircTab) / colSums(CircTab)
    tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
      quantile(x[x != 0], probs = .75))
    normFacts <- tmpNF/exp(mean(log(tmpNF)))
    method <- "UQ"
  } else {
    normFacts <- edgeR:::calcNormFactors(CircTab, method = method)
  }# END - ifelse: upperquartile only of non-zero counts
  #VERY IMPORTANT: multiply by library sizes and renormalize. 
  #edgeR calculates scaling factors, which still have to be multiplied by library sizes to get to the size 
  #factors of effective sequencing depth, i.e. robust estimates of the library sizes
  #normFacts = normFacts*sample_sums(physeq)
  #normFacts = normFacts/exp(mean(log(normFacts)))
  if (all(is.na(normFacts))) #Resort to proportion normalization in case of failure for all samples
  {
    normFacts = sample_sums(physeq)
  }
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts)) 
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}# END - function: normEdgeR

### function that apply different normalisations and build *DESeqDataSet* object
### for DESeq2 analysis
normDESeq2 <- function(physeq, whichOTUs = NULL, method = c("poscounts","ratio", "RLE"))
{
  # require(DESeq2)
  method <- match.arg(method)
  
  ### Coerce count data to vanilla matrix of integers and check if there are zeroes
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## select which OTUs to analyse
  # if (!missing(whichOTUs) || !is.null(whichOTUs))
  # {
  #   physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], physeq)
  # } else {}# END - if: whichOTUs
  
  #   otu_table(physeq) <- otu_table(otuTab, taxa_are_rows = TRUE)
  
  ## Calculate size factors
  if (method == "poscounts")
  {
    obj <- phyloseq_to_deseq2(physeq,design = ~grp)
    normFacts <- sizeFactors(DESeq2::estimateSizeFactors(obj,type = "poscounts"))
  } else {
    otuTab <- as(otu_table(physeq), "matrix")
    if (any(otuTab == 0))
    {
      otuTab <- otuTab + 1L
    } else {}
    normFacts <- DESeq2::estimateSizeFactorsForMatrix(otuTab)
  }
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}# END - function: normDESeq2


### Perform EdgeR, robust version for overdispersion estimation.
### edgeR_QLFTest_robust_3.6.4
#   function (counts, group, design = NULL, prior.df = 10) 
edgeR_robust <- function(physeq, 
                         design = as.formula("~ grp"), 
                         prior.df = 10, 
                         normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  # require(edgeR)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
    if (any(otu_table(physeq) == 0))
    {
      otu_table(physeq) <- otu_table(physeq) + 1L
    } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  dgeW <- estimateGLMRobustDisp(y = dge, design, prior.df = prior.df, maxit = 10)
  glmFit <- glmQLFit(y = dgeW, dispersion = dgeW$tagwise.dispersion, robust = TRUE,
                     design = design)
  glmRes <- glmQLFTest(glmFit, coef = 2)
  pval <- glmRes$table$PValue
  padj <- p.adjust(pval, "BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = taxa_names(physeq)
  statInfo <- cbind("logFC" = glmRes$table$logFC, "logCPM" = glmRes$table$logCPM, "F" = glmRes$table$F)
  list("pValMat" = pValMat, "dispEsts" = dgeW$tagwise.dispersion, "statInfo" = statInfo)
}# END: edgeRRobust

edgeR_standard <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE)
{
  # require(edgeR)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  dge <- estimateDisp(dge, design)
  glmFit <- glmFit(dge, design)
  glmRes <- glmLRT(glmFit, coef = 2)
  pval <- glmRes$table$PValue
  padj <- p.adjust(pval, "BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = taxa_names(physeq)
  statInfo <- cbind("logFC" = glmRes$table$logFC, "logCPM" = glmRes$table$logCPM, "F" = glmRes$table$F)
  list("pValMat" = pValMat, "dispEsts" = dge$tagwise.dispersion, "statInfo" = statInfo)
}# END: edgeR - Standard

edgeR_zinbweights <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE, weights)
{
  # require(edgeR)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  dge$weights <- weights
  dge = estimateDisp(dge,design)
  #plotBCV(dge)
  glmFit <- glmFit(dge, design)
  glmlrt <- glmWeightedF(glmFit,coef = 2)
  pval <- glmlrt$table$PValue
  padj <- p.adjust(pval, "BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = taxa_names(physeq)
  statInfo <- cbind("logFC" = glmlrt$table$logFC, "logCPM" = glmlrt$table$logCPM, "F" = glmlrt$table$F)
  list("pValMat" = pValMat, "dispEsts" = dge$tagwise.dispersion, "statInfo" = statInfo)
}# END: edgeR - ZINBWaVE weights

limma_voom <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  # require(edgeR)
  # require(limma)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }# END: limma-voom
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  v <- voom(dge, design, plot=FALSE, lib.size = colSums(counts)*NFs)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, n = nrow(dge), sort.by="none")
  pval <- tt$P.Value
  padj <- p.adjust(pval,method="BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(tt)
  statInfo <- cbind("logFC" = tt$logFC, "t" = tt$t)
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END: limma voom

limma_voom_zinbweights <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"), weights)
{
  # require(edgeR)
  # require(limma)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  v <- voom(dge, design, plot = FALSE, weights = weights, lib.size = colSums(counts)*NFs)
  v$weights <- v$weights * weights
  fit <- lmFit(v, design, weights = v$weights)
  fit$df.residual <- rowSums(weights) - ncol(design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, n = nrow(dge), sort.by="none")
  pval <- tt$P.Value
  padj <- p.adjust(pval,method="BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(tt)
  statInfo <- cbind("logFC" = tt$logFC, "t" = tt$t)
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END: limma voom - ZINBWaVE weights

### performs negative binomial two-sample test of *DESeq2* to detect Diff. Abund.
negBinTestDESeq2 <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                             scRNAseq = FALSE,
                             normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  if (any(otu_table(physeq) == 0))
  {
   otu_table(physeq) <- otu_table(physeq) + 1L
  } else {}
  dds <- phyloseq_to_deseq2(physeq, design = design)
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
    }
  if(normFacts %in% c("NF.TMM","NF.TSS")){
    NFs = NFs/exp(mean(log(NFs)))
    }
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### Run DESeq
  if(scRNAseq == TRUE){
    # sizeFactors(dds) <- scran::computeSumFactors(assay(dds))
    ddsRes <- DESeq(object = dds,test = "LRT", reduced = ~1, parallel = FALSE, 
                    useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
    dispEsts <- dispersions(ddsRes)
    ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  } else {
    ddsRes <- DESeq(object = dds,test = "LRT", reduced = ~1,parallel = FALSE)
    dispEsts <- dispersions(ddsRes)
    ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  }
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ddsRes$log2FoldChange, "LRT" = ddsRes$stat)
  list("pValMat" = pValMat, "dispEsts" = dispEsts, "statInfo" = statInfo)
}# END - function: negBinTestDESeq2

### performs negative binomial two-sample test of *DESeq2* to detect Diff. Abund. accounting for zero inflation through zinb weights
negBinTestDESeq2_zinbweights <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                                         normFacts = c("TMM", "RLE", "poscounts", "ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE, weights)
{
  register(SerialParam())
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  # if (any(otu_table(physeq) == 0))
  # {
  #  otu_table(physeq) <- otu_table(physeq) + 0.0001L
  # } else {}
  
  dds <- phyloseq_to_deseq2(physeq, design = design)
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  
  if(normFacts == "NF.TMM"){
    NFs = NFs *sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### ZINB-WaVE weights
  counts <- as(otu_table(physeq), "matrix")
  weights[which(weights<1e-6)] <- 1e-06
  assays(dds)[["weights"]] = weights
  ### Run DESeq
  
  ddsRes <- DESeq(object = dds, test = "LRT", reduced = ~1, parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ddsRes$log2FoldChange, "LRT" = ddsRes$stat)
  list("pValMat" = pValMat,"dispEsts" = dispEsts, "statInfo" = statInfo)
}# END - function: negBinTestDESeq2 + ZINBWaVE

### performs negative binomial two-sample test of *DESeq2-ape.glm* to detect Diff. Abund.
negBinTestDESeq2apeglm <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                             normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #if (any(otu_table(physeq) == 0))
  #{
  #  otu_table(physeq) <- otu_table(physeq) + 1L
  #} else {}
  dds <- DESeqDataSetFromMatrix(physeq@otu_table@.Data, colData=physeq@sam_data, design=design)
  # dds <- phyloseq_to_deseq2(physeq, design = design)
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### Run DESeq
  ddsRes <- DESeq(object = dds,test = "LRT", reduced = ~1,parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ape.nb <- lfcShrink(ddsRes, coef=2, type="apeglm")
  # ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ape.nb$baseMean >= IndepFilter & !is.na(ape.nb$pvalue)
    ddsResFilt <- ape.nb[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ape.nb <- as(ddsResFilt, "data.frame")
    ape.nb[order(ape.nb$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ape.nb[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ape.nb$log2FoldChange)
  list("pValMat" = pValMat, "dispEsts" = dispEsts, "statInfo" = statInfo)
}# END - function: negBinTestDESeq2 using ape.glm

### performs negative binomial two-sample test of *DESeq2-ape.glm-zinb* to detect Diff. Abund.
negBinTestDESeq2apeglmZINB <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                                   normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #if (any(otu_table(physeq) == 0))
  #{
  #  otu_table(physeq) <- otu_table(physeq) + 1L
  #} else {}
  dds <- DESeqDataSetFromMatrix(physeq@otu_table@.Data, colData=physeq@sam_data, design=design)
  # dds <- phyloseq_to_deseq2(physeq, design = design)
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### Run DESeq
  ddsRes <- DESeq(object = dds,test = "LRT", reduced = ~1,parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  
  Y <- assay(dds)
  design <- model.matrix(design(dds), data=colData(dds))
  wts <- assays(dds)[["weights"]]
  # combine dispersion and wts into a parameter matrix,
  # which will be passed to apeglm
  param <- cbind(dispEsts, 1 - wts)
  offset <- matrix(log(sizeFactors(dds)), nrow=nrow(dds),
                   ncol=ncol(dds), byrow=TRUE)
  # need to put to natural log scale for apeglm
  mle <- log(2) * cbind(ddsRes$log2FoldChange, ddsRes$lfcSE)
  logLikZINB <- function (y, x, beta, param, offset) {
    xbeta <- x %*% beta + offset
    mean.hat <- exp(xbeta)
    k <- 1/param[1]
    omega <- param[-1]
    ZIM::dzinb(y, k=k, lambda=mean.hat, omega=omega, log=TRUE)
  }
  # run apeglm with a ZINB likelihood and zinbwave weights
  # used to define the probability of an excess zero
  fit <- apeglm(Y=Y, x=design, log.lik=logLikZINB, param=param,
                coef=2, mle=mle, offset=offset)
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(c(fit$svalue))
  colnames(pValMat) <- c("adjP")
  statInfo <- cbind("logFC" = log2(exp(1)) * fit$map[,2])
  list("pValMat" = pValMat, "dispEsts" = dispEsts, "statInfo" = statInfo)
}# END - function: negBinTestDESeq2 using ape.glm with zinb distribution


Seuratmodel <- function(physeq, design = as.formula("~ grp"),
                        normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  # ## add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  groupVar <- get_variable(physeq, "grp")
  counts <- as(otu_table(physeq), "matrix")
  names(groupVar) <- colnames(counts)
  
  sobj <- CreateSeuratObject(raw.data = counts,min.cells = 1,
                             normalization.method = "LogNormalize",
                             do.scale = TRUE,do.center = TRUE)
  sobj <- AddMetaData(sobj,metadata = groupVar,col.name = "grp")
  sobj@ident <- as.factor(groupVar)
  # Gene selection for input to CCA
  sobj <- FindVariableGenes(sobj, do.plot = F)
  response <- FindMarkers(sobj, ident.1 = "grp2", ident.2 = "grp1",print.bar = FALSE)
  rawP <- response$p_val
  adjP <- p.adjust(rawP,method = "BH")
  pValMat <- as.matrix(cbind(rawP=rawP,adjP = adjP))
  rownames(pValMat) <- rownames(response)
  otu_na <- which(is.na(match(rownames(counts),rownames(pValMat))))
  if(!is.null(otu_na))
  {
    pValMat <- rbind(pValMat,matrix(NA,ncol = 2,nrow = length(otu_na),dimnames = list(c(rownames(counts)[otu_na]),c("rawP","adjP"))))
  }
  statInfo <- response[,2:4]
  return(list("pValMat" = pValMat, "statInfo" = statInfo))
}# END - function: Seurat Wilcoxon test

SigEMD <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                             normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), 
                   returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  sizeFactor <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  norm.counts <- t(t(physeq@otu_table@.Data) / sizeFactor)
  
  cond <- as.character(physeq@sam_data$grp)
  names(cond) <- colnames(norm.counts)
  databin <- databin(norm.counts)
  
  Hur_gene <- idfyImpgene(norm.counts, databin, cond)
  genes_use <- idfyUsegene(norm.counts, databin, cond, ratio = 0.8) 
  # ratio is default as 0.8, but it can be set by users to obtain appropriate number of genes_use.
  
  if(length(Hur_gene)>2){
    datimp <- FunImpute(object = norm.counts, genes_use = (genes_use), genes_fit = (Hur_gene), 
                        dcorgene = NULL) 
    data <- datimp$alldat
    results <- calculate_single(data =  data, condition =  cond, Hur_gene = Hur_gene, binSize=0.2,nperm=5) 
  } else {
    results <- calculate_single(data =  norm.counts, condition =  cond, Hur_gene = NULL, binSize=0.2, nperm=5)
  }

  # results$emdall
  pValMat <- as.matrix(results$emdall[,"pvalue"])
  colnames(pValMat) <- c("adjP")
  statInfo <- cbind("emd" = results$emdall[,"emd"])
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END - function: SigEMD

BaySeq <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                   normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), 
                   returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  sizeFactor <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  norm.counts <- t(t(physeq@otu_table@.Data) / sizeFactor)
  
  cl <- snow:::makeCluster(4, "SOCK")
  group <- as.numeric(physeq@sam_data$grp)
  cd <- new("countData", data = norm.counts, replicates = group, groups = list(NDE = rep(1, length(group)), 
                                                                               DE = group))
  # cd@libsizes <- getLibsizes(cd)
  cd <- getPriors.NB(cd, equalDispersions = TRUE, estimation = "QL", cl = cl)
  cd <- getLikelihoods.NB(cd, pET = "BIC", cl = cl)
  try(snow:::stopCluster(cl), silent = TRUE)
  cd.table <- topCounts(cd, group = "DE", number = nrow(norm.counts))
  id <- match(row.names(norm.counts), row.names(cd.table))

  pValMat <- as.matrix(cd.table[id,]$FDR)
  colnames(pValMat) <- c("adjP")
  statInfo <- NULL
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END - function: baySeq

ShrinkBayes <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                   normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), 
                   returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  myoffsets <- log(NFs)
  
  g <- as.factor(physeq@sam_data$grp)
  form <- y ~ 1 + g + offset(myoffsets)
  form0 <- y ~ 1
  shrinksimul <- ShrinkSeq(form = form, dat = physeq@otu_table@.Data, shrinkfixed = "g", fams = "zinb")
  fitall <- FitAllShrink(form, dat = norm.counts, fams = "zinb",shrinksimul = shrinksimul)
  fitall0 <- FitAllShrink(form0, dat = norm.counts, fams = "zinb",shrinksimul = shrinksimul)  
  npprior <-MixtureUpdatePrior(fitall = fitall, fitall0 = fitall0, shrinkpara="g", ncpus = mc.cores)
  nppostshr <- MixtureUpdatePosterior(fitall, npprior, fitall0)
  lfdr <- SummaryWrap(nppostshr)
  fdr <- BFDR(lfdr)
  pGlobal <- as.vector(fdr)
  rm(list = c("ptm", "g", "form", "form0", "shrinksimul", "fitall", "fitall0", "npprior" ,"nppostshr", "lfdr", "fdr"))
  
  pValMat <- as.matrix(pGlobal)
  colnames(pValMat) <- c("adjP")
  statInfo <- NULL
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END - function: ShrinkBayes


DeSingle <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                   normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), 
                   returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  sizeFactor <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  # norm.counts <- t(t(physeq@otu_table@.Data) / sizeFactor)
  
  # Detecting the DE genes
  source("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/reference/DEsingle.R")
  if(is.null(normFacts)){
    results <- DEsingle_AFB(counts = physeq@otu_table@.Data, group = as.factor(physeq@sam_data$grp), sf = NULL)
  } else {
    results <- DEsingle_AFB(counts = physeq@otu_table@.Data, group = as.factor(physeq@sam_data$grp), sf = sizeFactor)
  }
  
  # Dividing the DE genes into 3 categories at threshold of FDR < 0.05
  # results.classified <- DEtype(results = results, threshold = 1)
  
  # results$emdall
  pValMat <- as.matrix(results$pvalue)
  colnames(pValMat) <- c("adjP")
  statInfo <- cbind("foldChange" = results$foldChange) 
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END - function: DeSingle


scdemodel <- function(physeq,design = as.formula("~ grp"))
{
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  groupVar <- get_variable(physeq, "grp")
  counts <- as(otu_table(physeq), "matrix")
  names(groupVar) <- colnames(counts)
  counts <- apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
  o.ifm <- scde.error.models(counts = counts, groups = groupVar, n.cores = 4, 
                             threshold.segmentation = TRUE, save.crossfit.plots = FALSE, 
                             save.model.plots = FALSE, 
                             #verbose = 1, 
                             min.size.entries = nrow(counts)*0.2)
  # filter out cells (samples) that don't show positive correlation with
  # the expected expression magnitudes (very poor fits)
  valid.samples <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.samples, ]
  # estimate gene expression (OTU abundance) prior
  o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
  # Differential abundace analysis
  # define two groups of cells
  groups <- groupVar[valid.samples]
  # run differential expression tests on all genes.
  ediff <- scde.expression.difference(o.ifm, counts, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
  
  rawP <- 1-pnorm(ediff$Z)
  adjP <- p.adjust(rawP,method = "BH")
  pValMat <- as.matrix(cbind(rawP=rawP,adjP = adjP))
  rownames(pValMat) <- rownames(ediff)
  statInfo <- ediff[,1:4]
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END - function: Single-Cell Differential Expression Analysis scde

runGFOLD <- function(physeq) {
  e <- ExpressionSet(physeq@otu_table@.Data, AnnotatedDataFrame(physeq@sam_data))
  
  gfold.samples.dir <- "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/GFOLD"
  gfold.out.dir <- "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/GFOLD"
  m <- ncol(e)
  n <- nrow(e)
  sample.fls <- replicate(m, tempfile(pattern = "sample",
                                      tmpdir = gfold.samples.dir, fileext = ""))
  out.fl <- tempfile(pattern = "out",
                     tmpdir = gfold.out.dir, fileext = "")
  for (i in seq_len(m)) {
    out <- cbind(1:n, rep(0,n), exprs(e)[,i], rep(0,n), rep(0,n))
    write.table(out, file=sample.fls[i],
                quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  }
  samples1 <- paste(sample.fls[pData(e)$grp == "grp1"], collapse=",")
  samples2 <- paste(sample.fls[pData(e)$grp == "grp2"], collapse=",")
  system(paste0("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/gfold.V1.1.4/gfold diff -s1 ", samples1, " -s2 ", samples2, " -o ", out.fl),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  res <- read.table(out.fl)
  
  # clean up
  for (fl in c(sample.fls, out.fl, paste0(out.fl,".ext"))) {
    system(paste0("rm -f ",fl))
  }
  colnames(res) <- c("GeneSymbol","GeneName","GFOLD(0.01)",
                     "E-FDR","log2fdc","1stRPKM","2ndRPKM")
  
  pValMat <- as.matrix(cbind(rawP=res$"E-FDR",adjP = res$"E-FDR"))
  statInfo <- cbind("logFC" = res$log2fdc)
  list("pValMat" = pValMat, "statInfo" = statInfo)
}

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
  zinbwg
}

oneSimRunGSOwn <- function(physeq, beta, true_weights = NULL, epsilon = 1e10) { 
  # Prevent NA when converting to integer due to some outlier generation during simulation
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  ## all normalisations
  physeq <- normEdgeR(physeq = physeq, method = "TMM")
  physeq <- normDESeq2(physeq = physeq, method = 'RLE') 
  # physeq <- normEdgeR(physeq = physeq, method = 'upperquartile')
  physeq <- normDESeq2(physeq = physeq, method = "poscounts")  # poscounts, similar to RLE
  # physeq <- normCSS(physeq = physeq)
  # physeq <- normTSS(physeq = physeq)
  cat("Normalisations: DONE\n")
  #zinb model estimation
  zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                       X = model.matrix(~ physeq@sam_data$grp), K = 0,
                       epsilon = epsilon, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
  weights <- computeExactWeights(model = zinbmodel,x = physeq@otu_table@.Data)
  colnames(weights) <- colnames(physeq@otu_table)
  rownames(weights) <- rownames(physeq@otu_table)
  cat("ZINB model estimation: DONE\n")
  returnList = list()
  #returnList$physeq = physeq
  returnList = within(returnList, {
    Y <- physeq@otu_table
    ##true lfc
    truemodel <- list("pValMat" = NULL, "statInfo" = cbind("logFC" = beta))
    ## edgeR Robust
    edgeR_TMM_robustDisp <- edgeR_robust(physeq, normFacts = "TMM")
    cat("EdgeR robust tests: DONE\n")
    ## edgeR Standard
    edgeR_TMM_standard <- edgeR_standard(physeq, normFacts = "TMM")
    cat("EdgeR standard tests: DONE\n")
    ## edgeR with zinbwave weights
    edgeR_TMM_zinbwave <- edgeR_zinbweights(physeq, normFacts = "TMM", weights = weights)
    cat("EdgeR with ZINB-WaVE weights tests: DONE\n")
    ## edgeR with true weights
    # if(!is.null(true_weights)){
    #   edgeR_TMM_trueweights <- edgeR_zinbweights(physeq, normFacts = "TMM", weights = true_weights)
    #   cat("EdgeR with true weights tests: DONE\n")
    # }
    edgeR_poscounts_standard <- edgeR_standard(physeq, normFacts = "poscounts")
    cat("EdgeR with RLE normalisation tests: DONE\n")
    
    ## limma-voom
    limma_voom_TMM <- limma_voom(physeq, normFacts = "TMM")
    cat("Limma Voom tests: DONE\n")
    ## limma-voom zinbwave weights
    limma_voom_TMM_zinbwave <- limma_voom_zinbweights(physeq, normFacts = "TMM", weights = weights)
    cat("Limma Voom ZINB-WaVE weights tests: DONE\n")
    ## limma-voom with true weights
    # if(!is.null(true_weights)){
    #   limma_voom_TMM_trueweights <- limma_voom_zinbweights(physeq, normFacts = "TMM", weights = true_weights)
    #   cat("Limma Voom true weights tests: DONE\n")
    # }
    
    ## NB test from DESeq2
    DESeq2_poscounts <- negBinTestDESeq2(physeq, normFacts = "poscounts")
    cat("NB DESeq2 tests: DONE\n")
    ## NB test from DESeq2 with zinbwave weights
    DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",
                                                              weights = weights)
    cat("NB DESeq2 with ZINB-WaVE weights tests: DONE\n")
    ## NB test from DESeq2 with true weights
    # if(!is.null(true_weights)){
    #   DESeq2_poscounts_trueweights <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",
    #                                                                weights = true_weights)
    #   cat("NB DESeq2 with true weights tests: DONE\n")
    # }
    DESeq2_TMM <- negBinTestDESeq2(physeq, normFacts = "TMM")
    cat("NB DESeq2 with TMM normalisation tests: DONE\n")
    ## NB test from DESeq2 with ape.glm
    DESeq2_poscounts_apeglm <- negBinTestDESeq2apeglm(physeq, normFacts = "poscounts")
    cat("NB DESeq2 with ape.glm tests: DONE\n")
    ## NB test from DESeq2 with ape.glm
    # DESeq2_poscounts_apeglmzinb <- negBinTestDESeq2apeglmZINB(physeq, normFacts = "poscounts")
    # cat("NB DESeq2 with ape.glm-ZINB tests: DONE\n")
    
    # EMD
    SigEMD_TMM <- SigEMD(physeq, normFacts = "TMM")
    cat("EMD with TMM normalization tests: DONE\n")
    SigEMD_RLE <- SigEMD(physeq, normFacts = "poscounts")
    cat("EMD with RLE normalization tests: DONE\n")
    
    # Seurat
    seurat_wilcoxon <- Seuratmodel(physeq)
    cat("Seurat Wilcoxon tests: DONE\n")
    # NODES
    nodes <- NODESmodel(physeq)
    cat("NODES Wilcoxon tests: DONE\n")
    
    ## baySeq
    # baySeq_TMM <-  BaySeq(physeq, normFacts = "TMM")
    # cat("BaySeq with TMM normalization tests: DONE\n")
    # baySeq_RLE <- BaySeq(physeq, normFacts = "RLE")
    # cat("BaySeq with RLE normalization tests: DONE\n")
    
    # DEsingle
    DeSingle_TMM <- DeSingle(physeq, normFacts = "TMM")
    cat("DeSingle with TMM normalization tests: DONE\n")
    DeSingle_RLE <- DeSingle(physeq, normFacts = "poscounts")
    cat("DeSingle with RLE normalization tests: DONE\n")
    
    ## ShrinkBayes
    # ShrinkBayes_TMM <-  ShrinkBayes(physeq, normFacts = "TMM")
    # cat("ShrinkBayes with TMM normalization tests: DONE\n")
    # ShrinkBayes_RLE <- ShrinkBayes(physeq, normFacts = "RLE")
    # cat("ShrinkBayes with RLE normalization tests: DONE\n")
    
    # scde single cell differential expression
    scde <- scdemodel(physeq)
    cat("scde single cell differential expression tests: DONE\n")
  
  })
  return(returnList)
}


oneSimRunscde <- function(physeq) { 
  # Prevent NA when converting to integer due to some outlier generation during simulation
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  
  returnList = list()
  #returnList$physeq = physeq
  returnList = within(returnList, {
    ## SCDE model
    scde <- scdemodel(physeq)
    cat("scde test: DONE\n")
  })
  return(returnList)
}


oneSimRunGSOwnFastestMethod <- function(physeq, beta, true_weights = NULL, epsilon = 1e10) { 
  # Prevent NA when converting to integer due to some outlier generation during simulation
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  physeq <- normDESeq2(physeq = physeq, method = 'RLE') 
  physeq <- normDESeq2(physeq = physeq, method = "poscounts")  # poscounts, similar to RLE
  physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
  cat("Normalisations: DONE\n")
  #zinb model estimation
  zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                       X = model.matrix(~ physeq@sam_data$grp), K = 0,
                       epsilon = epsilon, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
  weights <- computeExactWeights(model = zinbmodel,x = physeq@otu_table@.Data)
  colnames(weights) <- colnames(physeq@otu_table)
  rownames(weights) <- rownames(physeq@otu_table)
  cat("ZINB model estimation: DONE\n")
  returnList = list()
  #returnList$physeq = physeq
  returnList = within(returnList, {
    Y <- physeq@otu_table
    ##true lfc
    truemodel <- list("pValMat" = NULL, "statInfo" = cbind("logFC" = beta))
    ## NB test from DESeq2 with scRNAseq paremeters
    DESeq2_scRNAseq <- negBinTestDESeq2(physeq, scRNAseq = TRUE, normFacts = "poscounts")
    cat("NB DESeq2 with scRANseq parameters tests: DONE\n")
    ## NB test from DESeq2
    DESeq2 <- negBinTestDESeq2(physeq, normFacts = "RLE")
    cat("NB DESeq2 tests: DONE\n")
    ## NB test from DESeq2 with zinbwave weights
    DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",
                                                              weights = weights)
    cat("NB DESeq2 with ZINB-WaVE weights tests: DONE\n")
    ## NB test from DESeq2 with ape.glm
    DESeq2_poscounts_apeglm <- negBinTestDESeq2apeglm(physeq, normFacts = "poscounts")
    cat("NB DESeq2 with ape.glm tests: DONE\n")
  })
  return(returnList)
}

oneSimRunGSOwn_time <- function(physeq, true_weights = NULL, epsilon = 1e10) { 
  # Prevent NA when converting to integer due to some outlier generation during simulation
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  ## all normalisations
  
  #zinb model estimation
  zinb_start <- proc.time()
  zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                       X = model.matrix(~ physeq@sam_data$grp), K = 0,
                       epsilon = epsilon, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
  weights <- computeExactWeights(model = zinbmodel,x = physeq@otu_table@.Data)
  colnames(weights) <- colnames(physeq@otu_table)
  rownames(weights) <- rownames(physeq@otu_table)
  zinb_end <- proc.time()
  zinb_time <- zinb_end-zinb_start
  
  returnList = list()
  returnList = within(returnList, {
  
    edgeR_TMM_robustDisp <- system.time(
      {## edgeR Robust
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      edgeR_TMM_robustDisp <- edgeR_robust(physeq, normFacts = "TMM")
      cat("EdgeR robust tests: DONE\n")}
    )
    
    edgeR_TMM_standard <- system.time(
      {## edgeR Standard
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      edgeR_TMM_standard <- edgeR_standard(physeq, normFacts = "TMM")
      cat("EdgeR standard tests: DONE\n")}
    )
    
    edgeR_TMM_zinbwave <- system.time(
      {## edgeR with zinbwave weights
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      edgeR_TMM_zinbwave <- edgeR_zinbweights(physeq, normFacts = "TMM", weights = weights)
      cat("EdgeR with ZINB-WaVE weights tests: DONE\n")}
    ) + zinb_time
    
    edgeR_poscounts_standard <- system.time(
      {physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
      edgeR_poscounts_standard <- edgeR_standard(physeq, normFacts = "poscounts")
      cat("EdgeR with RLE normalisation tests: DONE\n")}
    )
    
    limma_voom_TMM <- system.time(
      {## limma-voom
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      limma_voom_TMM <- limma_voom(physeq, normFacts = "TMM")
      cat("Limma Voom tests: DONE\n")}
    )
    
    limma_voom_TMM_zinbwave <- system.time(
      {## limma-voom zinbwave weights
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      limma_voom_TMM_zinbwave <- limma_voom_zinbweights(physeq, normFacts = "TMM", weights = weights)
      cat("Limma Voom ZINB-WaVE weights tests: DONE\n")}
    ) + zinb_time
    
    DESeq2_poscounts <- system.time(
      {## NB test from DESeq2
        physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
        DESeq2_poscounts <- negBinTestDESeq2(physeq, normFacts = "poscounts")
        cat("NB DESeq2 tests: DONE\n")}
    )
    
    DESeq2_TMM <- system.time(
      {## DESeq2 TMM
        physeq <- normEdgeR(physeq = physeq, method = "TMM")
        DESeq2_TMM <- negBinTestDESeq2(physeq, normFacts = "TMM")
        cat("NB DESeq2 with TMM normalisation tests: DONE\n")}
    )
    
    DESeq2_poscounts_zinbwave <- system.time(
      {physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
        ## NB test from DESeq2 with zinbwave weights
        DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",weights = weights)
        cat("NB DESeq2 with ZINB-WaVE weights tests: DONE\n")}
    ) + zinb_time
    
    mgsZig_CSS <- system.time(
      {## metagenomeSeq Zero-Inflated Gaussian
        physeq <- normCSS(physeq = physeq)
        mgsZig_CSS <- metagenomeSeqZIG(physeq, normFacts = "CSS")
        cat("MetagenomeSeq ZIG tests: DONE\n")}
    )
    
    ALDEx2 <- system.time(
      {## ALDEx2 Wilcoxon test
       ALDEx2 <- ALDEx2model(physeq)
       cat("ALDEx2 Wilcoxon test: DONE\n")}
    )
    
    MAST <- system.time(
      {## MAST hurdle models
       MAST <- MASTmodel(physeq)
       cat("MAST lrt tests: DONE\n")}
    )
    
    scde <- system.time(
      {## scde single cell differential expression
        scde <- scdemodel(physeq)
        cat("scde single cell differential expression tests: DONE\n")}
    )
    
    seurat_wilcoxon <- system.time(
      {## Seurat 
        seurat_wilcoxon <- Seuratmodel(physeq)
        cat("Seurat Wilcoxon tests: DONE\n")}
    )
  })
  return(returnList)
}

