
runDESeq2 <- function(e, retDDS=FALSE,w=NULL) {
  library(DESeq2)
  col.dat <- DataFrame(pData(e))
  dds <- DESeqDataSetFromMatrix(countData = ceiling(exprs(e)[,rownames(col.dat)]),
                                          colData = col.dat,
                                          design = ~ condition)
  dds <- DESeq(dds, quiet=TRUE) 
  res <- results(dds, independentFiltering = F)
  beta <- res$log2FoldChange
  names(beta) <- rownames(exprs(e))
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals, padj=padj, beta=beta))
  cat("NB DESeq2 test: DONE\n")
  
}

runDESeq2.ZI <- function(e, retDDS=FALSE,w=NULL) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE, sfType = "poscounts", useT = TRUE, 
               minmu = 1e-6, minReplicatesForReplace = Inf)
  res <- results(dds)
  beta <- res$log2FoldChange
  names(beta) <- rownames(exprs(e))
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals, padj=padj, beta=beta))
  cat("NB DESeq2 with scRNA-seq parameters tests: DONE\n")
}

DESeq_zinbweights <- function(e, w){
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  w[w < 0.00001] = 0.00001
  assays(dds, withDimnames = F)[["weights"]] = w
  dds <- DESeq(dds, sfType="poscounts", minmu=1e-6, minRep=Inf)
  res <- results(dds)
  beta <- res$log2FoldChange
  names(beta) <- rownames(exprs(e))
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals, padj=padj, beta=beta))
  cat("NB DESeq2 with ZINB-WaVe weights tests: DONE\n")
}# END: DESeq - ZINBWaVE weights

runDESeq2_gampoi <- function(e, w){
  library(DESeq2)
  library(glmGamPoi)
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  deg <- DESeq(dds, fitType="glmGamPoi")
  res <- results(deg)
  beta <- res$log2FoldChange
  names(beta) <- rownames(exprs(e))
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals, padj=padj, beta=beta))
  cat("NB DESeq2 with glmGamPoi: DONE\n")
}# END: DESeq2 - GammaPoisson

runDESeq2.apeglm <- function(e, retDDS=FALSE,w=NULL) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE, sfType = "poscounts", useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
  res <- lfcShrink(dds, coef=2, type="apeglm")
  beta <- res$log2FoldChange
  names(beta) <- rownames(exprs(e))
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals, padj=padj, beta=beta))
  cat("NB DESeq2 with ape-glm shrinkage tests: DONE\n")
}


runEdgeR <- function(e,w=NULL) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel, design)
  dgel <- estimateGLMTrendedDisp(dgel, design)
  dgel <- estimateGLMTagwiseDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"], predbeta10=predbeta10[,"pData(e)$conditionB"]))
  cat("NB EdgeR tests: DONE\n")
}

runEdgeRRobust <- function(e,w=NULL) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"]))
  cat("NB EdgeR-Robust tests: DONE\n")
}

edgeR_zinbweights <- function(e, w){
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel$weights <- w
  dgel <- calcNormFactors(dgel)
  
  dgel = estimateDisp(dgel,design)
  #plotBCV(dgel)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmWeightedF(edger.fit, coef = 2)
  pvals <- edger.lrt$table$PValue
  padj <- p.adjust(pvals, "BH")
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))

  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)

  return(list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"]))
  cat("NB EdgeR with ZINB-WaVe weights tests: DONE\n")
}# END: edgeR - ZINBWaVE weights

runVoom <- function(e,w) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  beta = tt$logFC
  names(beta) <- rownames(exprs(e))
  names(beta) <- rownames(beta)
  pvals <- tt$P.Value 
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  list(pvals=pvals, padj=padj, beta=beta)
  cat("Limma + Voom tests: DONE\n")
}# END: limma+voom

runSAMseq <- function(e,w) {
  set.seed(1)
  x <- exprs(e)
  y <- pData(e)$condition
  capture.output({samfit <- SAMseq(x, y, resp.type = "Two class unpaired")})
  beta <- log2(samfit$samr.obj$foldchange)
  names(beta) <- rownames(exprs(e))
  pvals <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  names(pvals) <- rownames(exprs(e))
  padj[is.na(padj)] <- 1
  names(padj) <- rownames(exprs(e))
  return(list(pvals=pvals,padj=padj,beta=beta))
  # cat("SAMseq tests: DONE\n")
  
}

runEBSeq <- function(e,w) {
  sizes <- MedianNorm(exprs(e))
  out <- capture.output({
    suppressMessages({
      res <- EBTest(Data = exprs(e),
                    Conditions = pData(e)$condition,
                    sizeFactors = sizes,
                    maxround = 5)
    })
  })
  padj <- rep(1, nrow(exprs(e)))
  # we use 1 - PPDE for the FDR cutoff as this is recommended in the EBSeq vignette
  padj[match(rownames(res$PPMat), rownames(e))] <- res$PPMat[,"PPEE"]
  names(padj) <- rownames(exprs(e))
  beta <- rep(0, nrow(exprs(e)))
  names(beta) <- rownames(exprs(e))
  return(list(pvals=padj, padj=padj, beta=beta))
  # cat("EBSeq tests: DONE\n")
}

# runSCDE <- function(e,w){
#   groupVar <- pData(e)$condition
#   counts <- exprs(e)
#   names(groupVar) <- colnames(counts)
#   counts <- apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
#   o.ifm <- scde.error.models(counts = counts, groups = groupVar, n.cores = 4, 
#                              threshold.segmentation = TRUE, save.crossfit.plots = FALSE, 
#                              save.model.plots = FALSE, 
#                              #verbose = 1, 
#                              min.size.entries = nrow(counts)*0.2)
#   # filter out cells (samples) that don't show positive correlation with
#   # the expected expression magnitudes (very poor fits)
#   valid.samples <- o.ifm$corr.a > 0
#   o.ifm <- o.ifm[valid.samples, ]
#   # estimate gene expression (OTU abundance) prior
#   o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
#   # Differential abundace analysis
#   # define two groups of cells
#   groups <- groupVar[valid.samples]
#   # run differential expression tests on all genes.
#   ediff <- scde.expression.difference(o.ifm, counts, o.prior, groups  =  groups, 
#                                       n.randomizations  =  100, n.cores  =  1, verbose  =  1)
#   rawP <- 1-pnorm(ediff$Z)
#   adjP <- p.adjust(rawP,method = "BH")
#   list(pvals=rawP, padj=adjP, beta=beta)
# 
# }
