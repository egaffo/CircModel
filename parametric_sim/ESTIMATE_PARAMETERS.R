#' Given a design matrix and some counts, we estimate various parameters for simulation:
#' - average count.
#' - spread of library sizes.
#' - dispersion trend for counts (NB).
#' - dispersion trend for counts (ZINB).

by.group
refdesign
count.matrix <- as.matrix(counts)

# colnames(count.matrix) <- rownames(coldata)

#############################################
## Filter & Normalize data               ####
#############################################
library(edgeR)
by.sample <- DGEList(count.matrix)
# by.sample <- by.sample[rowMeans(by.sample$counts) >= 1L, ]
# library(scran)
# sf <- computeSumFactors(by.sample$counts)
# by.sample$samples$norm.factors <- sf/by.sample$samples$lib.size
# by.sample$samples$norm.factors <- by.sample$samples$norm.factors/exp(mean(log(by.sample$samples$norm.factors)))
# hist(log2(sf), xlab = expression(log[2]~"size factor"), cex.lab = 1.4, cex.axis = 1.2,
#                 ylab = "Number of samples", main = "", col = "grey80")


# Normalizing using deseq2 method
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = ceiling(count.matrix[,rownames(coldata)]),
                              colData = coldata,
                              design = ~ Condition)
dds <- estimateSizeFactors(dds, type = "poscounts")
sf <- sizeFactors(dds)
pdf("../../reference/result/libesizes.pdf")
hist(log2(sf), xlab = expression(log[2]~"size factor"), cex.lab = 1.4, cex.axis = 1.2,
               ylab = "Number of samples", main = "", col = "grey80")
dev.off()

# Estimate the NB dispersion using deseq2 method
dds <- estimateDispersions(dds, fitType = "local",
                           minmu = 1e-6)

# Estimate coefficient 
dds <- nbinomWaldTest(dds)
coef(dds)

# Estimate the log-overall mean
centered.off <- getOffset(by.sample)
centered.off <- centered.off - mean(centered.off)
logmeans <- mglmOneGroup(by.sample$counts, offset = centered.off, dispersion = dispersions(dds))

pdf("/../../reference/result/avecounts.pdf")
hist(logmeans/log(2), xlab = expression(log[2]~"average count"), cex.lab = 1.4,
     cex.axis = 1.2, ylab = "Number of samples", main = "", col = "grey80")
dev.off()

pdf("../../reference/result/sampledisp.pdf")
plotDispEsts(dds)
dev.off()

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

##################################
## sampling function          ####
##################################
#' circRNA Count Data Simulation from Zero Inflated Negative-Binomial Distribution
#' 
#' This function simulates count data from Zero Inflated Negative-Binomial distribution
#' for two-sample RNA-seq experiments with given mean, dispersion 
#' and fold change. 
#' A count data matrix is generated.
#' 
#' 
#' @param nGenes total number of genes, the default value is \code{1000}.
#' @param pi0 proportion of non-differentially expressed genes, 
#'            the default value is \code{0.8}.
#' @param m sample size per treatment group.
#' @param mu a vector (or scalar) of mean counts in control group 
#'           from which to simulate.
#' @param disp a vector (or scalar) of dispersion parameter 
#'             from which to simulate.
#' @param fc a vector (or scalar, or a function that takes an integer n 
#'                     and generates a vector of length n)
#'           of fold change for differentially expressed (DE) genes.  
#' @param up proportion of up-regulated genes among all DE genes, 
#'           the default value is \code{0.5}.
#' @param replace sample with or without replacement from given parameters. 
#'                See Details for more information.
#'                
#' @return \item{counts}{circRNA count data matrix}
#' @return \item{group}{treatment group vector}
#' @return \item{lambda}{mean counts for each circRNA}
#' @return \item{phi0}{dispersion parameter for each circRNA}
#' @return \item{de}{differentially expressed circRNAs indicator: 
#'                   \code{0} for non-differentially expressed, 
#'                   \code{1} for up-regulated, 
#'                   \code{-1} for down-regulated}
#' @return \item{delta}{log2 fold change for each circRNA between 
#'                      treatment group and control group}
#'                      
#'                      

#' @details If the total number of genes \code{ngenes} is larger 
#'          than length of \code{sample.mu} or \code{sample.disp}, 
#'          \code{replace} always equals \code{TRUE}.
#'          

COUNT_FUN_CIRC <- function(sample.mu = logmeans, sample.disp = dispersion, #estimated param.s for NB
                          zi.mu = zinb.mean, zi.disp = zinb.disp, zi.prop = zinb.prop, #estimated param.s for ZINB
                          relative.size){
    # ngenes <- number of genese to be simulated
    # fc <- fold change 
    # m <- sample size per treatment group
    # pi0 <- proportion of non-differentially expressed genes
    # up <- proportion of up-regulated genes among all DE genes
    # ingroup <- sample of control group
  
    relative.size <- relative.size/mean(relative.size) # mean centering.
    
    if(length(unique(c(length(sample.mu), length(sample.disp),
                       length(zi.mu), length(zi.disp), length(zi.prop)))) != 1){
      stop("vector lengths must be the same")}
    
    function(ncircular = 10000, #n.circRNAs fixed 10.000
             pzero = sim[[8]], #SparcityEffect c(0.15,0.30)
             nde = NULL, #n.DECs 
             pi0 = 1-sim[[5]], #1-TPR c(0.9, 0.5) scenario 3
             m = sim[[4]], #n.samples.per.group c(3, 5)
             fc = sim[[6]], #FC c(2, 5) scenario 2
             zinb = sim[[3]], #distribution c("NB","ZINB")
             mod.lib = 1, mod.shape = 1, 
             up = 0.5, #symmetry
             replace = TRUE, 
             ingroup=which(sample.of.origin=="A"),
             simName = NULL){
      
    if(!is.null(pi0)){
      nde = (1-pi0)*ncircular
    }
    
    counts <- matrix(0, nrow = ncircular, ncol = 2 * m)
    
    chosen <- true_means <- true_disps <- rep(0, ncircular)
    chosen <- sample(length(sample.mu), ncircular, replace = replace)
    
      if(zinb!="ZINB"){
        true_means <- exp(sample.mu[chosen])
        true_disps <- sample.disp[chosen]
        tot_zinb_mu <- sum(true_means)
        # link function is log
        zinb_mu_rel <- as.numeric(true_means)/tot_zinb_mu
      } else {
        true_means <- zi.mu[chosen]
        true_disps <- zi.disp[chosen]
        zi.prop <- zi.prop[chosen]
        tot_zinb_mu <- sum(true_means)
        # link function is log
        zinb_mu_rel <- as.numeric(true_means)/tot_zinb_mu
      }
      
    lambda <- true_means * matrix(rep(1, ncircular), ncol = m * 2, nrow = ncircular)
  
    # Adding DE
    if (nde > 0 | !is.null(pi0)) {
      # chosen.de <- sample(ncircular, nde)
      # is.up <- rep(c(TRUE, FALSE), length.out=nde)
      # fc <- sqrt(fc)
      # lambda[chosen.de[is.up],ingroup] <- lambda[chosen.de[is.up],ingroup] * fc
      # lambda[chosen.de[is.up],-ingroup] <- lambda[chosen.de[is.up],-ingroup]/fc
      # lambda[chosen.de[!is.up],ingroup] <- lambda[chosen.de[!is.up],ingroup]/fc
      # lambda[chosen.de[!is.up],-ingroup] <- lambda[chosen.de[!is.up],-ingroup]*fc
      
      
      ## expected false positives
      FP <- round(ncircular * (1-nde/ncircular))
      TP <- ncircular - FP 
      
      ## types of true positives
      TP_up <- round(TP * up)
      TP_down <- TP - TP_up 
      
      de <- c(rep(0, FP), rep(1, TP_up), rep(-1, TP_down))
      de <- de[sample.int(length(de))] ## resample
      chosen.de <- which(de!=0)
      is.de <- logical(ncircular)
      is.de[chosen.de] <- TRUE
      
      # log fold change, approximately half positive, half negative
      delta <- rep(0, ncircular)
      if (is.function(fc)) {
        lfc <- log(fc(TP))
      } else {
        lfc <- fc
      }
      
      delta[de != 0] <- lfc * de[de != 0]
      
      lambda <- matrix(true_means, ncol = 1) %*%
        matrix(rep(1, 2 * m), nrow = 1) * #intercept
        cbind(matrix(rep(1, ncircular * m), ncol = m),
              matrix(rep(exp(delta), m), ncol = m))
      ## mean of counts
      
      phi <- matrix(rep(true_disps, 2 * m), ncol = 2 * m)
      ## dispersion of counts
      
      delta <- delta / log(2)
      
    } else {
      #simulation without differntially expressed circRNAs
      chosen.de <- integer(0)
      de <- logical(ncircular)
      de[chosen.de] <- TRUE
      delta = NULL
      is.de <- logical(ncircular)
      is.de[chosen.de] <- TRUE
    }

    nsamples <- 2 * m
    
    # Expanding and adding variable library sizes.
    # expanded <- as.integer(plates)
    # cell.means <- plate.means[,expanded]
    # ncells <- length(plates)
    # chosen.lib.mult <- sample(relative.size, ncells, replace=TRUE) * mod.lib
    # cell.means <- t(t(cell.means)*chosen.lib.mult)
    
    # Simulating counts.
    counts <- matrix(rnbinom(ncircular * 2 * m, mu = lambda, 
                             size = 1/true_disps), ncol = 2 * m, nrow = ncircular)

    # Information about True Positives
    CIRC_names <- paste("circ_",1:nrow(counts),sep = "") 
    newSampleNamesUp <- paste0(CIRC_names[de == 1], "-TPup")
    newSampleNamesDown <- paste0(CIRC_names[de == -1], "-TPdown")
    CIRC_names[de == 1] <- newSampleNamesUp
    CIRC_names[de == -1] <- newSampleNamesDown
    
    ### Plots mu
    # par(mfrow = c(2,2))
    # plot(true_means, col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "up"))+1,main = "Original Up",ylim = c(-5,5))
    # plot(true_means, col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "down"))+1,main = "Original Down",ylim = c(-5,5))
    # plot(log(zinb_mu_fc),col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "up"))+1,main = "Added fold effect Up",ylim = c(-5,5))
    # plot(log(zinb_mu_fc),col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "down"))+1,main = "Added fold effect Down",ylim = c(-5,5))
    
    # if (zinb=="ZINB"){
      # if(!is.null(pzero)){
      #   
      #   zi.prop[grepl(x = CIRC_names, pattern = "up")] <- unlist(sapply(zi.prop[grepl(x = CIRC_names,
      #                                                                                 pattern = "up")],
      #                                                                   function(x) max(x-as.numeric(pzero),
      #                                                                                   1e-08)))
      #   zi.prop[grepl(x = CIRC_names, pattern = "down")] <- unlist(sapply(zi.prop[grepl(x = CIRC_names,
      #                                                                                   pattern = "down")],
      #                                                                     function(x) min(x+as.numeric(pzero),
      #                                                                                     1-1e-8)))
      #   is.zero <- matrix(rbinom(ncircular * 2 * m, 1, zi.prop), ncol = 2 * m, 
      #                     nrow = ncircular)==1L
      #   counts[is.zero] <- 0
      # 
      #   } else {
        
          # is.zero <- matrix(rbinom(ncircular * 2 * m, 1, zi.prop), ncol = 2 * m, 
          #                   nrow = ncircular)==1L
          # counts[is.zero] <- 0
        # }
      # }
    
    # Sample names
    sample_names <- c(paste("Sample_",1:as.numeric(m),"_grp1",sep=""),paste("Sample_",(as.numeric(m)+1):(2*as.numeric(m)),"_grp2",sep=""))
    colnames(counts) <- sample_names
    rownames(counts) <- CIRC_names
    # Trim too rare OTUs
    counts_filtered <- simpleTrimGen(counts)
    counts_filtered_NA <- counts_filtered[!is.na(rownames(counts_filtered)),]
    zinb_mu_rel_fc <- data.frame(apply(lambda, 1, mean))
    rownames(zinb_mu_rel_fc) <- CIRC_names
    
    obj <- list(ps = phyloseq(otu_table(counts_filtered_NA,taxa_are_rows = TRUE),
                              sample_data(data.frame(grp = rep(c("grp1","grp2"),
                                                               each = as.numeric(m)),
                                                     row.names = sample_names))),
                name = simName,
                beta = log(zinb_mu_rel_fc[rownames(counts_filtered_NA),]))
    return(obj)
    }

}

COUNT_CIRC <- COUNT_FUN_CIRC(sample.mu = logmeans, sample.disp = dispersion, 
                            zi.mu = zinb.mean, zi.disp = zinb.disp, zi.prop = zinb.prop, 
                            relative.size = by.sample$samples$lib.size
                            )
indir <- "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/results_EstParRealdata"

saveRDS(COUNT_CIRC, file = file.path(indir, "function_ZhengEst.rds"))

                            

















