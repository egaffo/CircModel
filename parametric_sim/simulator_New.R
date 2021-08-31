library(zinbwave)
library(edgeR)
library(EDASeq)
library(ggplot2)
library(gridExtra)
library(BiocParallel)
library(digest)
library(data.table)
source("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/reference/sampling_func_glmm.R")

## Load only one of these:

### Li
load(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/LiData_list.RData") 
load(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/LiZinbFit_models.RData") #create in datasets_and_models.R
data_list = LiData_list
models = LiData_models

### Zheng
load(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/ZhengData_list.RData") 
load(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/ZhengZinbFit_detmet_models.RData") #create in datasets_models.R
load(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/ZhengZinbFit_models.RData") #create in datasets_and_models.R
data_list = ZhengData_list
models = ZhengData_models

# # Trim by prevalence and total circRNA reads
simpleTrimGen <- function(circTab, minReads = 2, minPrev = 1) {
  # `prevalence` is the fraction of samples in which a circ is observed at
  # least `minReads` times.
  # prevalence <- rowMeans(circTab > minReads)
  absolute <- rowSums(circTab > minReads)
  # cv_index <- apply(circTab,1,function(row) stats::sd(row)/mean(row))
  ## Will only keep Circulars that appear in more than minPrev samples 
  # indCIRCs2Keep <- (prevalence > minPrev)
  indCIRCs2Keep <- (absolute > minPrev)
  return(circTab[indCIRCs2Keep, ])
}  # END - function: simpleTrim general


dataset <- names(data_list) # circRNA detection methods
# dataset <- "ccp2"
distribution <- c("ZINB", "NB")
simulation = 1:30
sampleSize = c(3,5)
TPR = c(0.1,0.5)
foldEffect <- c(2,5)
compensation <- "no"
sparsityEffect <- c(0.15, 0.30)
niter = 30

simulation_flow <- data.frame(expand.grid(simulation = simulation,
                                          dataset = dataset,
                                          distribution = distribution,
                                          sampleSize = sampleSize,
                                          TPR = TPR,
                                          foldEffect = foldEffect, 
                                          compensation = compensation,
                                          sparsityEffect = sparsityEffect
                                          ))
# Removing senseless simulations: 
# 1. nb generated dataset with sparsity_effect!=0
# simulation_flow <- simulation_flow[-which(simulation_flow$sparsityEffect != 0 & simulation_flow$distribution == "NB"),]
rownames(simulation_flow) <- 1:nrow(simulation_flow)
# Unique seed for each mix of variables
for(i in 1:nrow(simulation_flow)){
  simulation_flow$seed[i] <- strtoi(paste0("0x",substr(digest::digest(simulation_flow[i,1:ncol(simulation_flow)],algo = "sha256"),start = 1,stop = 7)))
}

### Li
save(simulation_flow,file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Li_dcc_simulation_flow.RData")

### Zheng
save(simulation_flow,file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng_detmet_simulationflow.RData")

simulation_flow.glmm <- data.frame(expand.grid(simulation = simulation,
                                          distribution = distribution,
                                          sampleSize = sampleSize,
                                          TPR = TPR,
                                          foldEffect = foldEffect, 
                                          compensation = compensation,
                                          sparsityEffect = sparsityEffect
))

rownames(simulation_flow.glmm) <- 1:nrow(simulation_flow.glmm)
# Unique seed for each mix of variables
for(i in 1:nrow(simulation_flow.glmm)){
  simulation_flow.glmm$seed[i] <- strtoi(paste0("0x",substr(digest::digest(simulation_flow.glmm[i,1:ncol(simulation_flow.glmm)],algo = "sha256"),start = 1,stop = 7)))
}
save(simulation_flow.glmm,file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng_glmm_simulation_flow.RData")

sims <- apply(simulation_flow, 1, function(sim){
  
  # sim = simulation_flow[1,]
  
  simName <- paste(colnames(simulation_flow),sim[1:ncol(simulation_flow)],
                   sep = ":",
                   collapse = "_")
  cat(simName,"\n")
  set.seed(sim[[ncol(simulation_flow)]])
  # tot_zinb_mu <- sum(exp(true_model@beta_mu))
  # link function is log
  # zinb_mu_rel <- as.numeric(exp(true_model@beta_mu))/tot_zinb_mu 
  # zinb_mu <- as.numeric(true_model@beta_mu)
  ncircular = 10000
  true_model <- models[[sim[[2]]]]
  range <- min(rbindlist(lapply(models, function(x) data.frame(length(x$sample.mu)))))
  chosen <- sample(range, ncircular, replace = TRUE)
  COUNT_CIRC <- COUNT_FUN_CIRC(sample.mu = true_model$sample.mu, 
                               sample.disp = true_model$sample.disp, 
                               zi.mu = true_model$zi.mu, 
                               zi.disp = true_model$zi.disp, 
                               zi.prop = true_model$zi.prop, 
                               relative.size = true_model$relative.size, 
                               chosen = chosen)
  obj.sim <- COUNT_CIRC(ncircular = ncircular, #n.circRNAs fixed 10.000
                        pzero = sim[[8]], #SparcityEffect c(0,0.15) scenario 4 
                        nde = NULL, #n.DECs 
                        pi0 = 1-as.numeric(sim[[5]]), #1-TPR c(0.9, 0.5) scenario 3
                        m = as.numeric(sim[[4]]), #n.samples.per.group c(3, 5)
                        fc = as.numeric(sim[[6]]), #FC c(2, 5) scenario 2
                        zinb = sim[[3]], #distribution c(NB,ZINB)
                        mod.lib = 1, mod.shape = 1, 
                        up = 0.5, #symmetry
                        replace = TRUE, 
                        ingroup=which(sample.of.origin=="grp2"),
                        simName = simName)
  return(obj.sim)
})

for(i in 1:length(sims)){ sims[[i]]$number <- i }
names(sims) <- lapply(sims,function(sim) sim$name)

## Depending on which dataset you chose at the beginning
### Zheng
save(sims, file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/ZhengData_detmet_simulationsZINB.RData")


## simulate data for glmm model
# Generating simulations
simulation_flow = simulation_flow[order(simulation_flow$simulation, simulation_flow$sampleSize, simulation_flow$TPR, simulation_flow$foldEffect), ]
simulation_flow$sims.glmm <- rep(rep(seq(1,8,1), each = 8),30)
simulation_flow <- simulation_flow[order(simulation_flow$sims.glmm),]

sims.glmm <- list()
for (k in seq(1,nrow(simulation_flow)-1, by = 8)) {
  # k = 1
  range <- min(rbindlist(lapply(models, function(x) data.frame(length(x$sample.mu)))))
  chosen <- sample(range, ncircular, replace = TRUE)
  ncircular = 10000
  
  nde = (1-(1-simulation_flow[k,5]))*ncircular
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
  lfc <- log(simulation_flow[k,6])
  
  delta[de != 0] <- lfc * de[de != 0]
  
  sims.glmm[[k]] <- apply(simulation_flow[k:(k+8-1),], 1, function(sim){
  
  # sim = simulation_flow[1,]
  
  simName <- paste(colnames(simulation_flow),sim[1:ncol(simulation_flow)],
                   sep = ":",
                   collapse = "_")
  cat(simName,"\n")
  set.seed(sim[[ncol(simulation_flow)]])
  # tot_zinb_mu <- sum(exp(true_model@beta_mu))
  # link function is log
  # zinb_mu_rel <- as.numeric(exp(true_model@beta_mu))/tot_zinb_mu 
  # zinb_mu <- as.numeric(true_model@beta_mu)
  true_model <- models[[sim[[2]]]]
  
  COUNT_CIRC <- COUNT_FUN_CIRC(sample.mu = true_model$sample.mu, sample.disp = true_model$sample.disp, 
                               zi.mu = true_model$zi.mu, zi.disp = true_model$zi.disp, zi.prop = true_model$zi.prop, 
                               relative.size = true_model$relative.size, 
                               chosen = chosen, de = de, chosen.de = chosen.de, delta = delta)
  
  obj.sim <- COUNT_CIRC(ncircular = ncircular, #n.circRNAs fixed 10.000
                    pzero = sim[[8]], #SparcityEffect 
                    nde = NULL, #n.DECs 
                    pi0 = 1-as.numeric(sim[[5]]), #1-TPR c(0.9, 0.5) scenario 3
                    m = as.numeric(sim[[4]]), #n.samples.per.group c(3, 5)
                    fc = as.numeric(sim[[6]]), #FC c(2, 5) scenario 2
                    zinb = sim[[3]], #distribution c(NB,ZINB)
                    mod.lib = 1, mod.shape = 1, 
                    up = 0.5, #symmetry
                    replace = TRUE, 
                    ingroup=which(sample.of.origin=="grp2"),
                    simName = simName)
  return(obj.sim)})
}
sims.glmm2 = sims.glmm[-which(sapply(sims.glmm, is.null))]

names(sims.glmm2) <- apply(simulation_flow.glmm, 1, function(sim) paste(colnames(simulation_flow.glmm),sim[1:ncol(simulation_flow.glmm)],
      sep = ":",
      collapse = "_"))
## Depending on which dataset you chose at the beginning
### Zheng
save(sims.glmm2, file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/ZhengData_glmm_simulationsZINB.RData")
length(sims.glmm2)
