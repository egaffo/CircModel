library(zinbwave)
library(edgeR)
library(EDASeq)
library(ggplot2)
library(gridExtra)
library(BiocParallel)
library(digest)
library(data.table)
source("/blackhole/alessia/CircModel/parametric_sim/sampling_func_glmm.R")

## Load only one of these:
### IPF
load(file = "/blackhole/alessia/CircModel/data/IPFData_list.RData") 
load(file = "/blackhole/alessia/CircModel/data/IPFZinb_nb_Fit_detmet_models.RData") #created in datasets_and_models.R
data_list = IPFData_list
models = IPFData_models

### ALZ
load(file = "/blackhole/alessia/CircModel/data/ALZData_list.RData") 
load(file = "/blackhole/alessia/CircModel/data/ALZZinb_nb_Fit_detmet_models.RData") #created in datasets_and_models.R
load("/blackhole/alessia/CircModel/data/ALZData_list.RData")
data_list = ALZDataFilt_list
models = ALZData_models

# # Trim by prevalence and total circRNA reads
simpleTrimCirc <- function(circTab, minReads = 2, minPrev = 1) {
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


# dataset <- names(data_list) # circRNA detection methods
dataset <- "ccp2"
distribution <- c("ZINB")
simulation = 1:30
sampleSize = c(3,10,20)
TPR = c(0.1)
foldEffect <- c(1.5)
# sparsityEffect <- c(0.15, 0.30)
niter = 30

simulation_flow <- data.frame(expand.grid(simulation = simulation,
                                          dataset = dataset,
                                          distribution = distribution,
                                          sampleSize = sampleSize,
                                          TPR = TPR,
                                          foldEffect = foldEffect 
                                          # sparsityEffect = sparsityEffect
                                          ))
# Removing senseless simulations: 
# 1. nb generated dataset with sparsity_effect!=0
# simulation_flow <- simulation_flow[-which(simulation_flow$sparsityEffect != 0 & simulation_flow$distribution == "NB"),]
rownames(simulation_flow) <- 1:nrow(simulation_flow)
dim(simulation_flow)
# Unique seed for each mix of variables
for(i in 1:nrow(simulation_flow)){
  simulation_flow$seed[i] <- strtoi(paste0("0x",substr(digest::digest(simulation_flow[i,1:ncol(simulation_flow)],algo = "sha256"),
                                                       start = 1, stop = 7)))
}

### IPF
save(simulation_flow,file = "/blackhole/alessia/CircModel/data/IPF_simulation_flow.RData")
### ALZ
save(simulation_flow,file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_simulation_flow.RData")

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
  ncircular = 5000
  true_model <- models[[sim[[2]]]]
  range <- min(rbindlist(lapply(models, function(x) data.frame(length(x$sample.mu)))))
  chosen <- sample(range, ncircular, replace = TRUE)
  COUNT_CIRC <- COUNT_FUN_CIRC(sample.mu = true_model$sample.mu, 
                               sample.disp = true_model$sample.disp, 
                               zi.mu = true_model$zi.mu, 
                               zi.disp = true_model$zi.disp, 
                               zi.prop = true_model$zi.prop, 
                               relative.size = true_model$relative.size)
  obj.sim <- COUNT_CIRC(ncircular = ncircular, #n.circRNAs fixed 10.000
                        #pzero = sim[[7]], #SparcityEffect c(0,0.15) scenario 4 
                        nde = NULL, #n.DECs 
                        pi0 = 1-as.numeric(sim[[5]]), #1-TPR c(0.9, 0.5) scenario 3
                        m = as.numeric(sim[[4]]), #n.samples.per.group c(3, 5)
                        fc = as.numeric(sim[[6]]), #FC c(2, 5) scenario 2
                        zinb = sim[[3]], #distribution c(NB,ZINB)
                        mod.lib = 1, mod.shape = 1, 
                        up = 0.5, #symmetry
                        replace = TRUE, 
                        # ingroup=which(sample.of.origin=="grp2"),
                        simName = simName)
  return(obj.sim)
})

for(i in 1:length(sims)){ sims[[i]]$number <- i }
names(sims) <- lapply(sims,function(sim) sim$name)

## Depending on which dataset you chose at the beginning
### ALZ
save(sims, file = "/blackhole/alessia/CircModel/parametric_sim/ALZData_detmet_parametricsimulations.RData")

## simulate data for glmm model
# Generating simulations
dataset <- names(data_list)[1:6] # circRNA detection methods
simulation = 1:30
sampleSize = c(3,10,20)
TPR = c(0.1)
foldEffect <- c(1.5)
#sparsityEffect <- c(0.15, 0.30)
niter = 30

simulation_flow.glmm <- data.frame(expand.grid(simulation = simulation,
                                          dataset = dataset,
                                          sampleSize = sampleSize,
                                          TPR = TPR,
                                          foldEffect = foldEffect 
                                          #sparsityEffect = sparsityEffect
))
rownames(simulation_flow.glmm) <- 1:nrow(simulation_flow.glmm)
dim(simulation_flow.glmm)

# Unique seed for each mix of variables
for(i in 1:nrow(simulation_flow.glmm)){
  simulation_flow.glmm$seed[i] <- strtoi(paste0("0x",substr(digest::digest(simulation_flow.glmm[i,1:ncol(simulation_flow.glmm)],algo = "sha256"),
                                                       start = 1, stop = 7)))
}

simulation_flow.glmm = simulation_flow.glmm[order(simulation_flow.glmm$simulation, simulation_flow.glmm$dataset, simulation_flow.glmm$sampleSize), ]
simulation_flow.glmm$sims.glmm <- rep(seq(1,3,1), 180)

simulation_flow.glmm <- simulation_flow.glmm[order(simulation_flow.glmm$sampleSize),]

### ALZ
save(simulation_flow.glmm,file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_simulation_flow.RData")

sims.glmm <- apply(simulation_flow.glmm, 1, function(sim){
  
  sim = simulation_flow.glmm[4,]
  
  simName <- paste(colnames(simulation_flow.glmm),sim[1:ncol(simulation_flow.glmm)],
                   sep = ":",
                   collapse = "_")
  cat(simName,"\n")
  set.seed(sim[[ncol(simulation_flow.glmm)]])
  # tot_zinb_mu <- sum(exp(true_model@beta_mu))
  # link function is log
  # zinb_mu_rel <- as.numeric(exp(true_model@beta_mu))/tot_zinb_mu 
  # zinb_mu <- as.numeric(true_model@beta_mu)
  ncircular = 5000
  true_model <- models[[sim[[2]]]]
  range <- min(rbindlist(lapply(models, function(x) data.frame(length(x$sample.mu)))))
  chosen <- sample(range, ncircular, replace = TRUE)
  COUNT_CIRC <- COUNT_FUN_CIRC(sample.mu = true_model$sample.mu, 
                               sample.disp = true_model$sample.disp, 
                               zi.mu = true_model$zi.mu, 
                               zi.disp = true_model$zi.disp, 
                               zi.prop = true_model$zi.prop, 
                               relative.size = true_model$relative.size)
  obj.sim <- COUNT_CIRC(ncircular = ncircular, #n.circRNAs fixed 10.000
                        #pzero = sim[[4]], #SparcityEffect c(0,0.15) scenario 4 
                        nde = NULL, #n.DECs 
                        pi0 = 1-as.numeric(sim[[4]]), #1-TPR c(0.9, 0.5) scenario 3
                        m = as.numeric(sim[[3]]), #n.samples.per.group c(3, 5)
                        fc = as.numeric(sim[[5]]), #FC c(2, 5) scenario 2
                        zinb = "ZINB", #distribution c(NB,ZINB)
                        mod.lib = 1, mod.shape = 1, 
                        up = 0.5, #symmetry
                        replace = TRUE, 
                        # ingroup=which(sample.of.origin=="grp2"),
                        simName = simName)
  return(obj.sim)
})

sims.glmm <- list()
for (k in seq(1,nrow(simulation_flow.glmm)-1, by = 6)) { #for each scenario
  # k = 181
  range <- min(rbindlist(lapply(models, function(x) data.frame(length(x$sample.mu)))))
  ncircular = 5000
  chosen <- sample(range, ncircular, replace = TRUE)
  
  nde = (1-(1-simulation_flow.glmm[k,4]))*ncircular
  ## expected false positives
  FP <- round(ncircular * (1-nde/ncircular))
  TP <- ncircular - FP 
  
  ## types of true positives
  TP_up <- round(TP * 0.5)
  TP_down <- TP - TP_up 
  
  de <- c(rep(0, FP), rep(1, TP_up), rep(-1, TP_down))
  de <- de[sample.int(length(de))] ## resample
  chosen.de <- which(de!=0)
  is.de <- logical(ncircular)
  is.de[chosen.de] <- TRUE
  
  # log fold change, approximately half positive, half negative
  delta <- rep(0, ncircular)
  lfc <- sim[[5]]
  
  delta[de != 0] <- lfc * de[de != 0]
  
  sims.glmm[[k]] <- apply(simulation_flow.glmm[k:(k+6-1),], 1, function(sim){
  
  # sim = simulation_flow.glmm[k:(k+6-1),][1,]
  
  simName <- paste(colnames(simulation_flow.glmm),sim[1:ncol(simulation_flow.glmm)],
                   sep = ":",
                   collapse = "_")
  cat(simName,"\n")
  set.seed(sim[[ncol(simulation_flow.glmm)]])
  # tot_zinb_mu <- sum(exp(true_model@beta_mu))
  # link function is log
  # zinb_mu_rel <- as.numeric(exp(true_model@beta_mu))/tot_zinb_mu 
  # zinb_mu <- as.numeric(true_model@beta_mu)
  true_model <- models[[sim[[2]]]]
  sample.mu = true_model$sample.mu 
  sample.disp = true_model$sample.disp
  zi.mu = true_model$zi.mu
  zi.disp = true_model$zi.disp
  zi.prop = true_model$zi.prop
  relative.size = true_model$relative.size
  relative.size <- relative.size/mean(relative.size) # mean centering
  m = as.numeric(sim[[3]]) #n.samples.per.group c(3, 5)
  fc = as.numeric(sim[[5]]) #FC c(2, 5) scenario 2
  zinb = "ZINB" #distribution c("NB","ZINB")
  mod.lib = 1
  mod.shape = 1
  #counts <- matrix(0, nrow = ncircular, ncol = 2 * m)
  
  true_means <- true_disps <- rep(0, ncircular)
  #chosen <- sample(length(sample.mu), ncircular, replace = replace)
  
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
  
  lambda <- true_means * matrix(rep(1, ncircular), ncol = as.numeric(m) * 2, nrow = ncircular)
  
  lambda <- matrix(true_means, ncol = 1) %*%
    matrix(rep(1, 2 * as.numeric(m)), nrow = 1) * #intercept
    cbind(matrix(rep(1, ncircular * as.numeric(m)), ncol = as.numeric(m)),
          matrix(rep(exp(delta), as.numeric(m)), ncol = as.numeric(m)))
  
  ## mean of counts
  phi <- matrix(rep(true_disps, 2 * as.numeric(m)), ncol = 2 * m)
  
  ## dispersion of counts
  delta <- delta / log(2)
  
  nsamples <- 2 * m
  
  # Simulating counts.
  counts <- matrix(rnbinom(ncircular * 2 * as.numeric(m), mu = lambda, 
                           size = 1/true_disps), ncol = 2 * as.numeric(m), nrow = ncircular)
  
  # Information about True Positives
  CIRC_names <- paste("circ_",1:nrow(counts),sep = "") 
  newSampleNamesUp <- paste0(CIRC_names[de == 1], "-TPup")
  newSampleNamesDown <- paste0(CIRC_names[de == -1], "-TPdown")
  CIRC_names[de == 1] <- newSampleNamesUp
  CIRC_names[de == -1] <- newSampleNamesDown

  # Sample names
  sample_names <- c(paste("Sample_",1:as.numeric(m),"_grp1",sep=""),paste("Sample_",(as.numeric(m)+1):(2*as.numeric(m)),"_grp2",sep=""))
  colnames(counts) <- sample_names
  rownames(counts) <- CIRC_names
  # Trim too rare circRNAs
  # counts_filtered <- simpleTrimCirc(counts)
  # counts_filtered_NA <- counts_filtered[!is.na(rownames(counts_filtered)),]
  
  mu_rel_fc <- data.frame(apply(lambda, 1, mean))
  rownames(mu_rel_fc) <- CIRC_names
  
  obj <- list(ps = phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                            sample_data(data.frame(grp = rep(c("grp1","grp2"),
                                                             each = as.numeric(m)),
                                                   row.names = sample_names))),
              name = simName,
              beta = log(mu_rel_fc[rownames(counts),]))
  return(obj)
  })

}

length(sims.glmm)
sims.glmm2 = sims.glmm[-which(sapply(sims.glmm, is.null))]
length(sims.glmm2)

names(sims.glmm2) <- apply(simulation_flow, 1, function(sim) paste(colnames(simulation_flow)[-c(2,3)],
                                                                        sim[-c(2,3)],
                                                                        sep = ":",
                                                                        collapse = "_"))

## Depending on which dataset you chose at the beginning
### ALZ
save(sims.glmm2, file = "/blackhole/alessia/CircModel/parametric_sim/ALZData_glmm_parametricsimulations.RData")
length(sims.glmm2)
