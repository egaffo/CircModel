### This code was supplied to sbatch, e.g.:
#                                       to be generated     output filedir     log                  errors
#                                             |                       |         |                     |
# Rscript ./power/eval_functions_call.R ./data/sims.RData ./data/evals.RDS > ./data/sims.log 2> ./data/sims.err

source("/blackhole/alessia/CircModel/parametric_sim/eval_functions.R")

register(SerialParam())
args = commandArgs(trailingOnly=TRUE)
args[1] <- "/blackhole/alessia/CircModel/parametric_sim/ALZData_detmet_parametricsimulations.RData"
args[2] <- "/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evalsDESeq2ZI_parametricsimulations_power.RDS"
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("At least two argument must be supplied (input_file and output_file)", call.=FALSE)
}  
### simulations are required in input
input_file <- args[1]
### name for the output is the second input to supply
output_file <- args[2]

load(file = input_file)
load(file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_simulation_flow.RData")
simulation_flow
start_time <- Sys.time()
evals <- lapply(X = sims, FUN = function(sim){
  # sim = sims$`simulation: 1_dataset:ccp2_distribution:ZINB_sampleSize: 3_TPR:0.1_foldEffect:1.5_seed:210192957`
  cat(green("Simulation",sim$number,"\n"))
  cat(green("Simulation name:",sim$name,"\n"))
  cat("Overall time:",Sys.time()-start_time,"\n")
  start <- Sys.time()
  eval <- oneSimRunDeseq2ZI(physeq = sim$ps, beta = sim$beta)
  end <- Sys.time()
  cat(magenta("dataset evaluation time:", end-start, "\n"))
  return(eval)
})

evals_all = readRDS(file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_parametricsimulations_power.RDS")

evals_all_ZI = list()
for(i in 1:length(evals_all)){
  evals_all_ZI[[i]] = append(evals_all[[i]], list(evals[[i]]$DESeq2_ZI))  
  names(evals_all_ZI)[[i]] = names(evals_all)[[i]]
}
#saveRDS(evals_all_ZI,file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_all_parametricsimulations_power.RDS")

########################################################
####### GLMM

register(SerialParam())
args = commandArgs(trailingOnly=TRUE)
args[1] <- "/blackhole/alessia/CircModel/parametric_sim/ALZData_glmm_parametricsimulations.RData"
args[2] <- "/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_evals_parametricsimulations_power.RDS"
if (length(args)!=2) {
  stop("At least two argument must be supplied (input_file and output_file)", call.=FALSE)
}  
### simulations are required in input
input_file <- args[1]
### name for the output is the second input to supply
output_file <- args[2]
load(file = input_file)
load(file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_simulation_flow.RData")
library(stringr)

start_time <- Sys.time()
det.tools <- levels(simulation_flow$dataset)

library(foreach)
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

evals.glmm <- list()
minutes <- c()
eval.glmm = foreach::foreach (i=c(1:30, 61:90)) %dopar% {#seq(1,length(sims)-1, by = 8)){
  # i = 1
  cat("Simulation name:",names(sims.glmm2)[[i]],"\n")
  sim.glmm <- list()
  
  for (k in seq(1, 6, 1)){
    # k = 1
    sim.glmm[[k]] <- as.data.frame(sims.glmm2[[i]][[k]][["ps"]]@otu_table, 
                                   row.names = rownames(sims.glmm2[[i]][[k]][["ps"]]@otu_table))
  }
  
  names(sim.glmm) <- det.tools
  sim.glmm.db <- do.call("rbind", sim.glmm)
  sim.glmm.db$method <- sub("\\..*", "", rownames(sim.glmm.db))
  sim.glmm.db$circ_id <- sub(".*\\.", "", rownames(sim.glmm.db))
  sim.glmm.db.melt <- invisible(reshape2::melt(sim.glmm.db, id.vars = c("method", "circ_id")))
  sim.glmm.db.melt$group <- stringr::str_extract(sim.glmm.db.melt$variable, '[^_]+$')
  sim.glmm.db.melt$sample.name.ext <- paste(sim.glmm.db.melt$variable, sim.glmm.db.melt$method, sep=".")
  meta.data <- data.frame(sample = c(paste("Sample_",1:as.numeric(simulation_flow[i,4]),sep=""),
                                     paste("Sample_",(simulation_flow[i,4]+1):(2*simulation_flow[i,4]),sep="")),
                          condition = c(rep("grp1", simulation_flow[i,4]), rep("grp2", simulation_flow[i,4])))
  colnames(sim.glmm.db.melt)[colnames(sim.glmm.db.melt) == 'group'] <- 'condition'
  colnames(sim.glmm.db.melt)[colnames(sim.glmm.db.melt) == 'variable'] <- 'sampleID'
  sim.glmm.db.melt <- data.table::as.data.table(sim.glmm.db.melt)
  colData.dt <- sim.glmm.db.melt[, .N, by = .(sampleID, method, 
                                              condition, 
                                              sample.name.ext)][, N := NULL][]
  
  colData <- data.frame(colData.dt, 
                        row.names = "sample.name.ext")
  library(dplyr)
  library(reshape2)
  library(glmmTMB)
  pheno <- colData %>%
    dplyr::rename(SampleID = sampleID) %>% 
    dplyr::rename(MethodID = method) %>% 
    dplyr::rename(condition = condition)
  count.glmm = reshape2::dcast(sim.glmm.db.melt, circ_id~sample.name.ext, fill=0)
  # count.filt.glmm = count.glmm[count.glmm$circ_id%in%rownames(sims[[i]]$ps@otu_table),]
  # count.filt.glmm = CREART::smallest_group_filter(x = as.data.table(count.glmm), 
  #                               cond = as.data.table(colData),
  #                               rthr = 1)
  count.matrix.glmm = as.matrix(count.glmm[,-1])
  rownames(count.matrix.glmm) = count.glmm$circ_id
  #count.matrix.glmm["circ_145-TPup",]
  circularcounts <- count.matrix.glmm[, rownames(pheno)]
  colnames(circularcounts) = sub("\\..*", "", colnames(circularcounts))
  # dim(circularcounts)
  # testIDx = rownames(circularcounts)
  dge <- edgeR::DGEList(circularcounts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  offsets <- dge$samples$norm.factors # norm.factors per samples
  pheno$condition <- as.factor(pheno$condition)
  pheno$MethodID <- as.factor(pheno$MethodID)
  pheno$ln.lib.size <- offsets
  # circularcounts <- log(sweep(circularcounts,2,apply(circularcounts,2,mean),'/'))
  circularcounts <- t(circularcounts)
  allsamples_test <- cbind(pheno,circularcounts) %>% as.data.frame()
  
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
  lfc <- as.numeric(unlist(lapply(summaries, function(x){x$coefficients$cond["conditiongrp2",1]})))
  pValMat = data.frame(pvalue = pvalues,
                              row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
  pValMat$padj = p.adjust(pValMat$pvalue,
                          method = "BH")
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- data.frame(LFC = lfc, row.names = colnames(allsamples_test)[5:ncol(allsamples_test)])
  
  evals.glmm[[i]] <- list("pValMat" = pValMat, "statInfo" = statInfo)
  return(evals.glmm[[i]])
}

names(eval.glmm) = names(sims.glmm2)[c(1:30, 61:90)]
# saveRDS(eval.glmm, file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_parametricsimulations_power_S1.RDS")
# saveRDS(eval.glmm, file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_parametricsimulations_power_S2.RDS")
saveRDS(eval.glmm, file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_parametricsimulations_power_S1S3.RDS")

eval.S1S3 = readRDS("/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_parametricsimulations_power_S1S3.RDS")

evals.glmm = list(eval.S1S3[c(1:30)],
                  readRDS("/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_parametricsimulations_power_S2.RDS"),
                  eval.S1S3[c(31:60)])


evals_all = readRDS("/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_all_parametricsimulations_power.RDS")
evals_all_glmm = list()
for(k in 1:3){
    if(k==1){
      # k=1
      for(i in c(1:30)){
        evals_all_glmm[[i]] = append(evals_all[[i]], list(evals.glmm[[k]][[i]]))  
        names(evals_all_glmm[[i]]) = c(names(evals_all[[i]])[-15], "DESeq2-ZI", "GLMM")
      }      
    }
    if(k==2){
      #k=2
      for(i in 31:60){
        evals_all_glmm[[i]] = append(evals_all[[i]], list(evals.glmm[[k]][[i-30]]))  
        names(evals_all_glmm[[i]]) = c(names(evals_all[[i]])[-15], "DESeq2-ZI", "GLMM")
      }
    }
    if(k==3){
        for(i in 61:90){
          evals_all_glmm[[i]] = append(evals_all[[i]], list(evals.glmm[[k]][[i-60]]))
          names(evals_all_glmm[[i]]) = c(names(evals_all[[i]])[-15], "DESeq2-ZI", "GLMM")
        }
    }
}

length(evals_all_glmm)
names(evals_all_glmm) = names(sims.glmm2)
saveRDS(evals_all_glmm, file = "/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_allGLMM_parametricsimulations_S123_power.RDS")
