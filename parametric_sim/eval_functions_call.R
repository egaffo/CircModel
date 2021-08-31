### This code was supplied to sbatch, e.g.:
#                                       to be generated     output filedir     log                  errors
#                                             |                       |         |                     |
# Rscript ./power/eval_functions_call.R ./data/sims.RData ./data/evals.RDS > ./data/sims.log 2> ./data/sims.err

source("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/reference/R/eval_functions.R")
source("/blackhole/alessia/circzi/SigEMD-master/FunImpute.R")
source("/blackhole/alessia/circzi/SigEMD-master/SigEMDHur.R")
source("/blackhole/alessia/circzi/SigEMD-master/SigEMDnonHur.R")
source("/blackhole/alessia/circzi/SigEMD-master/plot_sig.R")


list.function.files <- paste("/blackhole/alessia/circzi/glmm/NBZIMM/R/", 
                             c(list.files(path = "/blackhole/alessia/circzi/glmm/NBZIMM/R/", 
                                          pattern = "*.R"), list.files(path = "/blackhole/alessia/circzi/glmm/NBZIMM/R/", 
                                                                       pattern = "*.r")), sep = "")

lapply(list.function.files, source)

list.function.files <- paste("/blackhole/alessia/circzi/glmm_DE_circular/R/",
                             c(list.files(path = "/blackhole/alessia/circzi/glmm_DE_circular/R/",
                                        pattern = ".*R"), list.files(path = "/blackhole/alessia/circzi/glmm_DE_circular/R/", 
                                                                     pattern = "*.r")), sep = "")

lapply(list.function.files, source)

register(SerialParam())
args = commandArgs(trailingOnly=TRUE)
args[1] <- "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/ZhengData_detmet_simulationsZINB.RData"
args[2] <- "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/power/Zheng_detmet_evals_simulations_ZINB_power.RDS"
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("At least two argument must be supplied (input_file and output_file)", call.=FALSE)
}  
### simulations are required in input
input_file <- args[1]
### name for the output is the second input to supply
output_file <- args[2]

load(file = input_file)
load(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng_detmet_simulationflow.RData")
simulation_flow
start_time <- Sys.time()
evals <- lapply(X = sims[1:240], FUN = function(sim){
  cat(green("Simulation",sim$number,"\n"))
  cat(green("Simulation name:",sim$name,"\n"))
  cat("Overall time:",Sys.time()-start_time,"\n")
  start <- Sys.time()
  eval <- oneSimRunGSOwn(physeq = sim$ps, beta = sim$beta)
  end <- Sys.time()
  cat(magenta("dataset evaluation time:",end-start,"\n"))
  return(eval)
})
saveRDS(evals,file = output_file)

register(SerialParam())
args = commandArgs(trailingOnly=TRUE)
args[1] <- "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/ZhengData_glmm_simulationsZINB.RData"
args[2] <- "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/power/Zheng_glmm_evals_simulations_ZINB_power.RDS"
if (length(args)!=2) {
  stop("At least two argument must be supplied (input_file and output_file)", call.=FALSE)
}  
### simulations are required in input
input_file <- args[1]
### name for the output is the second input to supply
output_file <- args[2]
load(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng_glmm_simulation_flow.RData")
simulation_flow.glmm
library(stringr)
load(file = input_file)
start_time <- Sys.time()
load("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng_glmm_simulation_flow.RData")
det.tools <- levels(simulation_flow$dataset)

library(foreach)
evals.glmm <- list()
minutes <- c()
for (i in 1:30){#seq(1,length(sims)-1, by = 8)){
  # i = 2
  cat(green("Simulation name:",names(sims.glmm2)[[i]],"\n"))
  sim.glmm <- list()
  
  for (k in seq(1, 8, 1)){
    # k = 1
    sim.glmm[[k]] <- as.data.frame(sims.glmm2[[i]][[k]][["ps"]]@otu_table, 
                                   row.names = rownames(sims.glmm2[[i]][[k]][["ps"]]@otu_table))
  }
  names(sim.glmm) <- det.tools
  sim.glmm.db <- do.call("rbind", sim.glmm)
  sim.glmm.db$method <- sub("\\..*", "", rownames(sim.glmm.db))
  sim.glmm.db$circ_id <- sub(".*\\.", "", rownames(sim.glmm.db))
  sim.glmm.db.melt <- invisible(reshape2::melt(sim.glmm.db, id.vars = c("method", "circ_id")))
  sim.glmm.db.melt$group <- str_extract(sim.glmm.db.melt$variable, '[^_]+$')
  sim.glmm.db.melt$sample.name.ext <- paste(sim.glmm.db.melt$variable, sim.glmm.db.melt$method, sep=".")
  meta.data <- data.frame(sample = c(paste("Sample_",1:as.numeric(simulation_flow.glmm[1,3]),sep=""),
                                      paste("Sample_",(simulation_flow.glmm[1,3]+1):(2*simulation_flow.glmm[1,3]),sep="")),
                          condition = c(rep("grp1", simulation_flow.glmm[1,3]), rep("grp2", simulation_flow.glmm[1,3])))
  colnames(sim.glmm.db.melt)[colnames(sim.glmm.db.melt) == 'group'] <- 'condition'
  colnames(sim.glmm.db.melt)[colnames(sim.glmm.db.melt) == 'variable'] <- 'sampleID'
  sim.glmm.db.melt <- as.data.table(sim.glmm.db.melt)
  colData.dt <- sim.glmm.db.melt[, .N, by = .(sampleID, method, 
                                        condition, 
                                        sample.name.ext)][, N := NULL][]
  
  colData <- data.frame(colData.dt, 
                        row.names = "sample.name.ext")
  
  colData$condition <- factor(colData$condition)
  count.matrix.dt <- dcast(data = sim.glmm.db.melt, 
                           formula = circ_id ~ sample.name.ext, 
                           fill = 0, fun.aggregate = sum, 
                           value.var = "value")
  count.matrix.filtered.dt <- count.matrix.dt[rowSums(count.matrix.dt>5)>5,]
  count.matrix <- as.matrix(count.matrix.filtered.dt[,-1])[, rownames(colData)]
  rownames(count.matrix) <- count.matrix.filtered.dt$circ_id
  # head(count.matrix)
  se <- SummarizedExperiment(assays = list(counts = count.matrix),
                             colData = colData)
  count <- assay(se)
  sample <- colData
  
  glmmdb <- glmm.dataset(count = assay(se), sample = colData)
  
  #test
  start.time <- Sys.time()
  LR.test <- glmm.test(glmm.data = glmmdb, condition = "condition", test = "glmm", 
                       use.zeroinfl = TRUE, coldata = colData, family = "nbinom2")

  stop.time <- Sys.time()
  minutes[i] <- round(difftime(stop.time, start.time, units = "min"), 3)
  cat("Computational time:", minutes[i], "minutes \n")

  pValMat <- data.frame(pvalue = rbindlist(lapply(LR.test$results$PValue.wald, function(x) as.data.frame(x)))$x)
  pValMat$padj = p.adjust(pValMat$pvalue,
                             method = "BH")
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- data.frame(LRT = rbindlist(lapply(LR.test$results$statistic, function(x) as.data.frame(x)))$x)
  
  evals.glmm[[i]] <- list("pValMat" = pValMat, "statInfo" = statInfo)
  
}

saveRDS(evals.glmm, file = output_file)
evals.glmm

load(output_file)
