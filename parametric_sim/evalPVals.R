library(plyr)
library(AUC)

methods = c("DESeq2_scRNAseq","DESeq2",
            "DESeq2_poscounts_zinbwave", 
            "DESeq2_poscounts_apeglm")
evalPVals <- function(resi, alpha = 0.05, pvalsType = "adjP", rawPvalsType = "rawP") {
  # Rarely a fit has failed, then we return 0 for sens, 1 for specificity and
  # NA for FDR, AUC and the lib/cons areas
  if (!is.matrix(resi)) {
    cat("Resi is not a matrix! \n")
    return(c(Sensitivity = 0, Specificity = 1, FDR = 0, AUC = 0.5))
  }
  # Some DESeq2 results (for example) had NA adjusted p-values Replace NA
  # $adjP values to highest possible value (1.0)
  # resi[is.na(resi[, pvalsType]), pvalsType] <- 1
  # Or just count them
  NA_prop <- sum(is.na(resi[, pvalsType]))/nrow(resi)
  
  # Evaluate detection performance.
  wh.pred = (resi[, pvalsType] < alpha)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("TP", rownames(resi))
  # Calc number of differentially abundant circRNAs
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  # Sensitivity: True Positives divided by all positives (sum of true
  # positives and false negatives)
  Sensitivity = TPs/(TPs + FNs)
  # Specificity: True negatives divided by all negatives (sum of true
  # negatives and false positives)
  Specificity = TNs/(TNs + FPs)
  # false discovery rate: false positives divided by all detected positives
  FDR = if ((TPs + FPs) == 0) 
    0 else FPs/(TPs + FPs)
  
  # If no true positives, return NA's for irrelevant measures
  wh.truth = (1:nrow(resi) %in% wh.TP)
  
  # IF AUC cannot be calculated, return NA
  rocObj = try(AUC::roc(1 - resi[, pvalsType], factor(as.numeric(wh.truth))))
  return(c("NA_Proportion" = NA_prop, 
           TPs = length(wh.TP),
           DECs = length(wh.TP),
           Sensitivity = Sensitivity, 
           Specificity = Specificity, 
           FDR = FDR, 
           AUC = ifelse(class(rocObj)[1] == "try-error", NA, AUC::auc(rocObj))))
}  # END - function: evalPVals

df_creator <- function(evals_file, sim_flow_file, out_dir){
  cat("Reading evals","\n")
  evals <- readRDS(file = evals_file)
  evals_noTM <- lapply(evals, function(methods){methods %>% purrr::list_modify("truemodel" = NULL)})
  evals_noTM <- lapply(evals_noTM, function(methods){methods %>% purrr::list_modify("Y" = NULL)})
  load(file = sim_flow_file)
  cat("Creating data.frame from evals","\n")
  eval_stats <- ldply(.data = evals_noTM,.fun = function(methods){
    ldply(.data = methods,.fun = function(m){
      evalPVals(resi = m$pValMat, alpha = 0.05, pvalsType = "adjP", rawPvalsType = "rawP")
    })
  })
  colnames(eval_stats) <- c("method",colnames(eval_stats)[-1])
  nmethods <- length(unique(eval_stats$method))
  simulation_flow_df <- apply(simulation_flow[1:length(evals_noTM),], 2, 
                              function(col) sapply(col,function(cell) rep(cell,each = nmethods)))
  evals_stats_df <- data.frame(eval_stats,simulation_flow_df) 
  
  cat("Computing ROC from pVals","\n")
  eval_ROC <- ldply(.data = evals_noTM,.fun = function(methods){
    ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-m$pValMat[,"adjP"], labels = as.factor(grepl(pattern = "TP", x = rownames(m$pValMat))))
      AUC = pROC::roc(as.factor(grepl(pattern = "TP", x = rownames(m$pValMat))), 1-m$pValMat[,"adjP"])$auc
      cbind(fpr = ROC$fpr, tpr = ROC$tpr, auc = AUC)
    })
  })
  colnames(eval_ROC) <- c("method",colnames(eval_ROC)[-1])
  lengths_ROC <- ldply(.data = evals_noTM,.fun = function(methods){
    sum(ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-m$pValMat[,"adjP"], labels = as.factor(grepl(pattern = "TP",x = rownames(m$pValMat))))
      return(length(ROC$tpr))
    })$V1)
  })
  simulation_flow_ROC_df <- apply(simulation_flow[1:length(evals_noTM),], 2, 
                                  function(col) unlist(mapply(col,lengths_ROC$V1,FUN = function(cell,times) rep(x = cell,times)),
                                                       use.names = FALSE))
  evals_ROC_df <- cbind(eval_ROC,simulation_flow_ROC_df) 
  cat("Summarizing ROC values","\n")
  evals_ROC_summary_df <- ddply(.data = evals_ROC_df[,-ncol(evals_ROC_df)],.variables = ~ 
                                  method + 
                                  dataset + 
                                  distribution + 
                                  sampleSize +
                                  simulation +
                                  TPR +
                                  foldEffect, 
                                  # compensation + 
                                  # sparsityEffect,
                                .fun = function(x){
                                    support <- seq(0,1,length.out = 101)
                                    fpr_tpr <- data.frame(fpr = 0, tpr = 0)
                                    for(i in 2:length(support)){
                                      fpr_s <- support[i-1]
                                      if(sum(x$fpr>=support[i-1] & x$fpr<support[i]) > 0)
                                        tpr_s <- mean(x$tpr[x$fpr>=support[i-1] & x$fpr<support[i]])
                                      else tpr_s <- fpr_tpr[i-2,2]
                                      fpr_tpr[i-1,] <- c(fpr_s,tpr_s)
                                    }
                                    fpr_tpr[1,] <- c(0,0)
                                    fpr_tpr[length(support),] <- c(1,1)
                                    return(fpr_tpr)
                                  })
  evals_ROC_summary_mean_df <- ddply(.data = evals_ROC_summary_df,.variables = ~ 
                                       method + 
                                       dataset + 
                                       distribution + 
                                       sampleSize +
                                       simulation +
                                       TPR +
                                       foldEffect + 
                                       fpr,
                                     # compensation + 
                                     # sparsityEffect,
                                     .fun = function(x){
                                         tpr = mean(x$tpr)
                                         se = sqrt(var(x$tpr))
                                         return(data.frame(tpr = tpr, se = se))
                                       })
  cat("Saving data","\n")
  saveRDS(evals_stats_df,file = paste0(out_dir,"evals_statsFC5zinb_df.RDS"))
  saveRDS(evals_ROC_df,file = paste0(out_dir,"evals_ROCFC5zinb_df.RDS"))
  saveRDS(evals_ROC_summary_df,file = paste0(out_dir,"evals_ROCFC5zinb_summary_df.RDS"))
  saveRDS(evals_ROC_summary_mean_df,file = paste0(out_dir,"evals_ROCFC5zinb_summary_mean_df.RDS"))
}

### Example code to generate power data.frames 
### The simulation files are heavy, for this reason they are not saved in github
### However the final data.frame is available.
df_creator(evals_file="/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/power/Zheng_detmet_evalsDESEq2_simulations_ZINB_FC5_power.RDS", #from eval_function_call.R
           sim_flow_file="/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng_detmet_simulationFC5_flow.RData", #form simulator
           out_dir="/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng/")

evals_ROC_df <- readRDS("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng/evals_ROCFC5zinb_df.RDS")
evals_ROC_df$method <- factor(evals_ROC_df$method)
evals_ROC_df$method <- factor(evals_ROC_df$method, levels = levels(evals_ROC_df$method), 
                              labels = c("DESeq2", "ape-glm", "DESeq2 + ZINB-WaVe", "DESeq2 (scRNAseq opts.)"))
detection.labels <- c("ccp", "circFinder", "circExplorer2 - bwa", 
                      "circExplorer2 - STAR", "circExplorer2 - tophat",
                      "CIRI", "DCC",  "findcirc")
names(detection.labels) <- levels(evals_ROC_df$dataset)

library(RColorBrewer)
# display.brewer.all()
cols <- c(
  # DEseq
  brewer.pal(n = 9, "YlOrRd")[c(4,5,6,7,9)],
  # Edger
  brewer.pal(n = 9, "GnBu")[c(6,7,8,9)],
  # limma
  brewer.pal(n = 9, "RdPu")[c(5,7)]
)
methods <- c("DESeq2", "ape-glm", "DESeq2 + ZINB-WaVe", "DESeq2 (scRNAseq opts.)")
names(cols) <- methods
method_cols <- c(
  brewer.pal(n=8, "Dark2")
)

tiff(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/power/AUC_ZhengSimFC5zinb.tiff",
     width = 10, height = 10, units = "in", res = 300)
ggplot(evals_ROC_df, aes(x=method, y = auc, color = method)) + 
  geom_boxplot() + 
  ylab("AUC") +
  xlab("") +
  labs(fill = "DE methods") +
  scale_color_manual(values = cols) +
  facet_wrap( ~ dataset,  
             labeller = labeller(dataset = detection.labels)) +
  theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                   size = 7, angle = 45, vjust=.9, hjust=0.8),
        strip.text.x = element_text(face = "bold.italic", size = 7),
        strip.text.y = element_text(face = "bold.italic", size = 7), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1)) +
  ggtitle("FC=5")
dev.off()

evals_stat_df <- readRDS("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng/evals_statsFC5zinb_df.RDS")
evals_stat_df <- readRDS("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/data/Zheng/evals_statsFC2zinb_df.RDS")

evals_stat_df$method <- factor(evals_stat_df$method)
evals_stat_df$method <- factor(evals_stat_df$method, levels = levels(evals_stat_df$method),
                               labels = c("DESeq2", "ape-glm", "DESeq2 + ZINB-WaVe", "DESeq2 (scRNAseq opts.)"))
evals_stat_df$dataset <- factor(evals_stat_df$dataset)
evals_stat_df$dataset <- factor(evals_stat_df$dataset, levels = levels(evals_stat_df$dataset),
                               labels = c("ccp", "circFinder", "circExplorer2 - bwa", 
                                          "circExplorer2 - STAR", "circExplorer2 - tophat",
                                          "CIRI", "DCC",  "findcirc"))
p <- evals_stat_df %>% 
  dplyr::group_by(method, dataset) %>% 
  dplyr::summarise(sens.mean = mean(Sensitivity),
            spec.mean  = mean(Specificity)) %>% 
  ggplot(aes(y=sens.mean, x=1-spec.mean, color=method, shape=method))
tiff(file = "/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/power/sens_specZhengSimFC2zinb.tiff",
     width = 10, height = 10, units = "in", res = 300)
p + geom_point(size = 3) + theme_bw() + facet_wrap( ~ dataset) +
  scale_shape_manual(values=1:5) +
  scale_color_manual(values = cols) +
  ylab("Sensitivity") +
  xlab("1 - specificity (false positive rate)") + 
  # coord_cartesian(xlim=c(-.003,.035)) + 
  geom_vline(xintercept=.05) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(#face = "bold", color = "#993333",
                                   size = 7, angle = 45, vjust=.9, hjust=0.8),
        strip.text.x = element_text(face = "bold.italic", size = 7),
        strip.text.y = element_text(face = "bold.italic", size = 7), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1)) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + 
  ggtitle("FC = 2")
  # scale_x_continuous(breaks=c(0,.1))
dev.off()


