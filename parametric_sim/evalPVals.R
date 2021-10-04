library(plyr)
library(AUC)
library(dplyr)

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
  FDR = if ((TPs + FPs) == 0) 0 else FPs/(TPs + FPs)
  
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


methods_Sel = c("circMeta","DESeq2_poscounts_gampoi","DESeq2_poscounts_zinbwave","DESeq2_poscounts",
                "limma_voom_TMM_zinbwave","limma_voom_TMM",
                "edgeR_TMM_zinbwave","edgeR_TMM_standard","edgeR_TMM_robustDisp", "DESeq2-ZI", "GLMM")                         
library(purrr)
evals_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_allGLMM_parametricsimulations_S123_power.RDS" #from eval_function_call.R
sim_flow_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_simulation_flow.RData"
df_creator <- function(evals_file, sim_flow_file, out_dir){
  cat("Reading evals","\n")
  evals <- readRDS(file = evals_file)
  evals_noTM <- lapply(evals, function(methods){methods %>% purrr::list_modify("truemodel" = NULL)})
  evals_noTM <- lapply(evals_noTM, function(methods){methods %>% purrr::list_modify("Y" = NULL)})
  evals_noTM = purrr::map(evals_noTM, function(x) keep(x, .p = names(x)%in%methods_Sel))
  load(file = sim_flow_file)
  cat("Creating data.frame from evals","\n")
  
  # evals.glmm_df = foreach(i=1:3, .combine = rbind) %dopar% {
  #   # i=2
  #   eval.glmm <- evals.glmm[[i]]
  #   eval.glmm_noTM <- lapply(eval.glmm, function(methods){methods %>% purrr::list_modify("truemodel" = NULL)})
  #   eval.glmm_noTM <- lapply(eval.glmm_noTM, function(methods){methods %>% purrr::list_modify("Y" = NULL)})
  #   eval_stats.glmm <- ldply(.data = eval.glmm_noTM,.fun = function(m){
  #     evalPVals(resi = as.matrix(m$pValMat), alpha = 0.05, pvalsType = "rawP", rawPvalsType = "rawP")
  #     })
  #   eval_stats.glmm$method = "GLMM"
  #   if(i==1){evals_stats.glmm_df = data.frame(eval_stats.glmm[,-1], simulation_flow[i:length(eval.glmm_noTM),])}
  #   if(i==2){evals_stats.glmm_df = data.frame(eval_stats.glmm[,-1], simulation_flow[31:60,])}
  #   if(i==3){evals_stats.glmm_df = data.frame(eval_stats.glmm[,-1], simulation_flow[61:90,])}
  #                                
  #   evals_stats.glmm_df
  # }

  eval_stats <- ldply(.data = evals_noTM,.fun = function(methods){
    #e=evals_noTM$`simulation: 1_sampleSize: 3_TPR:0.1_foldEffect:1.5_seed:210192957`
    ldply(.data = methods,.fun = function(m){
      #m=e$GLMM
      evalPVals(resi = as.matrix(m$pValMat), alpha = 0.05, pvalsType = "rawP", rawPvalsType = "rawP")
    })
  })
  
  colnames(eval_stats) <- c("method",colnames(eval_stats)[-1])
  nmethods <- length(unique(eval_stats$method))
  simulation_flow_df <- apply(simulation_flow[1:length(evals_noTM),], 2, 
                              function(col) sapply(col,function(cell) rep(cell,each = nmethods)))
  evals_stats_df <- data.frame(eval_stats,simulation_flow_df) 

  evals_stats_df$method <- factor(evals_stats_df$method)
  evals_stats_df$method <- factor(evals_stats_df$method, levels = levels(evals_stats_df$method), 
                                labels = c("circMeta",
                                           "DESeq2",
                                           "DESeq2-GamPoi",
                                           "DESeq2-ZINB Wave",
                                           "DESeq2-ZI",
                                           "edgeR-robust",
                                           "edgeR",
                                           "edgeR-ZINB Wave",
                                           "GLMM",
                                           "voom",
                                           "voom-ZINB Wave"))
  cat("Computing ROC from pVals","\n")
  
  eval_ROC <- ldply(.data = evals_noTM,.fun = function(methods){
    ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-as.matrix(m$pValMat)[,"rawP"], labels = as.factor(grepl(pattern = "TP", x = rownames(as.matrix(m$pValMat)))))
      AUC = pROC::roc(as.factor(grepl(pattern = "TP", x = rownames(as.matrix(m$pValMat)))), 1-as.matrix(m$pValMat)[,"rawP"])$auc
      cbind(fpr = ROC$fpr, tpr = ROC$tpr, auc = AUC)
    })
  })
  
  colnames(eval_ROC) <- c("method",colnames(eval_ROC)[-1])
  eval_ROC$method <- factor(eval_ROC$method)
  eval_ROC$method <- factor(eval_ROC$method, levels = levels(eval_ROC$method), 
                                labels = c("circMeta",
                                           "DESeq2",
                                           "DESeq2-GamPoi",
                                           "DESeq2-ZINB Wave",
                                           "DESeq2-ZI",
                                           "edgeR-robust",
                                           "edgeR",
                                           "edgeR-ZINB Wave",
                                           "GLMM",
                                           "voom",
                                           "voom-ZINB Wave"))
  lengths_ROC <- ldply(.data = evals_noTM,.fun = function(methods){
    sum(ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-as.matrix(m$pValMat)[,"rawP"], labels = as.factor(grepl(pattern = "TP",
                                                                                              x = rownames(as.matrix(m$pValMat)))))
      return(length(ROC$tpr))
    })$V1)
  })
  simulation_flow_ROC_df <- apply(simulation_flow[1:length(evals_noTM),], 2, 
                                  function(col) unlist(mapply(col,lengths_ROC$V1,FUN = function(cell,times) rep(x = cell,times)),
                                                       use.names = FALSE))
  evals_ROC_df <- cbind(eval_ROC, simulation_flow_ROC_df)
  
  cat("Summarizing ROC values","\n")
  evals_ROC_summary_df <- ddply(.data = evals_ROC_df[,-ncol(evals_ROC_df)],.variables = ~ 
                                  method + 
                                  dataset + 
                                  #distribution + 
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
                                       #distribution + 
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
  write.csv(evals_stats_df,file = paste0(out_dir,"evals_stats_ALZ_df.csv"))
  write.csv(evals_ROC_df,file = paste0(out_dir,"evals_ROC_ALZ_df.csv"))
  write.csv(evals_ROC_summary_df,file = paste0(out_dir,"evals_ROCF_ALZ_summary_df.csv"))
  write.csv(evals_ROC_summary_mean_df,file = paste0(out_dir,"evals_ROC_ALZ_summary_mean_df.csv"))
}

### Example code to generate power data.frames 
### The simulation files are heavy, for this reason they are not saved in github
### However the final data.frame is available.
df_creator(evals_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_detmet_evals_all_parametricsimulations_power.RDS", #from eval_function_call.R
           sim_flow_file="/blackhole/alessia/CircModel/parametric_sim/ALZ_glmm_simulation_flow.RData", #form simulator
           out_dir="/blackhole/alessia/CircModel/parametric_sim/")

evals_ROC_df <- read.csv("/blackhole/alessia/CircModel/parametric_sim/evals_ROC_ALZ_df.csv")
evals_ROC_df$method <- factor(evals_ROC_df$method)
evals_ROC_df$method <- factor(evals_ROC_df$method, 
                              levels = c("circMeta",
                                         "DESeq2",
                                         "DESeq2-GamPoi",
                                         "DESeq2-ZINB Wave",
                                         "DESeq2-ZI",
                                         "edgeR-robust",
                                         "edgeR",
                                         "edgeR-ZINB Wave",
                                         "GLMM",
                                         "voom",
                                         "voom-ZINB Wave"),
                              ordered = T)
detection.labels <- c("ccp2")#, "circFinder", "circExplorer2 - bwa", 
                      # "circExplorer2 - STAR", "circExplorer2 - tophat",
                      # "CIRI", "DCC",  "findcirc")
names(detection.labels) <- levels(evals_ROC_df$dataset)

library(RColorBrewer)
cols <- c(
  # circMeta
  brewer.pal(n = 9, "OrRd")[c(8)],
  # DEseq
  brewer.pal(n = 9, "YlOrRd")[c(3,4,5,6)],
  # Edger
  brewer.pal(n = 9, name = "GnBu")[c(5,6,7)],
  # GLMM
  brewer.pal(n = 9, "BuPu")[c(5)],
  # limma
  brewer.pal(n = 9, "RdPu")[c(5,7)]
)

methods2 <- c("circMeta",
              "DESeq2",
              "DESeq2-GamPoi",
              "DESeq2-ZINB Wave",
              "DESeq2-ZI",
              "edgeR-robust",
              "edgeR",
              "edgeR-ZINB Wave",
              "GLMM",
              "voom",
              "voom-ZINB Wave")
names(cols) <- methods2
# png(file = "/blackhole/alessia/CircModel/parametric_sim/Figure/AUC_ALZSim.png",
#      width = 12, height = 10, units = "in", res = 300)
p1 = ggplot(evals_ROC_df, aes(x=reorder(method, -auc), y = auc, color = method)) + 
  geom_boxplot() + 
  ylab("AUC") +
  xlab("") +
  labs(fill = "DE methods") +
  scale_color_manual(values = cols) +
  facet_wrap( ~ sampleSize) + #,
             # labeller = labeller(dataset = detection.labels)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 8.5, angle = 45, vjust=.97, hjust=1.1),
        strip.text.x = element_text(face = "bold.italic", size = 9),
        strip.text.y = element_text(face = "bold.italic", size = 9), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1),
        title = element_text(size = 12)) +
  ggtitle("AUC across sample size - Dataset ALZ\n(n.sim,30; ZINB distr.; FC,1.5; TPR,0.1)")
# dev.off()

evals_stats_df <- readRDS("/blackhole/alessia/CircModel/parametric_sim/evals_stats_ALZ_df.RDS")

evals_stats_df$method <- factor(evals_stats_df$method)
evals_stats_df$method <- factor(evals_stats_df$method, 
                               levels = c("circMeta",
                                          "DESeq2",
                                          "DESeq2-GamPoi",
                                          "DESeq2-ZINB Wave",
                                          "DESeq2-ZI",
                                          "edgeR-robust",
                                          "edgeR",
                                          "edgeR-ZINB Wave",
                                          "GLMM",
                                          "voom",
                                          "voom-ZINB Wave"),
                               ordered = T)
p <- evals_stats_df %>% 
  dplyr::group_by(method, sampleSize) %>% 
  dplyr::summarise(sens.mean = mean(Sensitivity),
            spec.mean  = mean(Specificity)) %>% 
  ggplot(aes(y=sens.mean, x=1-spec.mean, color=method))

# png(file = "/blackhole/alessia/CircModel/parametric_sim/Figure/sens_spec_ALZ.png",
     # width = 5, height = 5, units = "in", res = 300)
p2 = p + geom_point(size = 3) + 
  theme_bw() + 
  facet_wrap( ~ sampleSize) +
  # scale_shape_manual(values=1:5) +
  scale_color_manual(values = cols) +
  ylab("Sensitivity") +
  xlab("1 - specificity (false positive rate)") + 
  # coord_cartesian(xlim=c(-.003,.035)) + 
  geom_vline(xintercept=.05) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(#face = "bold", color = "#993333",
                                   size = 9, angle = 45, vjust=.9, hjust=0.8),
        strip.text.x = element_text(face = "bold.italic", size = 9),
        strip.text.y = element_text(face = "bold.italic", size = 9), #angle = 75),
        strip.background = element_rect(#fill = "lightblue",
          colour = "grey", size = 1)) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) + 
  ggtitle("Power across sample size - Dataset ALZ\n(n.sim,30; ZINB distr.; FC,1.5; TPR,0.1)")
  # scale_x_continuous(breaks=c(0,.1))
# dev.off()

legend = get_legend(p + geom_point(size = 3) + 
                      theme_bw() + 
                      facet_wrap( ~ sampleSize) +
                      # scale_shape_manual(values=1:5) +
                      scale_color_manual(values = cols) +
                      ylab("Sensitivity") +
                      xlab("1 - specificity (false positive rate)") + 
                      # coord_cartesian(xlim=c(-.003,.035)) + 
                      geom_vline(xintercept=.05) +
                      theme_classic() +
                      theme(legend.position = "bottom",
                            axis.text.x = element_text(#face = "bold", color = "#993333",
                              size = 9, angle = 45, vjust=.9, hjust=0.8),
                            strip.text.x = element_text(face = "bold.italic", size = 9),
                            strip.text.y = element_text(face = "bold.italic", size = 9), #angle = 75),
                            strip.background = element_rect(#fill = "lightblue",
                              colour = "grey", size = 1)) +
                      guides(shape=guide_legend(nrow=2,byrow=TRUE)) + 
                      ggtitle("Power across sample size - Dataset ALZ\n(n.sim,30; ZINB distr.; FC,1.5; TPR,0.1)"))

prow <- plot_grid( p1 ,
                   p2 ,
                   align = 'vh',
                   labels = c("a", "b"),
                   hjust = -1,
                   nrow = 1
)


p <- plot_grid( prow, legend, ncol = 1, rel_heights = c(1, .2))

pgrid.arrange(arrangeGrob(p1,
                         p2,
                         nrow=1),
             legend, nrow=2, heights=c(10, 1))
