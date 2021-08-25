library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(phyloseq)
library(vegan)
library(reshape2)
library(ffpe)
library(cowplot)
library(scales)
library(ggdendro)
library(metacoder)



# Data retrieval
generate_dataset <- function(A,B,seed,prop){
  # create a long data 
  count.data.file <- "/blackhole/circrna/analyses/VanVlierberghe/merge_T/circRNA_expression_per_sample.csv"
  count.data <- fread(count.data.file)
  ps.ccp2 <- as.matrix(count.data[,c(1,7:ncol(count.data)),with=FALSE],
                       rownames = "circ_id")
  # Keep one sample for each Random Subject IDentifier
  # ps.ccp2 <- ps.ccp2[rowSums(ps.ccp2>=2)>=2,]
  # colnames(ps.ccp2) <- gsub("-", "_", colnames(ps.ccp2))
  # count.data.dt <- melt(ps.ccp2, id = "circ_id", variable.name = "sample")
  colnames(count.data.dt) <- c("circ_id", "sample", "row.count")
  meta.data <- read.csv("/blackhole/alessia/circzi/checkCircRNAnormalizationdistribution/robustness_glmm/TALL/meta.csv")
  coldata <- DataFrame(subgroup = meta.data$condition,
                       condition = meta.data$groups,
                       sample = meta.data$sample_id,
                       row.names = gsub("-", "_", meta.data$sample_id))
  coldata$condition <- factor(gsub("-", "_", coldata$condition))
  coldata$subgroup <- factor(coldata$subgroup)
  coldata$sample <- as.character(coldata$sample)
  
  # data_all <<- count.data.dt %>% merge(coldata, by = "sample")
  
  ps_all = phyloseq(otu_table(ps.ccp2, taxa_are_rows = TRUE),
                    sample_data(data.frame(coldata[coldata$condition %in% c(A,B),], 
                                           row.names = coldata$sample[coldata$condition%in%c(A,B)])))
  # group B
  psB <<- subset_samples(ps_all, condition == B)
  
  # Group A
  psA <<- subset_samples(ps_all, condition == A)
  
  grp1 <<- grp2 <<- c()
  set.seed(seed)
  for(i in 1:length(psB@sam_data$sample)){
    if(psB@sam_data$sample[i] %in% psA@sam_data$sample){
      # if(sample(c(1,2),1,prob = c(0.37,0.63)) == 1){
      if(sample(c(1,2),1,prob = prop) == 1){
        grp1 <<- c(grp1,psB@sam_data$sample[i])
      } else grp2 <<- c(grp2,psB@sam_data$sample[i])
    } else grp1 <<- c(grp1,psB@sam_data$sample[i])
  }
  ps <<- merge_phyloseq(subset_samples(psB,sample %in% grp1),subset_samples(psA,sample %in% grp2))
  ps <<- filter_taxa(ps,function(x) sum(x>0)>0,1)
  return(ps)
}

#MDS
compute_MDS <- function(ps, normalization = "none" , method = "MDS", distance = "bray", color = "RUN_CENTER", names = "HMP_BODY_SUBSITE", ellipse = TRUE){
  if(normalization == "rarefy"){
    ps <- rarefy_even_depth(ps,rngseed = 123)
  } else if(normalization == "TSS"){
    ps@otu_table@.Data <- apply(ps@otu_table@.Data,2,function(col) col/sum(col)*median(colSums(ps@otu_table@.Data))) 
  } else {}
  
  ordination <- ordinate(ps, method = method, distance = distance) 
  ps@sam_data$A1 <- ordination$vectors[,1]
  ps@sam_data$A2 <- ordination$vectors[,2]
  p_o <- ggplot(data.frame(ps@sam_data), aes(x = A1, y = A2, color = eval(parse(text = color)))) +
    geom_point() +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle(label = paste0(method," on ",paste(unlist(unique(ps@sam_data[,names])),collapse = " vs ")),subtitle = paste0("Normalization: ",normalization,", Distance: ", distance)) +
    labs(color = color) + xlab(paste0("PCoA.1 (",round(ordination$values$Relative_eig[1]*100,digits = 2),"%)")) +
    ylab(paste0("PCoA.2 (",round(ordination$values$Relative_eig[2]*100,digits = 2),"%)"))
  if(ellipse) p_o <- p_o + stat_ellipse()
  
  p_l <- get_legend(ggplot(data.frame(ps@sam_data), aes(x = A1, y = A2, color = eval(parse(text = color)))) +
                      geom_point(size = 5) +
                      theme_minimal() +
                      theme(legend.position = "bottom") +
                      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
                      labs(color = color))
  
  return(list(plot = p_o,legend = p_l))
}

compute_concordance <- function(ps_fitted_list, maxrank){
  conc_df <- NULL
  for(i in 1:length(ps_fitted_list$test$lfc)){ # i in 1:n comparisons
    # i = 1
    # pval extraction
    P_df1 <- ps_fitted_list$Heldout$pval[[i]] # first half data
    P_df2 <- ps_fitted_list$test$pval[[i]] # second half data
    lfc_df1 <- ps_fitted_list$test$lfc[[i]]
    lfc_df1[is.na(lfc_df1)] = 0
    lfc_df2 <- ps_fitted_list$Heldout$lfc[[i]]
    lfc_df2[is.na(lfc_df2)] = 0
    
    nmethods <- length(names(ps_fitted_list$test$pval[[i]]))
    for(j in 1:nmethods){ # j in method names
      # j = 4
      cat("Mehod", names(ps_fitted_list$test$pval[[i]])[j],"with itself")
      vec1 = P_df1[-abs(lfc_df1[,j]),j]
      names(vec1) <- rownames(P_df1[-abs(lfc_df1[,j]),])
      vec2 = P_df2[-abs(lfc_df2[,j]),j]
      names(vec2) <- rownames(P_df2[-abs(lfc_df2[,j]),])
      # for(k in 1:nmethods){ # k in method names again
        # k = 1
        # cat("\t",names(ps_fitted_list$test[[i]])[k],"\n")
        # if(j != k){ # BMC computation
        #   # BMC for Hedlout data
        #   conc_subset1 <- data.frame(CATplot(vec1 = P_df1[,j],vec2 = P_df1[,k],make.plot = FALSE,maxrank = 100), 
        #                              method1 = names(ps_fitted_list$Heldout[[i]])[j], 
        #                              method2 = names(ps_fitted_list$Heldout[[i]])[k],
        #                              nfeatures = length(ps_fitted_list$Heldout[[i]][[j]]),
        #                              comparison = i,
        #                              subset = "Heldout")
        #   # BMC for test data
        #   conc_subset2 <- data.frame(CATplot(vec1 = P_df2[,j],vec2 = P_df2[,k],make.plot = FALSE,maxrank = 100), 
        #                              method1 = names(ps_fitted_list$test[[i]])[j], 
        #                              method2 = names(ps_fitted_list$test[[i]])[k],
        #                              #ndisc_0.1_method1 = length(adjP_df2[[j]]),
        #                              #ndisc_0.1_method2 = length(adjP_df2[[k]]),
        #                              nfeatures = length(ps_fitted_list$test[[i]][[j]]),
        #                              comparison = i,
        #                              subset = "test")
        #   conc <- rbind(conc_subset1,conc_subset2)
        # } else {
          # WMC computed between Subset1 and Subset2
          conc <- data.frame(CATplot(vec1 = vec1,
                                     vec2 = vec2, make.plot = FALSE, maxrank = maxrank), 
                             method1 = names(ps_fitted_list$test$pval[[i]])[j], 
                             method2 = names(ps_fitted_list$test$pval[[i]])[j],
                             #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                             #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                             nfeatures = mean(length(ps_fitted_list$test$pval[[i]][,j]), 
                                              length(ps_fitted_list$Heldout$pval[[i]][,j])),
                             comparison = i,
                             subset = "HeldoutvsTest")
        # }
        conc_df <- rbind(conc_df,conc)
      }
    # }
  }
  return(conc_df)
}

compute_concordance_withbetw <- function(ps_fitted_list, maxrank){
  conc_df <- NULL
  for(i in 1:length(ps_fitted_list$test$lfc)){ # i in 1:n comparisons
    # i = 3
    # pval extraction
    P_df1 <- ps_fitted_list$Heldout$pval[[i]] # first half data
    P_df2 <- ps_fitted_list$test$pval[[i]] # second half data
    lfc_df1 <- ps_fitted_list$test$lfc[[i]]
    lfc_df1[is.na(lfc_df1)] = 0
    lfc_df2 <- ps_fitted_list$Heldout$lfc[[i]]
    lfc_df2[is.na(lfc_df2)] = 0
    
    nmethods <- length(names(ps_fitted_list$test$pval[[i]]))
    for(j in 1:nmethods){ # j in method names
      # j = 1
      cat("Mehod", names(ps_fitted_list$test$pval[[i]])[j], "with") #,"with GLMM")
      vec1_pdf1 = P_df1[-abs(lfc_df1[,j]),j]
      names(vec1_pdf1) <- rownames(P_df1[-abs(lfc_df1[,j]),])
      vec1_pdf2 = P_df2[-abs(lfc_df2[,j]),j]
      names(vec1_pdf2) <- rownames(P_df2[-abs(lfc_df2[,j]),])
      
      for(k in 1:nmethods){ # k in method names again
        # k = 2
        cat("\t",names(ps_fitted_list$test$pval[[i]])[k],"\n")
        vec2_pdf1 = P_df2[-abs(lfc_df2[,k]),k]
        names(vec2_pdf1) <- rownames(P_df2[-abs(lfc_df2[,k]),])
        vec2_pdf2 = P_df2[-abs(lfc_df2[,k]),k]
        names(vec2_pdf2) <- rownames(P_df2[-abs(lfc_df2[,k]),])
        if(j != k){ # BMC computation
          # BMC for Hedlout data
          conc_subset1 <- data.frame(CATplot(vec1 = vec1_pdf1,vec2 = vec2_pdf1,
                                             make.plot = FALSE, maxrank = maxrank),
                                     method1 = names(ps_fitted_list$test$pval[[i]])[j],
                                     method2 = names(ps_fitted_list$test$pval[[i]])[k],
                                     nfeatures = length(ps_fitted_list$Heldout$pval[[i]][[j]]),
                                     comparison = i,
                                     subset = "Heldout")
          # BMC for test data
          conc_subset2 <- data.frame(CATplot(vec1 = vec1_pdf2,vec2 = vec2_pdf2,
                                             make.plot = FALSE,maxrank = maxrank),
                                     method1 = names(ps_fitted_list$test$pval[[i]])[j],
                                     method2 = names(ps_fitted_list$test$pval[[i]])[k],
                                     #ndisc_0.1_method1 = length(adjP_df2[[j]]),
                                     #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                                     nfeatures = length(ps_fitted_list$test$pval[[i]][[k]]),
                                     comparison = i,
                                     subset = "test")
          conc <- rbind(conc_subset1,conc_subset2)
        } else {
          # WMC computed between Subset1 and Subset2
          conc <- data.frame(CATplot(vec1 = vec1_pdf1,
                                     vec2 = vec1_pdf2, make.plot = FALSE, maxrank = maxrank), 
                             method1 = names(ps_fitted_list$test$pval[[i]])[j], 
                             method2 =  names(ps_fitted_list$test$pval[[i]])[j],
                             #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                             #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                             nfeatures = mean(length(ps_fitted_list$test$pval[[i]][,j]), 
                                              length(ps_fitted_list$Heldout$pval[[i]][,j])),
                             comparison = i,
                             subset = "HeldoutvsTest")
        }
        conc_df <- rbind(conc_df,conc)
      }
    }
  }
  return(conc_df)
}


compute_concordance_withGLMM <- function(ps_fitted_list, maxrank){
  conc_df <- NULL
  for(i in 1:length(ps_fitted_list$test$lfc)){ # i in 1:n comparisons
    # i = 1
    # pval extraction
    P_df1 <- ps_fitted_list$Heldout$pval[[i]] # first half data
    P_df2 <- ps_fitted_list$test$pval[[i]] # second half data
    lfc_df1 <- ps_fitted_list$test$lfc[[i]]
    lfc_df1[is.na(lfc_df1)] = 0
    lfc_df2 <- ps_fitted_list$Heldout$lfc[[i]]
    lfc_df2[is.na(lfc_df2)] = 0
    
    nmethods <- length(names(ps_fitted_list$test$pval[[i]]))
    for(j in 1:nmethods){ # j in method names
      # j = 2
      cat("Mehod", names(ps_fitted_list$test$pval[[i]])[j],"with GLMM")
      vec1 = P_df1[-abs(lfc_df1[,j]),j]
      names(vec1) <- rownames(P_df1[-abs(lfc_df1[,j]),])
      vec2 = P_df2[-abs(lfc_df2[,"GLMM_NB"]),"GLMM_NB"]
      names(vec2) <- rownames(P_df2[-abs(lfc_df2[,"GLMM_NB"]),])
      # for(k in 1:nmethods){ # k in method names again
      # k = 1
      # cat("\t",names(ps_fitted_list$test$pval[[i]][,-1])[k],"\n")
      # if(j != k){ # BMC computation
      #   # BMC for Hedlout data
      #   conc_subset1 <- data.frame(CATplot(vec1 = P_df1[,j],vec2 = P_df1[,k],make.plot = FALSE,maxrank = maxrank),
      #                              method1 = names(ps_fitted_list$Heldout[[i]])[j],
      #                              method2 = names(ps_fitted_list$Heldout[[i]])[k],
      #                              nfeatures = length(ps_fitted_list$Heldout[[i]][[j]]),
      #                              comparison = i,
      #                              subset = "Heldout")
      #   # BMC for test data
      #   conc_subset2 <- data.frame(CATplot(vec1 = P_df2[,j],vec2 = P_df2[,k],make.plot = FALSE,maxrank = maxrank),
      #                              method1 = names(ps_fitted_list$test$pval[[i]][,-1])[j],
      #                              method2 = names(ps_fitted_list$test$pval[[i]][,-1])[k],
      #                              #ndisc_0.1_method1 = length(adjP_df2[[j]]),
      #                              #ndisc_0.1_method2 = length(adjP_df2[[k]]),
      #                              nfeatures = mean(length(ps_fitted_list$test$pval[[i]][,j]), 
      #                                               length(ps_fitted_list$Heldout$pval[[i]][,k])),
      #                              comparison = i,
      #                              subset = "test")
      #   conc <- rbind(conc_subset1,conc_subset2)
      # } else {
      # WMC computed between Subset1 and Subset2
      conc <- data.frame(CATplot(vec1 = vec1,
                                 vec2 = vec2, make.plot = FALSE, maxrank = maxrank), 
                         method1 = names(ps_fitted_list$test$pval[[i]])[j], 
                         method2 = "GLMM_NB",
                         #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                         #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                         nfeatures = mean(length(ps_fitted_list$test$pval[[i]][,j]), 
                                          length(ps_fitted_list$Heldout$pval[[i]][,j])),
                         comparison = i,
                         subset = "HeldoutvsTest")
      # }
      conc_df <- rbind(conc_df,conc)
    }
    # }
  }
  return(conc_df)
}

compute_concordance.glmm <- function(ps_fitted_list, maxrank){
  conc_df <- NULL
  for(i in 1:length(ps_fitted_list$test$lfc)){ # i in 1:10 comparisons
    # i = 10
    # pval extraction
    P_df1 <- ps_fitted_list$Heldout$pval[[i]] # first half data
    P_df2 <- ps_fitted_list$test$pval[[i]] # second half data
    lfc_df1 <- ps_fitted_list$test$lfc[[i]]
    lfc_df2 <- ps_fitted_list$Heldout$lfc[[i]]
    vec1 = P_df1[rownames(-abs(lfc_df1)),1]
    names(vec1) <- rownames(-abs(lfc_df1))
    vec2 = P_df2[rownames(-abs(lfc_df2)),1]
    names(vec2) <- rownames(-abs(lfc_df2))
      # for(k in 1:nmethods){ # k in method names again
      # k = 1
      # cat("\t",names(ps_fitted_list$test[[i]])[k],"\n")
      # if(j != k){ # BMC computation
      #   # BMC for Hedlout data
      #   conc_subset1 <- data.frame(CATplot(vec1 = P_df1[,j],vec2 = P_df1[,k],make.plot = FALSE,maxrank = 100), 
      #                              method1 = names(ps_fitted_list$Heldout[[i]])[j], 
      #                              method2 = names(ps_fitted_list$Heldout[[i]])[k],
      #                              nfeatures = length(ps_fitted_list$Heldout[[i]][[j]]),
      #                              comparison = i,
      #                              subset = "Heldout")
      #   # BMC for test data
      #   conc_subset2 <- data.frame(CATplot(vec1 = P_df2[,j],vec2 = P_df2[,k],make.plot = FALSE,maxrank = 100), 
      #                              method1 = names(ps_fitted_list$test[[i]])[j], 
      #                              method2 = names(ps_fitted_list$test[[i]])[k],
      #                              #ndisc_0.1_method1 = length(adjP_df2[[j]]),
      #                              #ndisc_0.1_method2 = length(adjP_df2[[k]]),
      #                              nfeatures = length(ps_fitted_list$test[[i]][[j]]),
      #                              comparison = i,
      #                              subset = "test")
      #   conc <- rbind(conc_subset1,conc_subset2)
      # } else {
      # WMC computed between Subset1 and Subset2
      conc <- data.frame(CATplot(vec1 = vec1,
                                 vec2 = vec2, make.plot = FALSE, maxrank = maxrank), 
                         method1 = names(ps_fitted_list$test$pval[[i]])[1], 
                         method2 = names(ps_fitted_list$test$pval[[i]])[1],
                         #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                         #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                         nfeatures = mean(length(ps_fitted_list$test$pval[[i]]), 
                                          length(ps_fitted_list$Heldout$pval[[i]])),
                         comparison = i,
                         subset = "HeldoutvsTest")
      # }
      conc_df <- rbind(conc_df,conc)
    # }
  }
  return(conc_df)
}

### Compute AUC and AOC of concordance distribution 
# This method comes from p-value analysis
AucAocFun <- function(cVals, nfeatures, threshold = 0, plotIt = FALSE, ...) {
  ## sum of height differences between the curve and the y=x line
  MaxArea <- threshold # Total Area over = under the y=x line
  estimated <- cVals # Estimated values
  theoretic <- seq_along(cVals)/unique(nfeatures) # y=x values
  
  if (plotIt) {
    plot(theoretic, estimated, ...)
    abline(0, 1)
  } else {
  }
  
  diffPVals <- theoretic - estimated
  indConserv <- theoretic <= estimated # Values over the y=x line
  conservArea <- sum(-diffPVals[indConserv])/MaxArea # Area over the y = x
  liberalArea <- sum(diffPVals[!indConserv])/MaxArea # Area under the y = x
  
  c(conservArea = conservArea, liberalArea = liberalArea, totalArea = liberalArea + 
      conservArea)
  
}

# Concordance plot function

gheat <- function(AUC_AOC_between_methods,concordance_df_summary){
  # Filtering
  gheat.list = list()
  for(m in unique(AUC_AOC_between_methods$.id)){
    AUC_AOC_between_methods.F <- AUC_AOC_between_methods[AUC_AOC_between_methods$.id == m,]
    concordance_df_summary.F <- concordance_df_summary[concordance_df_summary$.id == m,]
    forlegend <- AUC_AOC_between_methods.F
    forlegend$method1 <- factor(forlegend$method1,levels = c("DESeq2",
                                                             "DESeq2.ZI",
                                                             "DESeq2.apeglm",
                                                             "DESeq2.ZINBWave",
                                                             "edgeR",
                                                             "edgeR.robust",
                                                             "edgeR.ZINBWave",
                                                             "voom",
                                                             "EBSeq"), ordered = TRUE)
    
    g_legend_dendrogram <- get_legend(ggplot() + 
                                        geom_point(data=forlegend, aes(x = method1, y = 1, color = method1),size = 5) +
                                        scale_color_manual(values = cols) +
                                        theme_minimal() +
                                        theme(legend.position = "bottom") +
                                        guides(color = guide_legend(title = "Methods:",title.position = "left",nrow = 3)))
    
    # Clustering
    dist_matrix <- dcast(data = AUC_AOC_between_methods.F, formula = method1 ~ method2,value.var = "conservArea")
    dist_df <- dist_matrix[,2:ncol(dist_matrix)]
    rownames(dist_df) <- colnames(dist_df)
    distances <- as.dist(1-dist_df)
    hc <- hclust(d = distances)
    # Area extraction
    area <- apply(concordance_df_summary.F,1,function(x){
      area <- AUC_AOC_between_methods.F$conservArea[AUC_AOC_between_methods.F$method1 == x["method1"] & AUC_AOC_between_methods.F$method2 == x["method2"]]
      return(1-area)
    })
    concordance_df_summary.F_area <- cbind(concordance_df_summary.F,area = area)
    # As factor
    concordance_df_summary.F_area$method1 <- factor(concordance_df_summary.F_area$method1,levels = levels(concordance_df_summary.F_area$method1)[hc$order]) 
    
    concordance_df_summary.F_area$method2 <- factor(concordance_df_summary.F_area$method2,levels = levels(concordance_df_summary.F_area$method2)[hc$order]) 
    # edges
    edges <- data.frame(x = c(0,0,100,100),
                        xend = c(0,100,100,0),
                        y = c(0,1,1,0),
                        yend = c(1,1,0,0))
    # heatmap
    g_heat <- ggplot(concordance_df_summary.F_area,aes(x = rank, y = concordance)) +
      #geom_line(size = 1) +
      facet_grid(method1 ~ method2,scales = "free_x",switch = "y") +
      xlab("Rank") + # ylab("Concordance") +
      theme_pubr() + 
      theme(axis.text = element_blank(),
            #axis.text.x = element_text(hjust = 1, angle = 45),
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.line.y.right = element_blank(),
            # strip.text = element_text(hjust = 100, vjust = 100),
            # strip.background = element_rect(fill = "gray",linetype = 1,color = "white")) + 
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.spacing = unit(0,"cm"),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
      #geom_abline(mapping = aes(intercept = 0,slope = 1/nfeatures),color = "red",lty = 2) +
      coord_cartesian(xlim = c(0,100), ylim = c(0,1)) +
      
      geom_ribbon(aes(ymin = rank/nfeatures, ymax = concordance, fill = area)) +
      geom_segment(concordance_df_summary.F_area[concordance_df_summary.F_area$method1 == concordance_df_summary.F_area$method2,],mapping = aes(x = 0, xend = 0, y = 0, yend = 1, color = "red")) +
      
      geom_segment(concordance_df_summary.F_area[concordance_df_summary.F_area$method1 == concordance_df_summary.F_area$method2,],mapping = aes(x = 100, xend = 100, y = 1, yend = 0, color = "red")) +
      geom_segment(concordance_df_summary.F_area[concordance_df_summary.F_area$method1 == concordance_df_summary.F_area$method2,],mapping = aes(x = 0, xend = 100, y = 1, yend = 1, color = "red")) +
      geom_segment(concordance_df_summary.F_area[concordance_df_summary.F_area$method1 == concordance_df_summary.F_area$method2,],mapping = aes(x = 100, xend = 0, y = 0, yend = 0, color = "red")) +
      #scale_fill_gradientn(colours = c("red","yellow","turquoise"),limits = c(-0.01,1)) +
      scale_fill_distiller(palette = "RdYlBu",limits = c(-0.01,1)) +
      #scale_color_gradientn(colours = c("red","yellow","turquoise"),limits = c(-0.1,1)) +
      scale_y_continuous(breaks = c(0,0.5,1),position = "right") +
      scale_x_continuous(breaks = c(0,50,100))
    
    g_vertical_dendrogram <- ggplot() + 
      geom_segment(data=dendro_data(hc)$segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_label(data=dendro_data(hc)$labels, aes(x=x, y=y, label=label, hjust=1,color=label), nudge_y = 0) +
      coord_flip() + scale_y_reverse(expand = c(0,0,0,0)) + scale_x_reverse() +
      scale_color_manual(values = cols) +
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank(),
            legend.position = "none",
            panel.spacing = unit(0, "lines"),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
      ggtitle(label = paste0("Detection method: ",m),
              subtitle = "Concordance heatmap")
    
    g_horizontal_dendrogram <- ggplot() + 
      geom_segment(data=dendro_data(hc)$segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_point(data=dendro_data(hc)$labels, aes(x=x, y=y,color=label),size = 5) +
      scale_y_continuous() +
      #scale_y_reverse(expand=c(2,1)) + scale_x_reverse(expand=c(2,1)) +
      scale_color_manual(values = cols) +
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank(),
            legend.position = "none",
            panel.spacing = unit(0, "lines"),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
    
    addline_format <- function(x,...){
      gsub(':\\s',':\n',x)
    }
    g_heat_w_legend <- get_legend(ggplot(concordance_df_summary.F_area,aes(x = rank, y = concordance)) +
                                    facet_grid(method1 ~ method2,scales = "free_x",switch = "y") +
                                    labs(fill = addline_format("Rescaled Area from Rank 1 to 100 between: bisector and concordance")) +
                                    theme_minimal() + 
                                    theme(legend.position = "bottom") +
                                    guides(fill = guide_colorbar(title.position = "top",barwidth = 15)) +
                                    geom_ribbon(aes(ymin = rank/nfeatures, ymax = concordance, fill = area),alpha = 0.8) +
                                    scale_fill_distiller(palette = "RdYlBu",limits = c(0,1)) +
                                    scale_y_continuous(breaks = c(0,0.5,1),position = "right") +
                                    scale_x_continuous(breaks = c(0,50,100)))
    
    
    
    a <- plot_grid(plotlist = list(g_horizontal_dendrogram,g_horizontal_dendrogram,
                                   g_vertical_dendrogram,g_heat),align = "hv",axis = "lrtb")
    b <- g_heat_w_legend
    c <- g_legend_dendrogram
    gheat.list[[m]] <- list(plot = a,legend_heat = b, legend_dendro = c)
  }
  return(gheat.list)
}


# div = diversity: high, mid, low
# tech = data type: 16S or WMS
g_AUC <- function(conc_df,div,tech, rank = 100){
  conc_df_sub <- conc_df[conc_df$rank == rank & conc_df$subset == "1vs2",]
  case <- conc_df_sub[conc_df_sub$tech == tech & conc_df_sub$diversity == div,]
  ord <- order(ddply(case,.variables = ~ method1, function(x) median(x[,"concordance"]))$V1)
  ggplot(case,aes(x = method1, y = concordance, color = method1)) +
    geom_boxplot() +
    coord_flip() +
    scale_x_discrete(limits = unique(case$method1)[rev(ord)]) +
    xlab("Method") + ylab("Concordance") +
    ggtitle(label = paste(strsplit(as.character(unique(case$comp)),split = "_")[[1]],collapse = " vs "),
            subtitle = paste(unique(case$tech),"-",unique(case$diversity),"diversity")) +
    theme_minimal() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none",
          panel.spacing = unit(1,"lines")) +
    scale_color_manual(values = cols) +
    scale_y_continuous(limits = c(0,1))
}