library(zinbwave)
# library(EDASeq)
library(edgeR)
# library(MGLM)
library(plyr)

### Extraction of ZINB coefs
computeExp <- function(zinbModel){
  (1 - t(getPi(zinbModel))) * t(getMu(zinbModel))
}
computeVar <- function(zinbModel){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zinbModel))
  (1 - pi) * mu * (1 + mu*(phi + pi))
}
computeP0 <- function(zinbModel){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zinbModel))
  pi + (1 - pi) * (1 + phi * mu) ^ (-1/phi)
}

zinb.loglik.matrix <- function(model, x) {
  mu <- getMu(model)
  theta <- getTheta(model)
  theta_mat <- matrix(rep(theta, each = nrow(x)), ncol = ncol(x))
  pi <- getPi(model)
  lik <- pi * (x == 0) + (1 - pi) * dnbinom(x, size = theta_mat, mu = mu)
  lik[lik == 0] <- min(lik[lik != 0]) #to avoid log lik to be infinite
  log(lik)
}

fitZINB <- function(counts)
{
  cat("Model: Zero-Inflated Negative Binomial \n")
  fit <- zinbwave::zinbFit(Y = counts[c(1:20),],
                           epsilon = 1e10, # Regularization parameter fixed
                           commondispersion = TRUE,
                           verbose = TRUE,
                           BPPARAM = BiocParallel::SerialParam())
  mu = t(zinbwave::getMu(fit))
  pi = t(zinbwave::getPi(fit))
  phi = zinbwave::getPhi(fit)

  Y = log1p(rowMeans((1 - pi) * mu))
  Y0 = rowMeans(pi + (1 - pi) * (1 + phi * mu) ^ (-1/phi))

  
  theta = getTheta(fit)
  logipi = getLogitPi(fit)
  zinbllik = c()

  llik <- lapply(1:nrow(counts), function(i){
    y = counts[i,]
    mu.t = t(mu)[,i]
    zinbllik = zinb.loglik(Y = y, mu = mu.t, theta = theta[i], logitPi = logipi[,i])
    return(zinbllik)
    })
  
  zinbwave.llik = unlist(llik)
  
  return(data.frame("Y" = Y, "Y0" = Y0, "zinbloglik" = zinbwave.llik, row.names = names(Y)))
}

# Negative Binomial fitting
fitNB <- function(counts){
  cat("Model: Negative Binomial \n")
  # Default normalization
  normFacts <- edgeR::calcNormFactors(counts)
  # DGEList object creation
  dge <- edgeR::DGEList(counts = counts, norm.factors = normFacts)
  # Dispersion estimate
  disp <- edgeR::estimateDisp(y = dge, tagwise = TRUE)
  # GLM
  fit <- edgeR::glmFit(dge$counts,
                       dispersion = disp$tagwise.dispersion)
  nbloglik <- rowSums(dnbinom(x = counts,
                              size=1/disp$tagwise.dispersion,
                              mu=rowMeans(fit$fitted.values),
                              log = TRUE))
  # Fitted values extraction
  # Return the log(average fitted values + 1) to
  Y = log1p(rowMeans(fit$fitted.values))
  Y0 = rowMeans((1 + fit$fitted.values * disp$tagwise.dispersion)^(-1/disp$tagwise.dispersion))
  return(data.frame("Y" = Y, "Y0" = Y0, "nbloglik" = nbloglik, row.names = names(Y)))
}

# Negative Binomial fitting
fitNB_TMM <- function(counts, design){
  normFacts <- edgeR::calcNormFactors(counts, method = "TMM")
  dge <- DGEList(counts = counts)
  dge$samples$norm.factors <- normFacts
  disp <- estimateDisp(y = dge,design = design,tagwise = TRUE)
  fit <- glmFit(dge$counts,dispersion = disp$tagwise.dispersion, design = design)
  list(fitted = fit$fitted.values, disp = disp$tagwise.dispersion)
}

prepareObserved <- function(counts,
                            scale = NULL){
  if(!is.null(scale)){
    if(scale == "median"){
      counts <- counts * stats::median(colSums(counts)) / colSums(counts)
    } else if(scale == "default"){
      counts <- counts * 1e6 / colSums(counts)
    } else stop("When specified, 'scale' must be 'median' or 'default'")
  }
  Y <- log1p(rowMeans(counts))
  Y0 <- rowMeans(counts == 0)
  return(data.frame("Y" = Y, "Y0" = Y0))
}

meanDifferences <- function(estimated,
                            observed){
  if(nrow(estimated) != nrow(observed))
    stop("Estimated and Observed data.frames have different number of rows")
  MD <- estimated$Y - observed$Y
  ZPD <- estimated$Y0 - observed$Y0
  return(data.frame("MD" = MD, "ZPD" = ZPD))
}

RMSE <- function(differences){
  sqrt(mean(differences^2,na.rm = TRUE))
}

fitModels <- function(counts,
                      models = c("NB","ZINB"), 
                      scale = NULL, design=NULL){
  
  keep = raster::rowSums(x = counts) >= 2
  counts = counts[keep, ]
  nCovariates = 0
  if (!is.null(design)) {
    covariates = as.matrix(design)
    p = dim(design)[2]
    if(p==1 && length(unique(design)) == 1){
      # only intercept
      nCovariates = 1
    } else {
      nCovariates = p + 1
    }
  }
  
  fittedModels <- list()
  observed <- prepareObserved(counts, scale = scale)
  if("NB" %in% models){
    fitted <- fitNB(counts)
    MD <- meanDifferences(estimated = fitted,
                          observed = observed)
    aicNB = -2 * fitted$nbloglik + 2 * (nCovariates + 1)
    fittedModels$NB <- data.frame(observed,MD,EY=fitted$Y,EY0=fitted$Y0,loglik=fitted$nbloglik, aic = aicNB)
  }
  if("ZINB" %in% models){
    # fitted <- fitZINB(counts)

    by.sample <- DGEList(counts)
    # Normalizing using deseq2 method
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = ceiling(counts[,rownames(colData)]),
                                  colData = colData,
                                  design = design)
    dds <- estimateSizeFactors(dds, type = "poscounts")
    sf <- sizeFactors(dds)
    # Estimate the NB dispersion using deseq2 method
    dds <- estimateDispersions(dds, fitType = "local",
                               minmu = 1e-6)
    dds = DESeq(dds)
    # Estimate the log-overall mean
    centered.off <- getOffset(by.sample)
    centered.off <- centered.off - mean(centered.off)
    logmeans <- mglmOneGroup(by.sample$counts, offset = centered.off, dispersion = dispersions(dds))
    dispersion = dispersions(dds)
    # Estimate of dispersion with ZINB model
    zinb.prop <- rep(-Inf, nrow(counts))
    zinb.disp <- dispersions(dds)
    zinb.mean <- exp(logmeans)
    zinb.loglik <- rep(NA, nrow(counts))
    
    library(pscl)
    for (i in 1:nrow(counts)){
      # i=1
      if(sum(by.sample$counts[i,]==0)>0){
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
          # zinb.mean[i] <- mean(zfit$fitted.values)
          zinb.prop[i] <- zfit$coefficients$zero
          zinb.disp[i] <- 1/zfit$theta
          zinb.loglik[i] <- zfit$loglik
        })  
      } else {
        zfit <- glm.nb(by.sample$counts[i,] ~ 1)
        
        zinb.mean[i] <- mean(exp(zfit$coefficients))
        # zinb.mean[i] <- mean(zfit$fitted.values)
        zinb.prop[i] <- -Inf
        zinb.disp[i] <- 1/zfit$theta
        zinb.loglik[i] <- zfit$twologlik/2
      }
      
    }
    zinb.prop <- exp(zinb.prop)/(1+exp(zinb.prop))

    EY <- log1p((1 - zinb.prop) * zinb.mean)
    
    EY0 = zinb.prop + (1 - zinb.prop) * (1 + zinb.disp * zinb.mean) ^ (-1/zinb.disp)
    
    
    aicZINB = -2 * zinb.loglik + 2 * (nCovariates + 1)
    fitted.zinb = data.frame("Y" = EY, "Y0" = EY0)
    MD <- meanDifferences(estimated = fitted.zinb,
                          observed = observed)
    fittedModels$ZINB <- data.frame(observed,MD,EY=EY,EY0=zinb.prop,loglik=zinb.loglik, aic = aicZINB)
  }
  return(fittedModels)
}

extract_values <- function(site,distributions = c("ZINB","NB"), varnames = c("EY","EY0")){
  cbind(ldply(lapply(site[distributions],function(dist){
    as.data.frame(dist[varnames])
  }),.id = "Models"),
  ldply(lapply(site["OBSERVED"],function(dist){
    observed_df <- as.data.frame(dist)
    return(observed_df)
  })))
}

# From model list to data.frame
model_to_data.frame <- function(model_list){
  # df_list <- lapply(model_list,extract_values)
  # df <- ldply(df_list)
  df_list = lapply(model_list, function(x) {
    temp = rbindlist(x, idcol = "Models")
    temp$circ_id = c(rownames(x$NB),rownames(x$ZINB))
    temp
  })
  df = ldply(df_list)
  colnames(df) <- c("Method", "Model","Y","Y0","MD", "ZPD", "EY","EY0", "loglik", "AIC", "circ_id")
  return(df)
}


# GOF indexes
# Root mean square errors
compute_RMSE <- function(df){
  ddply(df,.variables = ~ Models + Methods, function(x){
    MD <- sqrt(mean(x$mean_diff^2,na.rm = TRUE))
    # MD <- sqrt(x$mean_diff^2)
    MD_sd <- sd(MD)
    ZPD <- sqrt(mean(x$zero_prob^2,na.rm = TRUE))
    # ZPD <- sqrt(abs(x$zero_prob)^2)
    ZPD_sd <- sd(ZPD)
    return(data.frame("MD" = MD, "ZPD" = ZPD))
  })
}

compute_MAE <- function(df){
  ddply(df,.variables = ~ Models + Methods, function(x){
    MD <- mean(x$mean_diff,na.rm = TRUE)
    # MD <- sqrt(x$mean_diff^2)
    MD_sd <- sd(MD)
    ZPD <- mean(x$zero_prob,na.rm = TRUE)
    # ZPD <- sqrt(abs(x$zero_prob)^2)
    ZPD_sd <- sd(ZPD)
    return(data.frame("MD" = MD, "ZPD" = ZPD))
  })
}

RMSE.func <- function(differences){
  sqrt(mean(differences^2,na.rm = TRUE))
}


MDPlot <- function(data, difference = NULL, split = TRUE, method_cols){
  if(difference == "MD"){
    if("Model" %in% colnames(data) & 
       "Method" %in% colnames(data) &
       "Y" %in% colnames(data) &
       "MD" %in% colnames(data)){
      
      RMSE_MD <- plyr::ddply(.data = data,
                             .variables = ~ Method + Model,
                             .fun = function(m) cbind("RMSE" = RMSE.func(m$MD)))
      
      gobj <- ggplot(data = data, aes(x = Y, y = MD, color = Model)) +
        scale_color_manual(values = method_cols) +
        ggtitle(label = "Mean Differences plot",
                subtitle = paste0("Observed = log(mean(counts*)+1)",
                                  "\n",
                                  "Estimated = log(mean(fitted*)+1)"))
      
      if(split){
        gobj <- gobj +
          geom_text(data = RMSE_MD, color = "black",
                    aes(x = mean(data$Y)+2,
                        y = max(data$MD,na.rm = TRUE),
                        label = paste0("RMSE:",round(RMSE,3))), size = 3)
      }
      
    } else {stop("data should contains 'Model', 'Y', and 'MD' columns for model
              name, observed values and mean difference values respectively.")}
  } else if(difference == "ZPD"){
    if("Model" %in% colnames(data) &
       "Method" %in% colnames(data) &
       "Y0" %in% colnames(data) &
       "ZPD" %in% colnames(data)){
      
      RMSE_ZPD <- ddply(.data = data,
                              .variables = ~ Method + Model,
                              .fun = function(m) cbind("RMSE" = RMSE.func(m$ZPD)))
      
      gobj <- ggplot(data = data, aes(x = Y0,  y = ZPD,  color = Model)) +
        scale_color_manual(values = method_cols) +
        ggtitle(label = "Zero Probability Differences plot", subtitle =
                  "Observed = mean(counts=0)\nEstimated = mean(P(Y=0))")
      
      if(split){
        gobj <- gobj +
          geom_text(data = RMSE_ZPD, color = "black",
                    aes(x = 0.25,
                        y = max(data$ZPD,na.rm = TRUE),
                        label = paste0("RMSE:",round(RMSE,2))), size = 3)
      }
      
    } else {stop("df should contains 'Model', 'Y0', and 'ZPD' columns for model
              name, zero rate observed values and zero probability difference
              values respectively.")}
  } else stop("Difference must be 'MD' or 'ZPD'")
  
  gobj <- gobj +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    theme(legend.position = "bottom") +
    xlab("Observed") +
    ylab("Estimated-Observed")
  
  if(length(unique(data$Model))>1){
    if(split){
      gobj <- gobj +
        facet_grid(Model ~ Method, labeller = label_wrap_gen(width=10)) +
        geom_point(pch = 21) +
        geom_smooth(color = "black") +
        theme_classic() +
        theme(legend.position = "bottom",
              strip.text.x = element_text(face = "bold.italic", size = 7.5),
              strip.text.y = element_text(face = "bold.italic", size = 10),
              strip.background = element_rect(colour = "grey", size = 1),
              axis.text.x = element_text(hjust = 1, angle = 45, size = 14),
              axis.text.y = element_text(size = 14),
              text = element_text(size=14),
              plot.margin = unit(c(0,0,0,0), "cm"))
    } else {
      gobj <- gobj +
        geom_smooth()
    }
  }
  
  return(gobj)
}


RMSEPlot <- function(data, difference = NULL, method_cols){
  if(difference == "MD"){
    if("Model" %in% colnames(data) & "Y" %in% colnames(data) &
       "MD" %in% colnames(data)){
      
      RMSE <- plyr::ddply(.data = data,
                          .variables = ~ Model + Method,
                          .fun = function(m) cbind("RMSE" = RMSE.func(m$MD)))
      
      gobj <- ggplot(data = RMSE, aes(x = Model, y = RMSE, fill = Model)) +
        scale_fill_manual(values=method_cols) +
        geom_col() +
        geom_label(aes(x = Model, y = RMSE, label = round(RMSE,2)),
                   fill = "white") +
        ggtitle(label = "RMSE",
                subtitle = "Mean differences")
      
    } else {stop("data should contains 'Model', 'Y', and 'MD' columns for model
              name, observed values and mean difference values respectively.")}
  } else if(difference == "ZPD"){
    if("Model" %in% colnames(data) &
       "Y0" %in% colnames(data) &
       "ZPD" %in% colnames(data)){
      
      RMSE <- plyr::ddply(.data = data,
                          .variables = ~ Model + Method,
                          .fun = function(m) cbind("RMSE" = RMSE.func(m$ZPD)))
      
      gobj <- ggplot(data = RMSE, aes(x = Model, y = RMSE, fill = Model)) +
        geom_col() +
        scale_fill_manual(values=method_cols) +
        geom_label(aes(x = Model, y = RMSE, label = round(RMSE,4)),
                   fill = "white") +
        ggtitle(label = "RMSE",
                subtitle = "Zero probability difference")
      
    } else {stop("df should contains 'Model', 'Y0', and 'ZPD' columns for model
              name, zero rate observed values and zero probability difference
              values respectively.")}
    
  } else stop("Difference must be 'MD' or 'ZPD'")
  
  gobj <- gobj +
    facet_grid( ~ Method, labeller = labeller(.cols = label_both)) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_discrete() +
    theme_classic() +
    theme(legend.position = "bottom",
          strip.text.x = element_text(face = "bold.italic", size = 5),
          strip.text.y = element_text(face = "bold.italic", size = 8), #angle = 75),
          strip.background = element_rect(#fill = "lightblue",
            colour = "grey", size = 1),
          axis.text.x = element_text(hjust = 1, angle = 45),
          # axis.text.x = element_blank(),
          # axis.text.y = element_blank(),
          # panel.spacing = unit(0, "lines"),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  return(gobj)
}
