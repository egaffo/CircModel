load(file = "/blackhole/alessia/CircModel/power/IPF_sensitivityPrecision_CCP2_glmglmm_30rep.RData")
namesAlgos = c("circMeta", "voom")
resTest_voom_cMeta0 <- list()
resHeldout_voom_cMeta0 <- list()
res_voom_cMeta <- bplapply(1:30, function(i) {   
  
  cat(i," ")
  # i = 1
  
  testSet <- as.character(randomSubsets[i,c(1:7)])
  heldOutSet <- as.character(randomSubsets[i,-c(1:7)])
  
  eTest <- e[,testSet]
  
  eHeldout <- e[,heldOutSet]
  
  resIDx = rownames(resTes[[i]])
  
  eTest = eTest[resIDx,]
  eHeldout = eHeldout[resIDx,]
  
  resTest_voom_cMeta0$circMeta = runPois.ztest(e = eTest)
  resTest_voom_cMeta0$voom = runVoom(e = eTest)

  resHeldout_voom_cMeta0$circMeta = runPois.ztest(e = eHeldout)
  resHeldout_voom_cMeta0$voom = runVoom(e = eHeldout)
  
  resTest_voom_cMeta <- as.data.frame(c(lapply(resTest_voom_cMeta0, function(z) z$padj[resIDx])))
  resHeldout_voom_cMeta <- as.data.frame(c(lapply(resHeldout_voom_cMeta0, function(z) z$padj[resIDx])))
  rownames(resTest_voom_cMeta) <- resIDx
  rownames(resHeldout_voom_cMeta) <- resIDx
  
  list(resTestNew=resTest_voom_cMeta,resHeldoutNew=resHeldout_voom_cMeta)
}, BPPARAM =  MulticoreParam(workers = 5))
resTesNew <- lapply(res_voom_cMeta, "[[", "resTestNew") #lapply(res, function(x) lapply(x, "[[", "resTest"))
resHeldoutNew <- lapply(res_voom_cMeta, "[[", "resHeldoutNew")

for(i in 1:30){
  resTes[[i]]$voom = resTesNew[[i]]$voom
  resHeldout[[i]]$voom = resHeldoutNew[[i]]$voom
  resTes[[i]]$circMeta = resTesNew[[i]]$circMeta
  resHeldout[[i]]$circMeta = resHeldoutNew[[i]]$circMeta
}
