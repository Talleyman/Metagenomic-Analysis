TwoStage_Package <- function(X,Y,fileoutput, TMM.option){
  require(glmnet)   # Elastic-net
  require(MASS)     # GLMs-NB
  require(edgeR)    # TMM normalization
  
  threshold = 0.05  # FDR to control the type I error at significance level of 0.05
  
  #####################################################
  ###   Trimmed Mean (TMM) (Based on edgeR Paper)   ###
  #####################################################
  TMMNorm = function(X, Y, TMM.option){
    factors = calcNormFactors(X,method="TMM") #Calculate normaization factors
    if(TMM.option==1){
      eff.lib.size = colSums(X)*factors;
      ref.lib.size = mean(eff.lib.size); #Use the mean of the effective library sizes as a reference library size 
      X.output = sweep(X,MARGIN=2,eff.lib.size,"/")*ref.lib.size; #Normalized read counts
    } else if(TMM.option==2){ #Regenerate counts with a common dispersion.
      group = as.factor(Y);
      X.tmp = DGEList(counts=X, group=group);
      X.tmp = calcNormFactors(X.tmp);
      X.tmp = estimateCommonDisp(X.tmp);
      X.output = as.matrix(X.tmp$pseudo.counts);
    } else {
      stop("TMM.option must be either 1 or 2!");
    }
    return(X.output);
  }
  
  ############################################################################
  #### Negative binomial generalized linear model for two group comparison ###
  ############################################################################
  get.glmNB <- function(count,category){
    p.val = rep(NA,ncol(count))  #Declare initial pvalue
    
    for (i in 1:ncol(count)) {
      dat = data.frame(testing=count[,i], groups=as.factor(category))
      if (sum(dat$testing) != 0) { #in case that all counts are zero
        fit = glm.nb(testing ~ groups, data=dat)
        p.val[i] = summary(fit)$coefficients[2,4] 
      }
    }
    p.val[is.na(p.val)] = 1
    p.adj = p.adjust(p.val, method="fdr") #Adjusted p-value 
    res = cbind(p.val,p.adj)
    rownames(res) = colnames(count)
    return(res)  
  }
  
  #################################################################################
  #### Negative binomial generalized linear model for multiple group comparison ###
  #################################################################################
  get.multiglmNB <- function(count,df){
    
    category = as.character(df$Category)
    pheno = as.character(df$Phenotype)
    
    # A vector contains overall p-values
    p.overall = rep(NA, ncol(count)) 
    
    # Number of beta for hypothesis testing 
    nbeta.testing = nlevels(factor(category)) - 1  
    
    for (nb in 1:nbeta.testing){
      
      if (nb != 1){  
        base.old = which(category==1)
        base.new = which(category==nb)
        category[base.old] = nb
        category[base.new] = 1
      } 
      
      ref.group = unique(pheno[which(category==1)])
      other.group = unique(pheno[which(category!=1)])
      
      # Declare a matrix containg p-values
      p.val = matrix(NA,nrow=ncol(count), ncol=nbeta.testing) # A matrix contains all p-values for pairwise comparisons
      colnames(p.val) = paste(ref.group, "vs", other.group, sep="")
      rownames(p.val) = colnames(count)
      
      # Comparison with the base ref. group
      for (i in 1:ncol(count)) {
        dat = data.frame(testing=count[,i], groups=as.factor(category))
        if (sum(dat$testing) != 0) { #in case that all counts are zero
          fit = glm.nb(testing ~ groups, data=dat)
          if (nb == 1){
            fit0 = glm.nb(testing ~ 1, data=dat)
            p.overall[i] = anova(fit0,fit)$"Pr(Chi)"[2]  
          }
          p.val[i,] = summary(fit)$coefficients[-1,4] 
        }
      }
      
      if (nb == 1){
        res.pval = p.val
      } else{
        res.pval = cbind(res.pval, p.val)
      }
      
    } 
    
    # A matrix containg adjusted p-values
    res.padj = apply(res.pval, 2, function(x) p.adjust(x, method="fdr"))
    res.poverall.adj = p.adjust(p.overall, method="fdr")
    res = cbind(res.poverall.adj, res.padj) 
    colnames(res) = c("p.overall", colnames(res.pval))
    
    return(res)  
  }
  
  #########################################################
  #### Get model coefficients for glmnet-multinomial   ####
  #########################################################
  fgetcoef = function(coef){
    
    for (len in 2:length(coef)){
      temp = as.matrix(coef)[[len]]
      if (len == 2){
        mtx = as.matrix(temp[-1,])   #get rid of intercept 
      } else{
        temp = as.matrix(temp[-1,])  #get rid of intercept
        mtx = cbind(mtx, temp)
      }
    }
    
    return(mtx)
  }
  
  ##############################################################################
  ###                           Two-Stage Procedure                          ###
  ##############################################################################
  #################################
  ###   Two group comparison   ####
  #################################
  Twostage <- function(X, Y, fileoutput){
    ### A matrix contains cvm error information
    cv.mat = matrix(NA, nrow = length(alpha), ncol = 3)
    colnames(cv.mat) = c("Alpha", "Lambda", "CVM")
    cv.mat[,1] = alpha
    
    #### Cross-validation to determine the optimal alpha and lambda giving minimum cross-validation error 
    npop1 = table(Y)[1]
    npop2 = table(Y)[2]
    N = npop1 + npop2
    if (N <= 15) nfold = 3
    if (N > 15  & N < 50) nfold = 5
    if (N >= 50) nfold = 10  
    
    ### Fit a GLM with elastic net regularization
    for (iter in 1:3){
      set.seed(iter)
      
      #Select a random number of subjects
      newdata.ind = sample(1:nrow(X), floor(0.9*nrow(X)))
      
      #Subset the data
      X.new = X[newdata.ind,]
      Y.new = Y[newdata.ind]
      
      # Cross-varidation to determine lambda
      for (j in 1:length(alpha)){
        cv = cv.glmnet(X, Y, family = "binomial", alpha = alpha[j], nfolds=nfold)
        ind = match(cv$lambda.min, cv$lambda)
        cv.mat[j,2] = cv$lambda[ind]
        cv.mat[j,3] = cv$cvm[ind]
      }
      
      alpha.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),1]
      lambda.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),2]
      
      # Fit a GLM with elastic net regularization
      fit = glmnet(X.new, Y.new, family = "binomial", alpha = alpha.opt, lambda = lambda.opt)
      
      # Get model coefficients for glmnet
      coef = coef(fit)  
      
      if (iter == 1) {
        elast = as.matrix(coef)
      } else{
        elast = cbind(elast, as.matrix(coef))
      }
    } #end (iter)
    
    ### Get features for which coefficients are not 0
    elast = elast[-1, ] #get rid of intercept
    feature = rownames(elast)
    df = data.frame(elast, rowsum = rowSums(elast))
    ind.elast = which(df$rowsum !=0 ) #index of coefficients that are not zero 
    allFeatSel = as.character(feature[ind.elast])
    
    if (length(allFeatSel) != 0){
      ind = which(colnames(X) %in% allFeatSel)
      X.elast = as.matrix(X[,ind])
      colnames(X.elast) = allFeatSel
      # GLMs-NB Model
      pvalue.mat = get.glmNB(X.elast,Y) 
    } else{
      print("No significantly differentially abundant features!!!")
    }
    
    ### Summary of significantly differentially abundant features
    sig.mat = pvalue.mat[pvalue.mat[,"p.adj"] < threshold,] # Significantly differentially abundant features 
    sig = rownames(sig.mat)
    ind.sig = which(colnames(X) %in% sig)
    X.sig = t(as.matrix(X[,ind.sig]))
    X.sig.pop1 = X.sig[,1:npop1]
    X.sig.pop2 = X.sig[,(npop1+1):(npop1+npop2)]
    mean.pop1 = apply(X.sig.pop1, 1, function(x) mean(x))
    mean.pop2 = apply(X.sig.pop2, 1, function(x) mean(x))
    sd.pop1 = apply(X.sig.pop1, 1, function(x) sd(x))
    sd.pop2 = apply(X.sig.pop2, 1, function(x) sd(x))
    sig.df = data.frame(Annotation=as.character(rownames(X.sig)), 
                        mean_group1=mean.pop1, sd_group1=sd.pop1,
                        mean_group2=mean.pop2, sd_group2=sd.pop2,
                        p.val=sig.mat[,"p.val"], p.adj=sig.mat[,"p.adj"], stringsAsFactors=F)
    rownames(sig.df) = NULL
    sigsort.df = sig.df[order(sig.df$p.adj),] 
    write.csv(sigsort.df, file = fileoutput, row.names=F)
  }
  
  ######################################
  ###   Multiple group comparison   ####
  ######################################
  Twostage.multicomp <- function(X, Y, fileoutput, ngroups){
    ### A matrix contains cvm error information
    cv.mat = matrix(NA, nrow = length(alpha), ncol = 3)
    colnames(cv.mat) = c("Alpha", "Lambda", "CVM")
    cv.mat[,1] = alpha
    
    #### Cross-validation to determine the optimal alpha and lambda giving minimum cross-validation error 
    npop = as.integer(table(Y))    
    N = sum(npop)
    if (N <= 15) nfold = 3
    if (N > 15  & N < 50) nfold = 5
    if (N >= 50) nfold = 10  
    
    ### Fit a GLM with elastic net regularization
    for (iter in 1:3){
      set.seed(iter)
      
      #Select a random number of subjects
      newdata.ind = sample(1:nrow(X), floor(0.9*nrow(X)))
      
      #Subset the data
      X.new = X[newdata.ind,]
      Y.new = Y[newdata.ind]
      
      # Cross-validation to determine lambda
      for (j in 1:length(alpha)){
        cv = cv.glmnet(X, Y, family="multinomial", alpha=alpha[j], nfolds=nfold)
        ind = match(cv$lambda.min, cv$lambda)
        cv.mat[j,2] = cv$lambda[ind]
        cv.mat[j,3] = cv$cvm[ind]
      }
      
      alpha.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),1]
      lambda.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),2]
      
      # Fit a GLM with elastic net regularization
      fit = glmnet(X.new, Y.new, family = "multinomial", alpha = alpha.opt, lambda = lambda.opt)
      
      # Get model coefficients for glmnet
      coef = coef(fit)  
      
      if (iter == 1) elast1 = fgetcoef(coef)    
      if (iter == 2) elast2 = cbind(elast1, fgetcoef(coef))
      if (iter == 3) elast = cbind(elast2, fgetcoef(coef))
      
    } #end (iter)
    
    # Get features for which coefficients are not 0
    feature = rownames(elast)
    df = data.frame(elast, rowsum = rowSums(elast))
    ind.elast = which(df$rowsum !=0 ) #index of coefficients that are not zero 
    allFeatSel = as.character(feature[ind.elast])
    
    if (length(allFeatSel) != 0){
      ind = which(colnames(X) %in% allFeatSel)
      X.elast = as.matrix(X[,ind])
      colnames(X.elast) = allFeatSel
      
      # GLMs-NB Model
      Y.new = rep(NA, length(Y))
      for (ilevel in 1:nlevels(Y)){
        ind.level = which(Y == levels(Y)[ilevel])
        Y.new[ind.level] = ilevel
      }
      Y.df = data.frame(Phenotype=Y, Category=Y.new)
      df = Y.df[order(Y.df$Category),]
      pvalue.mat = get.multiglmNB(X.elast,df) 
      
    } else{
      sig = NA 
      print("No differentially abundant features!!!")
    }
    
    ### Summary of significantly differentially abundant features
    sig.mat0 = pvalue.mat[ pvalue.mat[,"p.overall"] < threshold , ]
    sig.mat = sig.mat0[apply(sig.mat0[,-1], 1, function(x) any(x < threshold)),]
    sig.sort = sig.mat[order(sig.mat[,"p.overall"]),] 
    
    fcutcol = function(ngroups){
      ntimes = ngroups - 2
      for (itime in 1:ntimes){
        ncuts = itime
        cutcols = rep(NA, ncuts)
        for (jcut in 1:ncuts){
          cutcols[jcut] = (itime * ngroups) - (jcut - 1)
        }
        if (itime == 1){
          res = cutcols
        } else{
          res = c(res,cutcols)
        }
      }
      return(res)
    }
    
    cutcols = fcutcol(ngroups) + 1 
    newsub = sig.sort[,-c(cutcols)]
    df = data.frame(Annotation=rownames(sig.sort), newsub)
    write.csv(df, file = fileoutput, row.names=F)
  }
  
  alpha = seq(0.01,0.1,by=0.01)
  
  ### TMM Normalization ###
  rownames(X) = X[,1]
  X = X[,-1]
  X = as.matrix(X)
  Y = factor(Y[,2])
  X.tmm = TMMNorm(X,Y,TMM.option)  # TMM normalization (edgeR)
  X = t(X.tmm)
  
  ngroups = nlevels(Y)
  
  if (ngroups == 2) Twostage(X, Y, fileoutput) # Two group comparison
  if (ngroups > 2) Twostage.multicomp(X, Y, fileoutput, ngroups) # Multiple group comparison
}