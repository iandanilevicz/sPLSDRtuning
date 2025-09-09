
# c-index 

sim_surv_c = function(lambda, n, M, eta, ncomp, nfold, maxt){
  ne = length(eta)
#  Eta = cbind(eta, rep(0, ne), rep(0, ne))
  Eta = matrix(0, ncol = 3, nrow = ne)
  Eta[1:ne,1] = eta
  count_max1  = rep(0, ne)
  count_max2  = rep(0, ne)
  cindex_eta1 = rep(0, ne)
  cindex_eta2 = matrix(0, ncol=nfold, nrow=ne)
  n_features  = matrix(0, ncol=ne, nrow=M)
  n_cases     = rep(0, M)
  for(m in 1:M){
    R  = matrix(c(1,   0.5, 0.5,
                  0.5, 1,   0.5,
                  0.5, 0.5, 1), 
                nrow = 3, ncol = 3, byrow = TRUE)
    mu1 = c(c(0.15, 0.1, 0.5))
    mu2 = c(c(-0.15, -0.1, -0.5))
    var1 = MASS::mvrnorm(n, mu = mu1, Sigma = R)
    var2 = MASS::mvrnorm(n, mu = mu2, Sigma = R)
    features = data.frame(feat1 = var1[,1], 
                          feat2 = var1[,2],
                          feat3 = var1[,3],
                          feat4 = var2[,1],
                          feat5 = var2[,2],
                          feat6 = var2[,3],
                          feat7 = stats::rnorm(n, 1, 0.5),
                          feat8 = stats::rnorm(n, 1, 0.5),
                          feat9 = stats::rnorm(n, 1, 0.5),
                          feat10 = stats::rnorm(n, 1, 0.5),
                          feat11 = stats::rnorm(n, 1, 0.5),
                          feat12 = stats::rnorm(n, 1, 0.5))
    sim = simsurv(lambdas = lambda,  betas = c(feat1 = 0.15, feat2=0.1, feat3=0.05, 
                                               feat4 = -0.15, feat5=-0.1, feat6=-0.05,
                                               feat7 = 0, feat8=0, feat9=0,
                                               feat10 = 0, feat11=0, feat12=0),
                  x = features, maxt = maxt, dist ="exponential")
    folds.fs = createFolds(y = sim$status, k = nfold)
    Y = cbind(sim$eventtime, sim$status); colnames(Y) = c("time", "event")
    n_cases[m] = sum(sim$status)
    for(et in 1:ne){
      mod_1 = splsdrcox_penalty(features, Y, n.comp = ncomp, penalty = eta[et], x.center = FALSE, x.scale = FALSE)
      cindex_eta1[et] = mod_1$survival_model$fit$concordance[6] # c-index 
      n_features[m, et] = dim(mod_1$X$W.star)[1]
      for(f in 1:nfold){
        subcases = folds.fs[[f]]
        Yf = Y[subcases,1:2] 
        features_f = features[subcases,1:12]
        mod_2 = splsdrcox_penalty(features_f, Yf, n.comp = ncomp, penalty = eta[et], x.center = FALSE, x.scale = FALSE)
        cindex_eta2[et, f] = mod_2$survival_model$fit$concordance[6] # c-index
       } 
    }
    pos_max1 = which.max(cindex_eta1)
    #cindex_eta2_max = rep(0, ne)
    cindex_eta2_max = rowMeans(cindex_eta2)
    #for(et in 1:ne){
    #  cindex_eta2_max[et] = max(cindex_eta2[et,])
    #}
    pos_max2 = which.max(cindex_eta2_max)
    count_max1 = count_max(eta, max1=pos_max1, count1=count_max1)
    count_max2 = count_max(eta, max1=pos_max2, count1=count_max2)
    print(paste0(m/M*100, " % of Monte Carlo done"))
  }
  Eta[,2] = count_max1
  Eta[,3] = count_max2
  colnames(Eta) = c("eta", "freq_1fold", "freq_nfold")
  #N_features = colMedians(n_features)
  #print(n_features)
  return(list(Eta=Eta, N_features = n_features, N_cases = n_cases))
}

sim_surv_auc = function(lambda, n, M, eta, ncomp, nfold, maxt){
  ne = length(eta)
  #  Eta = cbind(eta, rep(0, ne), rep(0, ne))
  Eta = matrix(0, ncol = 3, nrow = ne)
  Eta[1:ne,1] = eta
  count_max1  = rep(0, ne)
  count_max2  = rep(0, ne)
  cindex_eta1 = rep(0, ne)
  cindex_eta2 = matrix(0, ncol=nfold, nrow=ne)
  n_cases     = rep(0, M)
  for(m in 1:M){
    R  = matrix(c(1,   0.5, 0.5,
                  0.5, 1,   0.5,
                  0.5, 0.5, 1), 
                nrow = 3, ncol = 3, byrow = TRUE)
    mu1 = c(c(0.15, 0.1, 0.5))
    mu2 = c(c(-0.15, -0.1, -0.5))
    var1 = MASS::mvrnorm(n, mu = mu1, Sigma = R)
    var2 = MASS::mvrnorm(n, mu = mu2, Sigma = R)
    features = data.frame(feat1 = var1[,1], 
                          feat2 = var1[,2],
                          feat3 = var1[,3],
                          feat4 = var2[,1],
                          feat5 = var2[,2],
                          feat6 = var2[,3],
                          feat7 = stats::rnorm(n, 1, 0.5),
                          feat8 = stats::rnorm(n, 1, 0.5),
                          feat9 = stats::rnorm(n, 1, 0.5),
                          feat10 = stats::rnorm(n, 1, 0.5),
                          feat11 = stats::rnorm(n, 1, 0.5),
                          feat12 = stats::rnorm(n, 1, 0.5))
    sim = simsurv(lambdas = lambda,  betas = c(feat1 = 0.15, feat2=0.1, feat3=0.05, 
                                               feat4 = -0.15, feat5=-0.1, feat6=-0.05,
                                               feat7 = 0, feat8=0, feat9=0,
                                               feat10 = 0, feat11=0, feat12=0),
                  x = features, maxt = maxt, dist ="exponential")
    folds.fs = createFolds(y = sim$status, k = nfold)
    Y = cbind(sim$eventtime, sim$status); colnames(Y) = c("time", "event")
    n_cases[m] = sum(sim$status)
    for(et in 1:ne){
      mod_1 = splsdrcox_penalty(features, Y, n.comp = ncomp, penalty = eta[et], x.center = FALSE, x.scale = FALSE)
      yhat = predict(mod_1)
      cindex_eta1[et] = cenROC(Y=Y[,1], M=(yhat), censor=Y[,2], t=maxt) # AUC 
      for(f in 1:nfold){
        subcases = folds.fs[[f]]
        Yf = Y[subcases,1:2] 
        features_f = features[subcases,1:12]
        mod_2 = splsdrcox_penalty(features_f, Yf, n.comp = ncomp, penalty = eta[et], x.center = FALSE, x.scale = FALSE)
        cindex_eta2[et, f] = mod_2$survival_model$fit$concordance[6] # c-index
      } 
    }
    pos_max1 = which.max(cindex_eta1)
    #cindex_eta2_max = rep(0, ne)
    cindex_eta2_max = rowMeans(cindex_eta2)
    #for(et in 1:ne){
    #  cindex_eta2_max[et] = max(cindex_eta2[et,])
    #}
    pos_max2 = which.max(cindex_eta2_max)
    count_max1 = count_max(eta, max1=pos_max1, count1=count_max1)
    count_max2 = count_max(eta, max1=pos_max2, count1=count_max2)
    print(paste0(m/M*100, " % of Monte Carlo done"))
  }
  Eta[,2] = count_max1
  Eta[,3] = count_max2
  colnames(Eta) = c("eta", "freq_1fold", "freq_nfold")
  return(Eta)
}

# AIC  
sim_surv_aic = function(lambda, n, M, eta, ncomp, nfold, maxt){
  ne = length(eta)
  #  Eta = cbind(eta, rep(0, ne), rep(0, ne))
  Eta = matrix(0, ncol = 3, nrow = ne)
  Eta[1:ne,1] = eta
  count_max1  = rep(0, ne)
  count_max2  = rep(0, ne)
  cindex_eta1 = rep(0, ne)
  cindex_eta2 = matrix(0, ncol=nfold, nrow=ne)
  n_features  = matrix(0, ncol=ne, nrow=M)
  n_cases     = rep(0, M)
  for(m in 1:M){
    R  = matrix(c(1,   0.5, 0.5,
                  0.5, 1,   0.5,
                  0.5, 0.5, 1), 
                nrow = 3, ncol = 3, byrow = TRUE)
    mu1 = c(c(0.15, 0.1, 0.5))
    mu2 = c(c(-0.15, -0.1, -0.5))
    var1 = MASS::mvrnorm(n, mu = mu1, Sigma = R)
    var2 = MASS::mvrnorm(n, mu = mu2, Sigma = R)
    features = data.frame(feat1 = var1[,1], 
                          feat2 = var1[,2],
                          feat3 = var1[,3],
                          feat4 = var2[,1],
                          feat5 = var2[,2],
                          feat6 = var2[,3],
                          feat7 = stats::rnorm(n, 1, 0.5),
                          feat8 = stats::rnorm(n, 1, 0.5),
                          feat9 = stats::rnorm(n, 1, 0.5),
                          feat10 = stats::rnorm(n, 1, 0.5),
                          feat11 = stats::rnorm(n, 1, 0.5),
                          feat12 = stats::rnorm(n, 1, 0.5))
    sim = simsurv(lambdas = lambda,  betas = c(feat1 = 0.15, feat2=0.1, feat3=0.05, 
                                               feat4 = -0.15, feat5=-0.1, feat6=-0.05,
                                               feat7 = 0, feat8=0, feat9=0,
                                               feat10 = 0, feat11=0, feat12=0),
                  x = features, maxt = maxt, dist ="exponential")
    folds.fs = createFolds(y = sim$status, k = nfold)
    Y = cbind(sim$eventtime, sim$status); colnames(Y) = c("time", "event")
    n_cases[m] = sum(sim$status)
    for(et in 1:ne){
      mod_1 = splsdrcox_penalty(features, Y, n.comp = ncomp, penalty = eta[et], x.center = FALSE, x.scale = FALSE)
      cindex_eta1[et] = mod_1$survival_model$AIC # AIC 
      n_features[m, et] = dim(mod_1$X$W.star)[1]
      for(f in 1:nfold){
        subcases = folds.fs[[f]]
        Yf = Y[subcases,1:2] 
        features_f = features[subcases,1:12]
        mod_2 = splsdrcox_penalty(features_f, Yf, n.comp = ncomp, penalty = eta[et], x.center = FALSE, x.scale = FALSE)
        cindex_eta2[et, f] = mod_2$survival_model$AIC # AIC
      } 
    }
    pos_max1 = which.min(cindex_eta1)
    #cindex_eta2_max = rep(0, ne)
    cindex_eta2_max = rowMeans(cindex_eta2)
    #for(et in 1:ne){
    #  cindex_eta2_max[et] = max(cindex_eta2[et,])
    #}
    pos_max2 = which.min(cindex_eta2_max)
    count_max1 = count_max(eta, max1=pos_max1, count1=count_max1)
    count_max2 = count_max(eta, max1=pos_max2, count1=count_max2)
    print(paste0(m/M*100, " % of Monte Carlo done"))
  }
  Eta[,2] = count_max1
  Eta[,3] = count_max2
  colnames(Eta) = c("eta", "freq_1fold", "freq_nfold")
  return(list(Eta=Eta, N_features = n_features, N_cases = n_cases))
  return(Eta)
}

# aux 
count_max = function(eta, max1, count1){
  ne = length(eta)
  count1[max1] = count1[max1] + 1
  return(count1)
}



