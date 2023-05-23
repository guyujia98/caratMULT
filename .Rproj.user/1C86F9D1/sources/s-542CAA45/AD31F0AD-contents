
##---mean within stratum---##

meanS <- function(data,S){
  S_ls <- sort(unique(S))
  nS <- length(S_ls)
  n_col <- ncol(data)
  mean_s <- matrix(rep(0,nS*n_col),nS,n_col)
  for(s in S_ls){
    ind_s <- which(S==s)
    n_s <- length(S[ind_s])
    mean_s[s,] <- mean(data[ind_s,])
  }
  
  return(mean_s)
}

##---estimator for stratum common tau---##

tau_sc <- function(Y,X,S,A){
  nA <- length(unique(A))
  nS <- length(unique(S))
  p <- ncol(X)
  result <- rep(0,nA)
  pn_s <- rep(0,nS)
  for(s in sort(unique(S))){
    pn_s[s] <- length(which(S==s))/length(S)
  }
  meanXs <- meanS(X,S)
  for(i in 1:nA){
    ind <- which(A==(i-1))
    Ya <- Y[ind]
    Xa <- matrix(X[ind,],length(ind),p)
    Sa <- S[ind]
    beta_a <- coef_cal_sc(Ya,Sa,Xa)
    meanYas <- meanS(matrix(Ya),Sa)
    meanXas <- meanS(Xa,Sa)
    result[i] <- sum(pn_s*(meanYas - (meanXas - meanXs)%*%beta_a))
  }
  
  return(result)
}

##---estimator for stratum common sd, tau_a - tau_b---##

sd_sc_ne <- function(Y,X,S,A,a,b){
  sd_r <- 0
  sd_x <- 0
  sd_hy <- 0
  
  
  p <- ncol(X)
  sd_rx <- rep(0,p)
  
  ind_a <- which(A==a)
  ind_b <- which(A==b)
  
  pS <- c(table(S)/length(S))
  
  pi_as <- c(table(S[ind_a])/length(S[ind_a]))
  pi_bs <- c(table(S[ind_b])/length(S[ind_b]))
  
  beta_a <- coef_cal_sc_new(Y[ind_a],S[ind_a],matrix(X[ind_a,],length(ind_a),p),
                            pS,pi_as)
  beta_b <- coef_cal_sc_new(Y[ind_b],S[ind_b],matrix(X[ind_b,],length(ind_b),p),
                            pS, pi_bs)
  
  r_a <- Y[ind_a] - X[ind_a,] %*% beta_a
  r_b <- Y[ind_b] - X[ind_b,] %*% beta_b
  
  
  S_ls <- sort(unique(S))
  
  meanr_as <- meanS(r_a, S[ind_a])
  meanY_as <- meanS(matrix(Y[ind_a]), S[ind_a])
  
  meanr_bs <- meanS(r_b, S[ind_b])
  meanY_bs <- meanS(matrix(Y[ind_b]), S[ind_b])
  Xa <- matrix(X[ind_a,],ncol=p)
  Ya <- Y[ind_a]
  Xb <- matrix(X[ind_b,],ncol=p)
  Yb <- Y[ind_b]
  meanX_s <- meanS(X,S)
  for(s in S_ls){
    ind_as = which(S[ind_a]==s)
    ind_bs = which(S[ind_b]==s)
    pn_s = length(which(S==s))/length(S)
    n_as = length(which(S[ind_a]==s))
    n_bs = length(which(S[ind_b]==s))
    pi_as = n_as/length(which(S==s))
    pi_bs = n_bs/length(which(S==s))
    
    sd_r = sd_r + pn_s*(1/n_as*sum((r_a[ind_as] - meanr_as[s])^2)/pi_as + 
                          1/n_bs*sum((r_b[ind_bs]-meanr_bs[s])^2)/pi_bs)
    sd_hy = sd_hy + pn_s*((meanY_as[s]-sum(pS*meanY_as))-
                            (meanY_bs[s]-sum(pS*meanY_bs)))^2
    sd_rx = sd_rx + pn_s * (var(Xa[ind_as,]-mean(Xa[ind_as,]),Ya[ind_as]-mean(Ya[ind_as])) - 
                              var(Xb[ind_bs,]-mean(Xb[ind_bs,]),Yb[ind_bs]-mean(Yb[ind_bs])))
  }
  sd_rx = 2*t(beta_a-beta_b) %*% sd_rx
  sd_x = t(beta_a-beta_b)%*% var(X) %*% (beta_a-beta_b)
  
  sd = sd_r+sd_hy-sd_x + sd_rx 
  
  return(sd)
}

##---estimator for stratum specific tau ---##

tau_ss <- function(Y,X,S,A){
  nA <- length(unique(A))
  nS <- length(unique(S))
  p <- ncol(X)
  result <- rep(0,nA)
  pn_s <- rep(0,nS)
  for(s in sort(unique(S))){
    pn_s[s] <- length(which(S==s))/length(S)
  }
  beta_as <- coef_cal_ss(Y,S,A,X)
  meanXs <- meanS(X,S)
  for(i in 1:nA){
    ind <- which(A==(i-1))
    Ya <- Y[ind]
    Xa <- matrix(X[ind,],length(ind),p)
    Sa <- S[ind]
    meanYas <- meanS(matrix(Ya),Sa)
    meanXas <- meanS(Xa,Sa)
    result[i] <- sum(pn_s*(meanYas - diag((meanXas - meanXs)%*%beta_as[,((i-1)*nS+1):(i*nS)])))
  }
  
  return(result)
}

##---estimator for stratum specific sd---##
sd_ss <- function(Y,X,S,A,a,b){
  sd_r <- 0
  sd_x <- 0
  sd_hy <- 0
  N <- length(S)
  p <- ncol(X)
  S_ls <- sort(unique(S))
  nS <- length(S_ls)
  pS <- c(table(S)/length(S))
  
  ind_a <- which(A==a)
  ind_b <- which(A==b)
  
  beta <- coef_cal_ss(Y,S,A,X)
  
  beta_as <- matrix(beta[,(a*nS+1):((a+1)*nS)],p,nS)
  beta_bs <- matrix(beta[,(b*nS+1):((b+1)*nS)],p,nS)
  
  Ya <- Y[ind_a]
  Yb <- Y[ind_b]
  Xa <- matrix(X[ind_a,],length(ind_a),p)
  Xb <- matrix(X[ind_b,],length(ind_b),p)
  
  meanY_as <- meanS(matrix(Y[ind_a]), S[ind_a])
  meanY_bs <- meanS(matrix(Y[ind_b]), S[ind_b])
  
  
  for(s in S_ls){
    ind_s = which(S==s)
    ind_as = which(S[ind_a]==s)
    ind_bs = which(S[ind_b]==s)
    n_as = length(ind_as)
    n_bs = length(ind_bs)
    n_s = length(which(S==s))
    pi_as = n_as/n_s
    pi_bs = n_bs/n_s
    pn_s = n_s/N
    
    rss_sa = Ya[ind_as] - matrix(Xa[ind_as,],n_as,p) %*% beta_as[,s]
    rss_sb = Yb[ind_bs] - matrix(Xb[ind_bs,],n_bs,p) %*% beta_bs[,s]
    
    sd_r = sd_r + pn_s*(1/n_as*sum((rss_sa - mean(rss_sa))^2)/pi_as + 
                          1/n_bs*sum((rss_sb-mean(rss_sb))^2)/pi_bs)
    sd_hy = sd_hy + pn_s*((meanY_as[s]-sum(pS*meanY_as))-
                            (meanY_bs[s]-sum(pS*meanY_bs)))^2
    sd_x = sd_x + pn_s * t(beta_as[,s]-beta_bs[,s])%*% var(X[ind_s,]) %*% (beta_as[,s]-beta_bs[,s])
  }
  
  sd = sd_r+sd_hy+sd_x
  
  return(sd)
}

##---estimator of tau for fully saturated model---##
tau_sat <- function(Y,S,A){
  nA <- length(unique(A))
  nS <- length(unique(S))
  result <- rep(0,nA)
  pn_s <- rep(0,nS)
  for(s in sort(unique(S))){
    pn_s[s] <- length(which(S==s))/length(S)
  }
  for(i in 1:nA){
    ind <- which(A==(i-1))
    Ya <- Y[ind]
    Sa <- S[ind]
    meanYas <- meanS(matrix(Ya),Sa)
    result[i] <- sum(pn_s*(meanYas))
  }
  
  return(result)
  
}


##---estimator for fully saturated model---##
sd_sat<- function(Y,S,A,a,b){
  sd_y <- 0
  sd_hy <- 0
  
  ind_a <- which(A==a)
  ind_b <- which(A==b)
  
  Ya <- Y[ind_a]
  Yb <- Y[ind_b]
  
  pS <- c(table(S)/length(S))
  S_ls <- sort(unique(S))
  
  meanY_as <- meanS(matrix(Ya), S[ind_a])
  meanY_bs <- meanS(matrix(Yb), S[ind_b])
  
  for(s in S_ls){
    ind_as = which(S[ind_a]==s)
    ind_bs = which(S[ind_b]==s)
    pn_s = length(which(S==s))/length(S)
    n_as = length(which(S[ind_a]==s))
    n_bs = length(which(S[ind_b]==s))
    pi_as = n_as/length(which(S==s))
    pi_bs = n_bs/length(which(S==s))
    
    sd_y = sd_y + pn_s*(1/n_as*sum((Ya[ind_as] - meanY_as[s])^2)/pi_as + 
                          1/n_bs*sum((Yb[ind_bs]-meanY_bs[s])^2)/pi_bs)
    sd_hy = sd_hy + pn_s*((meanY_as[s]-sum(pS*meanY_as))-
                            (meanY_bs[s]-sum(pS*meanY_bs)))^2
  }
  
  sd = sd_y + sd_hy
  
  return(sd)
}