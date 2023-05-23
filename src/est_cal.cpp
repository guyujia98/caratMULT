#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// write functions to calculate the estimators. 
// X only contains additional covariates
// coefficients for stratum common regression
// [[Rcpp::export]]
arma::vec coef_cal_sc(arma::vec Y, arma::vec S, arma::mat X){
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  vec poss_S = unique(S);
  int n_s = poss_S.n_elem;
  
  vec beta = vec(ncol);
  
  mat meanX = mat(nrow,ncol,fill::zeros);
  vec meanY = vec(nrow,fill::zeros);
  mat SigmaX = mat(ncol,ncol,fill::zeros);
  vec SigmaXY = vec(ncol,fill::zeros);

  for(int i=0; i<n_s;i++){
    double s = poss_S(i);
    uvec ind = find(S == s);
    meanX.rows(ind) = repmat(mean(X.rows(ind)),ind.n_elem,1);
    meanY(ind) = ones<vec>(ind.n_elem) * mean(Y(ind));
  }

  for(int j=0;j<nrow;j++){
    SigmaX += (X.row(j) - meanX.row(j)).t() * (X.row(j)-meanX.row(j));
    SigmaXY += (X.row(j) - meanX.row(j)).t() * (Y(j) - meanY(j));
  }

  beta = solve(SigmaX, SigmaXY);

  return(beta);
 // return(meanY);
}

//coefficients for stratum specific regression
// [[Rcpp::export]]
arma::mat coef_cal_ss(arma::vec Y, arma::vec S, arma::vec A, arma::mat X){
  int np = X.n_cols;
  vec poss_A = sort(unique(A));
  vec poss_S = sort(unique(S));
  int n_a = poss_A.n_elem;
  int n_s = poss_S.n_elem;
  
  mat beta = mat(np,n_a*n_s,fill::zeros);
  
  
for(int k=0;k<n_a;k++){
//  int k=0;
  double a = poss_A(k);
  uvec ind_a = find(A==a);
  
  mat Xa = X.rows(ind_a);
  vec Ya = Y(ind_a);
  vec Sa = S(ind_a);

  for(int i=0; i<n_s;i++){
//    int i=0;
    double s = poss_S(i);
    uvec ind = find(Sa == s);
    
    mat Xas = Xa.rows(ind);
    vec Yas = Ya(ind);
    
    int nrow = Xas.n_rows;
    
    mat meanX = mat(nrow,np,fill::zeros);
    vec meanY = vec(nrow,fill::zeros);
    
    mat SigmasX = mat(np,np,fill::zeros);
    vec SigmasXY = vec(np,fill::zeros);
    
    meanX = repmat(mean(Xas),nrow,1);
    meanY = ones<vec>(nrow) * mean(Yas);
  
//    int j=0;
  for(int j=0;j<nrow;j++){
    SigmasX += (Xas.row(j) - meanX.row(j)).t() * (Xas.row(j)-meanX.row(j));
    SigmasXY += (Xas.row(j) - meanX.row(j)).t() * (Yas(j) - meanY(j));
  }
  
  beta.col(i+k*n_s) = solve(SigmasX,SigmasXY);
}
}  
  return(beta);
  // return(meanY);
}

//coefficients for stratum specific regression
// [[Rcpp::export]]
arma::mat coef_cal_sc_new(arma::vec Y, arma::vec S,arma::mat X, arma::vec pS, arma::vec pi_as){
  int np = X.n_cols;
  vec poss_S = sort(unique(S));
  int n_s = poss_S.n_elem;
  
  vec beta = vec(np);
  mat SigmaX = mat(np,np,fill::zeros);
  vec SigmaXY = vec(np,fill::zeros);
  
    for(int i=0; i<n_s;i++){
      //    int i=0;
      double s = poss_S(i);
      uvec ind = find(S == s);
      
      mat Xas = X.rows(ind);
      vec Yas = Y(ind);
      
      int nrow = Xas.n_rows;
      
      mat meanX = mat(nrow,np,fill::zeros);
      vec meanY = vec(nrow,fill::zeros);
      
      mat SigmasX = mat(np,np,fill::zeros);
      vec SigmasXY = vec(np,fill::zeros);
      
      meanX = repmat(mean(Xas),nrow,1);
      meanY = ones<vec>(nrow) * mean(Yas);
      
      //    int j=0;
      for(int j=0;j<nrow;j++){
        SigmasX += (Xas.row(j) - meanX.row(j)).t() * (Xas.row(j)-meanX.row(j))/nrow;
        SigmasXY += (Xas.row(j) - meanX.row(j)).t() * (Yas(j) - meanY(j))/nrow;
      }
      SigmaX += pS(i) * pi_as(i) * SigmasX;
      SigmaXY += pS(i) * pi_as(i) * SigmasXY;
    }
    beta = solve(SigmaX,SigmaXY);
  return(beta);
  // return(meanY);
}
