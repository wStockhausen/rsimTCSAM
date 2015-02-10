#'
#'@title Calculate the state transition matrices for a single sex.
#'
#'@description Function to calculate the state transition matrices for a single sex.
#'
#'@param mc     - model configuration list object
#'@param S1_msz - 3d array of pr(survival) from start of year to mating/molting by maturity state, shell condition, size class
#'@param P_sz   - 2d array of the probability by size class of molting for immature crab, by shell condition
#'@param Th_sz  - 2d array with pr(molt to maturity|size, molt) for immature crab by shell condition
#'@param T_szz  - 3d array with size transition matrix for growth by immature crab by shell condition
#'@param S2_msz - 3d array of pr(survival) from mating/molting to end of year by maturity state, shell condition, size class
#'
#'@return list with state transition matrices as named elements:
#'  A : imm, new -> imm, new
#'  B : imm, old -> imm, new
#'  C : imm, new -> imm, old
#'  D : imm, old -> imm, old
#'  E : imm, new -> mat, new
#'  F : imm, old -> mat, new
#'  G : mat, new -> mat, old
#'  H : mat, old -> mat, old
#'
#'@export
#'
calcStateTransitionMatrices<-function(mc,S1_msz,P_sz,Th_sz,T_szz,S2_msz){
    #create an identity matrix
    I <- diag(mc$dims$z$n);
    
    i<-1;#'immature'
    m<-2;#'mature';
    n<-1;#'new shell'
    o<-2;#'old shell'
    
    #immature new shell crab
    S2_in<-diag(S2_msz[i,n,]); #pr(survival|size) for immature new shell crab after molting occurs
    Tr_in<-T_szz[n,,];         #pr(Z_post|Z_pre, new shell, molt) for immature crab
    Th_in<-diag(Th_sz[n,]);    #pr(molt to maturity|pre-molt size,new shell, molting)
    Ph_in<-diag(P_sz[n,]);     #pr(molt|pre-molt size, new shell)
    S1_in<-diag(S1_msz[i,n,]); #pr(survival|size) for immature new shell crab before molting occurs
    #immature old shell crab
    S2_io<-diag(S2_msz[i,o,]); #pr(survival|size) for immature old shell crab after molting occurs (but they didn't molt)
    Tr_io<-T_szz[o,,];         #pr(Z_post|Z_pre, old shell, molt) for immature crab
    Th_io<-diag(Th_sz[o,]);    #pr(molt to maturity|pre-molt size,old shell, molting)
    Ph_io<-diag(P_sz[o,]);     #pr(molt|pre-molt size, old shell)
    S1_io<-diag(S1_msz[i,o,]); #pr(survival|size) for immature old shell crab before molting occurs
    #mature new shell crab
    S2_mn<-diag(S2_msz[i,n,]); #pr(survival|size) for mature new shell crab after molting occurs
    S1_mn<-diag(S1_msz[i,n,]); #pr(survival|size) for mature new shell crab before molting occurs (but they won't molt)
    #mature old shell crab
    S2_mo<-diag(S2_msz[i,o,]); #pr(survival|size) for mature old shell crab after molting occurs (but they didn't molt)
    S1_mo<-diag(S1_msz[i,o,]); #pr(survival|size) for mature old shell crab before molting occurs (but they won't molt)
    
    #full state transition matrices
    res<-list();
    res$A <- S2_in %*% t(Tr_in) %*% (I-Th_in) %*% Ph_in %*% S1_in;#imm, new -> imm, new
    res$B <- S2_in %*% t(Tr_io) %*% (I-Th_io) %*% Ph_io %*% S1_io;#imm, old -> imm, new
    res$C <- S2_io %*% (I-Ph_in) %*% S1_in;                    #imm, new -> imm, old
    res$D <- S2_io %*% (I-Ph_io) %*% S1_io;                    #imm, old -> imm, old
    res$E <- S2_mn %*% Tr_in %*% Th_in %*% Ph_in %*% S1_in;    #imm, new -> mat, new
    res$F <- S2_mn %*% Tr_io %*% Th_io %*% Ph_io %*% S1_io;    #imm, old -> mat, new
    res$G <- S2_mo %*% S1_mn;                                  #mat, new -> mat, old
    res$H <- S2_mo %*% S1_mo;                                  #mat, old -> mat, old
    
    return(res);
}
