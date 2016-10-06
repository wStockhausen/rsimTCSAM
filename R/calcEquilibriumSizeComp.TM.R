#'
#'@title Calculate the equilibrium size composition for one sex of a terminally-molting species with
#'two maturity states and two shell conditions
#'
#'@description Function to calculate the equilibrium size composition for one sex of a terminally-molting species with
#'two maturity states and two shell conditions.
#'
#'@param mc     - model configuration list object
#'@param R_z    - vector of equilibrium recruits-at-size
#'@param S1_msz - 3d array of pr(survival) from start of year to mating/molting by maturity state, shell condition, size class
#'@param P_sz   - 2d array of the probability by size class of molting for immature crab, by shell condition
#'@param Th_sz  - 2d array with pr(molt to maturity|size, molt) for immature crab by shell condition
#'@param T_szz  - 3d array with size transition matrix for growth by immature crab by shell condition
#'@param S2_msz - 3d array of pr(survival) from mating/molting to end of year by maturity state, shell condition, size class
#'
#'@return n_msz, a 3d array with equilibrium numbers-at-size
#'
#'@details None.
#'
#'@export
#'
calcEquilibriumSizeComp.TM<-function(mc,R_z,S1_msz,P_sz,Th_sz,T_szz,S2_msz){
    #create an identity matrix
    I <- diag(mc$dims$z$n);
    
    #calc the state transition matrices
    l<-calcStateTransitionMatrices(mc,S1_msz,P_sz,Th_sz,T_szz,S2_msz);
    
    #calculate inverses of matrix quantities
    iM1 <- solve(I - l$D);
    iM2 <- solve(I - l$A - l$B %*% iM1 %*% l$C);
    iM3 <- solve(I - l$H);
    
    #equilibrium solution
    imm.ns <- iM2 %*% R_z;                    #immature, new shell
    imm.os <- iM1 %*% l$C %*% imm.ns;         #immature, old shell
    mat.ns <- l$E %*% imm.ns + l$F %*% imm.os;#  mature, new shell
    mat.os <- iM3 %*% l$G %*% mat.ns;         #  mature, old shell
    
    n_msz<-dimArray(mc,'m.s.z');
    n_msz[1,1,] <- imm.ns;
    n_msz[1,2,] <- imm.os;
    n_msz[2,1,] <- mat.ns;
    n_msz[2,2,] <- mat.os;
    
    return(n_msz);
}
