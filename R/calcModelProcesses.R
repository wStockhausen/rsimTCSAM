#'
#'@title Calculate model processes
#'
#'@description Function to calculate model processes prior to running the model.
#'
#'@param mc - model configuration list object
#'@param showPlot - flag (T/F) to show plots
#'
#'@details None.
#'
#'@return mp - list object for model processes
#'
#'@export
#'
calcModelProcesses<-function(mc,showPlot=TRUE){
    #calculate weight-at-size
    W_list  <- calcWatZ(mc,showPlot=showPlot);
    W_cxmz  <- W_list$W_cxmz;
    W_yxmsz <- W_list$W_yxmsz;
    
    #calculate natural mortality
    M_list <- calcNaturalMortality(mc,showPlot=showPlot);
    M_cxm   <- M_list$M_cxm;
    M_yxmsz <- M_list$M_yxmsz;
    
    #calculate molting probabilities for immature crab
    prMolt_yxsz<-calcPrMolt(mc,showPlot=showPlot);
    
    #calculate pr(molt to maturity|size, molt) 
    prM2M_list <- calcPrM2M(mc,showPlot=showPlot);
    prM2M_cxz  <- prM2M_list$prM2M_cxz;
    prM2M_yxsz <- prM2M_list$prM2M_yxsz;
    
    #calculate size transition matrix for molting crab
    T_list <- calcZTM(mc,showPlot=showPlot);
    
    #calculate selectivity/retention functions
    sel_cz <- calcSelFcns(mc,showPlot=showPlot);
    
    #calculate fishing mortalities
    F_list <- calcFishingMortalities(mc,sel_cz,showPlot=showPlot);
    
    #calculate survey catchabilities (TODO: implement this)
    S_list <- calcSurveyCatchabilities(mc,sel_cz,showPlot=showPlot);
    
    #calculate total mortality schedules
    #--mortality BEFORE mating (assumes fishing midpoint happens before mating)
    Z1_yxmsz <- mc$params$mate.time*M_yxmsz + F_list$tmF_yxmsz;;
    #--mortality AFTER mating (assumes fishing midpoint happens before mating)
    Z2_yxmsz <- (1-mc$params$mate.time)*M_yxmsz;
    
    #calculate survival
    S1_yxmsz <- exp(-Z1_yxmsz);#from start of year to mating, includes fishing
    S2_yxmsz <- exp(-Z2_yxmsz);#from mating to end of year
    
    #calculate recruitment
    R_list <- calcRecruitment(mc,showPlot=showPlot);#note this won't work for Tier 1
    
    mp <- list(W_cxmz=W_cxmz,
               W_yxmsz=W_yxmsz,
               M_cxm=M_cxm,
               M_yxmsz=M_yxmsz,
               prMolt_yxsz=prMolt_yxsz,
               prM2M_cxz=prM2M_cxz,
               prM2M_yxsz=prM2M_yxsz,
               T_list=T_list,
               sel_cz=sel_cz,
               S1_yxmsz=S1_yxmsz,
               S2_yxmsz=S2_yxmsz,
               R_list=R_list,
               F_list=F_list,
               S_list=S_list)
    
    return(mp)
}
