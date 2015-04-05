#'
#'@title Calculate model processes
#'
#'@description Function to calculate model processes prior to running the model
#'
#'@param mc - model configuration list object
#'@param showPlot - flag (T/F) to show plots
#'
#'@return mp - list object for model processes
#'
#'@export
#'
calcModelProcesses<-function(mc,showPlot=TRUE){
    #calculate weight-at-size
    W_yxmsz <- calcWatZ(mc,showPlot=showPlot);
    
    #calculate natural mortality
    M_yxmsz <- calcNaturalMortality(mc,showPlot=showPlot);
    
    #calculate molting probabilities for immature crab
    prMolt_yxsz<-calcPrMolt(mc,showPlot=showPlot);
    
    #calculate pr(molt to maturity|size, molt) 
    prMolt2Mat_yxsz <- calcPrMoltToMaturity(mc,showPlot=showPlot);
    
    #calculate size transition matrix for molting crab
    T_yxszz <- calcZTM(mc,showPlot=showPlot);
    
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
    
    mp <- list(W_yxmsz=W_yxmsz,
               M_yxmsz=M_yxmsz,
               prMolt_yxsz=prMolt_yxsz,
               prMolt2Mat_yxsz=prMolt2Mat_yxsz,
               T_yxszz=T_yxszz,
               sel_cz=sel_cz,
               S1_yxmsz=S1_yxmsz,
               S2_yxmsz=S2_yxmsz,
               R_list=R_list,
               F_list=F_list,
               S_list=S_list)
    
    return(mp)
}
