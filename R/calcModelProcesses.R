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
    prMoltToMat_yxsz <- calcPrMoltToMaturity(mc,showPlot=showPlot);
    
    #calculate size transition matrix for molting crab
    T_yxszz <- calcZTM(mc,showPlot=showPlot);
    
    #calculate fishing mortalities
    F_list <- calcFishingMortalities(mc,showPlot=showPlot);
    
    #calculate survey catchabilities (TODO: implement this)
    S_list <- calcSurveyCatchabilities(mc,showPlot=showPlot);
    
    #calculate time-varying total mortality
    Z_yxmsz <- M_yxmsz;
    for (f in mc$dims$fisheries$nms){
        Z_yxmsz[,,,,] <- Z_yxmsz[,,,,] + (F_list$F_fyxmsz)[f,,,,,];
    }
    
    #calculate survival
    S_yxmsz <- exp(-Z_yxmsz);
    
    mp <- list(W_yxmsz=W_yxmsz,
               M_yxmsz=M_yxmsz,
               Z_yxmsz=Z_yxmsz,
               S_yxmsz=S_yxmsz,
               prMolt_yxsz=prMolt_yxsz,
               prMoltToMat_yxsz=prMoltToMat_yxsz,
               T_yxszz=T_yxszz,
               F_list=F_list,
               S_list=S_list)
    
    return(mp)
}
