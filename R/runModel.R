#'
#'@title Run the simulation model
#'
#'@description Function to run the simulation model
#'
#'@param mc - model configuration object
#'@param mp - model processes object
#'@param showPlot - flag to show plots
#'
#'@return List consisting of:
#'iN_xmsz - initial population abundance by sex/maturity/shell condition/size
#'P_list - list of population time series (see calcNatZ)
#'F_list - list of fisheries results (see calcNatZ.Fisheries)
#'S_list - list of survey results (see calcNatZ.Surveys)
#'
#'@export
#'
runModel<-function(mc,mp,showPlot=TRUE){
    
    #calculate initial numbers-at-size
    iN_xmsz<-calcInitSizeComps(mc,mp,showPlot=showPlot);
    
    #calculate time series of population abundance and mature biomass
    P_list<-calcNatZ(mc,mp,iN_xmsz,showPlot=showPlot);
    N_yxmsz <- P_list$N_yxmsz;#numbers-at-size
    
    #calculate fishery catches
    F_list<-calcNatZ.Fisheries(mc,mp,N_yxmsz,showPlot=showPlot);
    
    #calculate survey catches
    S_list<-calcNatZ.Surveys(mc,mp,N_yxmsz,showPlot=showPlot);
    
    return(list(iN_xmsz=iN_xmsz,P_list=P_list,F_list=F_list,S_list=S_list))
}