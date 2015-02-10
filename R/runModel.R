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
#'N_yxmsz - population abundance by year/sex/maturity/shell condition/size
#'F_list - list of fisheries results (see calcNatZ.Fisheries)
#'N_vyxmsz - survey catches (numbers) by survey/year/sex/maturity/shell condition/size
#'
#'@export
#'
runModel<-function(mc,mp,showPlot=TRUE){
    
    #calculate initial numbers-at-size
    iN_xmsz<-calcInitSizeComps(mc,mp,showPlot=showPlot);
    
    #calculate time series of population abundance and mature biomass
    pop.ts<-calcNatZ(mc,mp,iN_xmsz,showPlot=showPlot);
    N_yxmsz <- pop.ts$N_yxmsz;#numbers-at-size
    MB_yx   <- pop.ts$MB_yx;  #mature biomass
    
    #calculate fishery catches
    F_list<-calcNatZ.Fisheries(mc,mp,N_yxmsz,showPlot=showPlot);
    
    #calculate survey catches
    N_vyxmsz<-calcNatZ.Surveys(mc,mp,N_yxmsz,showPlot=showPlot);
    
    return(list(iN_xmsz=iN_xmsz,N_yxmsz=N_yxmsz,MB_yx=MB_yx,F_list=F_list,N_vyxmsz=N_vyxmsz))
}