#'
#'@title Run the TCSAM simulation model to produce a TCSAM2015 input file.
#'
#'@description Function to run the TCSAM simulation model to produce a TCASM2015 input file.
#'
#'@param fn - filename for output file (input file to gmacs)
#'@param showPlot - flag (T/F) to show plots
#'
#'@return list with elements:
#'mc - model configuration list object
#'mp - model processes list object
#'mr - model results list object
#'
#'@export
#'
runSim.TCSAM<-function(fn='tcsam.input.dat',showPlot=TRUE){
    #get model configuration for Tanner crab
    mc <- readModelConfiguration.TCSAM();
#    mc <- ModelConfiguration.TCSAM();
    
    #calculate model processes
    mp <- calcModelProcesses(mc,showPlot=showPlot);
    
    #run the model
    mr <- runModel(mc,mp,showPlot=showPlot);
    
    #output results to model files
    writeSim.TCSAM(mc,mp,mr,fn=fn,showPlot=showPlot);
    
    return(invisible(list(mp=mp,mc=mc,mr=mr)));
}