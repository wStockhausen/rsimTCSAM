#'
#'@title Run the TCSAM simulation model to produce a TCSAM2015 input file.
#'
#'@description Function to run the TCSAM simulation model to produce a TCASM2015 input file.
#'
#'@param fno - filename for output file (input file to TCSAM2015)
#'@param seed - seed for random number generator (NULL generates system seed)
#'@param fnd - filename for devs output
#'@param showPlot - flag (T/F) to show plots
#'@param pdf - if not NULL, a filename to output plots to in pdf format 
#'@param width - width of pdf output (inches)
#'@param height - height of pdf output (inches)
#'
#'@return list with elements:
#'seed - RNG seed used
#'mc - model configuration list object
#'mp - model processes list object
#'mr - model results list object
#'mo - model output list object
#'
#'@import reshape2
#'
#'@export
#'
runSim.TCSAM<-function(fno='tcsam.input.dat',
                       seed=NULL,
                       fnd='devs.csv',
                       showPlot=TRUE,
                       pdf=NULL,
                       width=8,
                       height=6){
    #set RNG seed
    set.seed(seed,kind='default',normal.kind='default');
    
    #crete pdf for pot output, if requested
    if (!is.null(pdf)){
        pdf(file=pdf,onefile=TRUE,width=width,height=height);
        on.exit(dev.off());
    }
    
    #get model configuration for Tanner crab
    mc <- readModelConfiguration.TCSAM();
#    mc <- ModelConfiguration.TCSAM();
    
    #calculate model processes
    mp <- calcModelProcesses(mc,showPlot=showPlot);
    
    #run the model
    mr <- runModel(mc,mp,showPlot=showPlot);
    
    #output results to model files
    mo <- writeSim.TCSAM(mc,mp,mr,fn=fno,showPlot=showPlot);

    #write out devs information to separate file
    mR <- melt(mp$R_list$devs_y,value.name='recdevs');
    mF <- melt(mp$F_list$devs_fy,value.name='fdevs');
    mQ <- melt(mp$S_list$devs_vy,value.name='qdevs');
    md <- cbind(mR,mF,mQ);
    write.csv(md,file=fnd,row.names=FALSE)
        
    return(invisible(list(seed=seed,mp=mp,mc=mc,mr=mr,mo=mo)));
}