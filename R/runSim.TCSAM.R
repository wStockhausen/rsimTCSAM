#'
#'@title Run the TCSAM simulation model to produce a TCSAM2015 input file.
#'
#'@description Function to run the TCSAM simulation model to produce a TCASM2015 input file.
#'
#'@param out.dir - directory path for output files (input files to TCSAM2015)
#'@param seed - seed for random number generator (NULL generates system seed)
#'@param fnd - filename for devs output
#'@param showPlot - flag (T/F) to show plots
#'@param pdf - if not NULL, a filename to output plots to in pdf format 
#'@param width - width of pdf output (inches)
#'@param height - height of pdf output (inches)
#'
#'@return list with elements:
#'seed - RNG seed used
#'mc - model configuration list object (see readModelConfiguration.TCSAM(...))
#'mp - model processes list object     (see calcModelProcesses(...))
#'mr - model results list object       (see runModel(...))
#'md - model data                      (see writeSim.TCSAM(...))
#'mo - model OFL list object           (see calcOFL(...))
#'
#'@import reshape2
#'
#'@export
#'
runSim.TCSAM<-function(out.dir='.',
                       seed=NULL,
                       fnd='devs.csv',
                       showPlot=TRUE,
                       pdf=NULL,
                       width=8,
                       height=6){
    #set RNG seed
    set.seed(seed,kind='default',normal.kind='default');
    
    #get model configuration for Tanner crab
    mc <- readModelConfiguration.TCSAM();
    if (is.null(mc)) {
        cat('Model configuration file not read\n');
        cat('Returning NULL\n');
        return(NULL);
    }
    
    #create pdf for plot output, if requested
    if (!is.null(pdf)){
        pdf(file=pdf,onefile=TRUE,width=width,height=height);
        on.exit(dev.off());
    }
    
    #calculate model processes
    mp <- calcModelProcesses(mc,showPlot=showPlot);
    
    #run the model
    mr <- runModel(mc,mp,showPlot=showPlot);
    
    #output results to model files
    md <- writeSim.TCSAM(mc,mp,mr,out.dir=out.dir,showPlot=showPlot);

    #compare initial, final size comps
    sizecomps<-list();
    sizecomps[['initial']]<-mr$P_list$N_yxmsz[as.character(mc$dims$y$mny),,,,];
    sizecomps[['final']]  <-mr$P_list$N_yxmsz[as.character(mc$dims$y$asy),,,,];
    compareSizeCompsGG(n_xmsz=sizecomps,title='Size Compositions')

    #write out devs information to csv file
    mR <- melt(mp$R_list$devs_y,value.name='recdevs');
    mF <- melt(mp$F_list$devs_fy,value.name='fdevs');
    mQ <- melt(mp$S_list$devs_vy,value.name='qdevs');
    mdvs <- cbind(mR,mF,mQ);
    write.csv(mdvs,file=file.path(out.dir,fnd),row.names=FALSE);
    
    #calc OFL
    ##TODO->implement this: 
    ##mo<-calcOFL(mc,mp,mr,showPlot=showPlot);
    mo<-NULL;

    #create ouput list
    rsim<-list(seed=seed,mc=mc,mp=mp,mr=mr,md=md,mo=mo);
    class(rsim)<-'rsimTCSAM';
        
    return(invisible(rsim));
}