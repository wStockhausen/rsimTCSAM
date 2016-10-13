#'
#'@title Get recruitment size distribution from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get recruitment size distribution from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format.
#'
#'@details Extracts mean growth increments.
#'
#'@export
#'
getMDFR.Pop.RecSizeDistribution<-function(rsims,verbose=FALSE){
    if (verbose) cat("--Getting recruitment size distribution.\n");

    if (inherits(rsims,'rsimTCSAM')){
        #rsims is a rsimTCSAM model object
        rsims<-list(rep=rsims);
        class(rsims)<-c("rsimTCSAM.resLst",class(rsims));
    }
    if (inherits(rsims,'rsimTCSAM.resLst')){
        #rsims is a rsimTCSAM.resLst model object
        rsims<-list(rsim=rsims);
    }
    
    path<-'mp/R_list/R_cz';
    mdfrp<-getMDFR(path,rsims,verbose);
    mdfrp$y<-'';
    ums<-as.character(unique(mdfrp$case))
    for (um in ums){
        idx<-(mdfrp$case==um);
        rsim<-rsims[[um]];
        if (inherits(rsim,"rsimTCSAM.resLst")) rsim<-rsim$rep;
        mdfrp$y[idx]<-reformatTimeBlocks(mdfrp$pc[idx],rsim$mc$dims);
    }
    mdfrp<-mdfrp[,c('case','pc','y','z','val')];

    mdfr<-rCompTCMs::getMDFR.CanonicalFormat(mdfrp);
    mdfr$process<-"population";
    mdfr$x<-"all";
    mdfr$m<-"immature";
    mdfr$s<-"all";


    if (verbose) cat("--Done. \n");
    return(mdfr);
}
