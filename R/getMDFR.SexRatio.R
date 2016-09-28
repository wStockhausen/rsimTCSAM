#'
#'@title Get recruitment sex ratio from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get recruitment sex ratio from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format.
#'
#'@details Extracts recruitment sex ratio.
#'
#'@export
#'
getMDFR.SexRatio<-function(rsims,verbose=FALSE){
    if (verbose) cat("--Getting recruitment sex ratio.\n");
    if (inherits(rsims,'rsimTCSAM')){
        rsims<-list(rsim=rsims);#wrap in list
    }
    
    path<-'mp/R_list/Rx_c';
    mdfrp<-getMDFR(path,rsims,verbose);
    mdfrp$y<-'';
    ums<-as.character(unique(mdfrp$case))
    for (um in ums){
        idx<-(mdfrp$case==um);
        mdfrp$y[idx]<-reformatTimeBlocks(mdfrp$pc[idx],rsims[[um]]$mc$dims);
    }
    mdfrp<-mdfrp[,c('case','pc','y','val')];

    mdfr<-getMDFR.CanonicalFormat(mdfrp);

    if (verbose) cat("--Done. \n");
    return(mdfr);
}
