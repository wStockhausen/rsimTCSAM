#'
#'@title Get growth transition matrices from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get growth transition matrices from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canconical format.
#'Note that 'z' is pre-molt size, 'zp' is post-molt size.
#'
#'@details Extracts growth transition matrices.
#'
#'@export
#'
getMDFR.Pop.GrowthTransitionMatrices<-function(rsims,verbose=FALSE){
    if (verbose) cat("--Getting growth transition matrices.\n");
    if (inherits(rsims,'rsimTCSAM')){
        rsims<-list(rsim=rsims);#wrap in list
    }
    
    mdfrp<-getMDFR('mp/T_list/T_cxzz',rsims,verbose);
    mdfrp$y<-'';
    ums<-as.character(unique(mdfrp$case))
    for (um in ums){
        idx<-(mdfrp$case==um);
        mdfrp$y[idx]<-reformatTimeBlocks(mdfrp$pc[idx],rsims[[um]]$mc$dims);
        mdfrp<-mdfrp[,c('case','pc','y','x','z','zp','val')];
    }
    #in mdfrp above, 'z' is post-molt size, 'zp' is pre-molt size
    #"transpose" matrices so 'z' represents pre-molt size, 'zp' post-molt size
    zp<-mdfrp$z;
    mdfrp$z<-mdfrp$zp;
    mdfrp$zp<-zp;
    
    mdfr<-getMDFR.CanonicalFormat(mdfrp);
    
    if (verbose) cat("--Done. \n");
    return(mdfr);
}
