#'
#'@title Get fishery selectivity functions from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get fishery selectivity functions from model results from TCSAM02 model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM object, rsimTCSAM.resLst object, or named list of the latter
#'@param verbose - flag (T/F) to print debug info
#'@param cast - casting formula for excluding x,m,s factor levels from an average-at-size across unspecified factors
#'
#'@return dataframe in canonical format
#'
#'@details Extracts fishery selectivity functions.
#'
#'@export
#'
getMDFR.Fisheries.SelFcns<-function(rsims,verbose=FALSE,cast="x"){
    if (verbose) cat("--rTCSAM02::Getting fishery selectivity functions.\n");

    path<-'mp/F_list/sel_fyxmsz';
    mdfr<-getMDFR(path,rsims,verbose);
    mdfr$fleet<-gsub("_"," ",mdfr$fleet,fixed=TRUE);#replace '_'s in fishery names with spaces
    mdfr$z<-as.numeric(mdfr$z);

    castform<-"case+type+fleet+category+pc+y&&cast+z~.";
    castform<-gsub("&&cast",paste0("+",cast),castform,fixed=TRUE);
    ddfr<-reshape2::dcast(mdfr,castform,fun.aggregate=mean,na.rm=TRUE,value.var='val',drop=TRUE)
    ddfr[['.']]<-ifelse(ddfr[['.']]==0,NA,ddfr[['.']]);
    ddfr<-ddfr[!is.na(ddfr[['.']]),];#remove NA's

    mdfr<-getMDFR.CanonicalFormat(ddfr);

    if (verbose) cat("--Done. \n");
    return(mdfr);
}
