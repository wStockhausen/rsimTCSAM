#'
#'@title Get population abundance time series from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get population abundance time series from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM.rep object, rsimTCSAM.resLst object, or named list of the latter
#'@param cast - casting formula for excluding y,x,m,s,z factor levels from an average-at-size across unspecified factors
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details Extracts the estimated population abundance time series.
#'
#'@export
#'
getMDFR.Pop.Abundance<-function(rsims,cast="y+x",verbose=FALSE){
    if (verbose) cat("--rsimTCSAM::Getting population abundance time series.\n");

    path<-'mr/P_list/N_yxmsz';
    mdfr<-getMDFR(path,rsims,verbose);
    mdfr$process<-'population';
    mdfr$type<-'predicted';
    mdfr<-removeImmOS(mdfr);

    castform<-"case+process+fleet+category+type+pc&&cast~.";
    castform<-gsub("&&cast",paste0("+",cast),castform,fixed=TRUE);
    ddfr<-reshape2::dcast(mdfr,castform,fun.aggregate=mean,na.rm=TRUE,value.var='val',drop=TRUE)
    ddfr[['.']]<-ifelse(ddfr[['.']]==0,NA,ddfr[['.']]);
    ddfr<-ddfr[!is.na(ddfr[['.']]),];#remove NA's

    mdfr<-rCompTCMs::getMDFR.CanonicalFormat(ddfr);

    if (verbose) cat("--Done. \n");
    return(mdfr);
}
