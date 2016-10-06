#'
#'@title Get survey abundance time series from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get survey abundance time series from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM.rep object, rsimTCSAM.resLst object, or named list of the latter
#'@param category - 'index' is only choice
#'@param cast - casting formula for excluding x,m,s factor levels from an average-at-size across unspecified factors
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details Extracts the estimated survey abundance time series.
#'
#'@export
#'
getMDFR.Surveys.Abundance<-function(rsims,category='index',cast="x",verbose=FALSE){
    if (verbose) cat("--rsimTCSAM::Getting survey abundance time series.\n");

    category<-category[1];

    path<-'mr/S_list/N_vyxmsz';
    mdfr<-getMDFR(path,rsims,verbose);
    mdfr$fleet<-gsub("_"," ",mdfr$fleet,fixed=TRUE);#replace '_'s in survey names with spaces
    mdfr$category<-category;
    mdfr<-removeImmOS(mdfr);

    castform<-"case+type+fleet+category+pc+y&&cast~.";
    castform<-gsub("&&cast",paste0("+",cast),castform,fixed=TRUE);
    ddfr<-reshape2::dcast(mdfr,castform,fun.aggregate=mean,na.rm=TRUE,value.var='val',drop=TRUE)
    ddfr[['.']]<-ifelse(ddfr[['.']]==0,NA,ddfr[['.']]);
    ddfr<-ddfr[!is.na(ddfr[['.']]),];#remove NA's

    mdfr<-getMDFR.CanonicalFormat(ddfr);

    if (verbose) cat("--Done. \n");
    return(mdfr);
}
