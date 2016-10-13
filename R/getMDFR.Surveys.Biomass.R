#'
#'@title Get survey biomass time series from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get survey biomass time series from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM.rep object, rsimTCSAM.resLst object, or named list of the latter
#'@param category - 'index' is only choice
#'@param cast - casting formula for excluding y,x,m,s,z factor levels from an average-at-size across unspecified factors
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details Extracts the estimated survey biomass time series.
#'
#'@export
#'
getMDFR.Surveys.Biomass<-function(rsims,category='index',cast="y+x",verbose=FALSE){
    if (verbose) cat("--rsimTCSAM::Getting survey biomass time series.\n");

    category<-category[1];

    path<-'mr/S_list/B_vyxmsz';
    mdfr<-getMDFR(path,rsims,verbose);
    mdfr$fleet<-gsub("_"," ",mdfr$fleet,fixed=TRUE);#replace '_'s in survey names with spaces
    mdfr$category<-category;
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
