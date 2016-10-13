#'
#'@title Get fishery catch biomass time series from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get fishery catch biomass time series from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM.rep object, rsimTCSAM.resLst object, or named list of the latter
#'@param category - 'captured','discarded','retained', 'discard mortality', 'total mortality', or 'index'
#'@param cast - casting formula for excluding y,x,m,s,z factor levels from an average-at-size across unspecified factors
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details Extracts the estimated fishery biomass time series.
#'
#'@export
#'
getMDFR.Fisheries.CatchBiomass<-function(rsims,
                                         category=c('captured','discarded','retained','discard mortality','total mortality','index'),
                                         cast="y+x",
                                         verbose=FALSE){
    if (verbose) cat("--rsimTCSAM::Getting fishery catch biomass time series.\n");

    category<-category[1];

    if (category=='captured'){
        path<-'mr/F_list/cpB_fyxmsz'; #total captured
    } else if (category=='discarded'){
        path<-'mr/F_list/dsB_fyxmsz'; #total discarded
    } else if (category=='retained'){
        path<-'mr/F_list/rmB_fyxmsz'; #total retained
    } else if (category=='discard mortality'){
        path<-'mr/F_list/dmB_fyxmsz'; #discard MORTALITY
    } else {
        cat("Category '",category,"' not recognized!\nReturning NULL...\n");
        return(NULL);
    }
    mdfr<-getMDFR(path,rsims,verbose);
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
