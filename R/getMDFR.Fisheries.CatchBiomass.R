#'
#'@title Get fishery catch biomass time series from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get fishery catch biomass time series from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM object, rsimTCSAM.resLst object, or named list of the latter
#'@param category - 'captured','discarded','retained', 'discard mortality', or 'index'
#'@param cast - casting formula for excluding x,m,s factor levels from an average-at-size across unspecified factors
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details Extracts the estimated fishery biomass time series.
#'
#'@export
#'
getMDFR.Fisheries.CatchBiomass<-function(rsims,
                                         category=c('captured','discarded','retained','discard mortality','index'),
                                         cast="x",
                                         verbose=FALSE){
    if (verbose) cat("--rsimTCSAM::Getting fishery catch biomass time series.\n");

    category<-category[1];

    if (category=='captured'){
        path<-'mr/F_list/cpB_fyxms'; #total captured
    } else if (category=='discarded'){
        path<-'mr/F_list/dsB_fyxms'; #total discarded
    } else if (category=='retained'){
        path<-'mr/F_list/rmB_fyxms'; #total retained
    } else if (category=='discard mortality'){
        path<-'mr/F_list/dmB_fyxms'; #discard MORTALITY
    } else {
        cat("Category '",category,"' not recognized!\nReturning NULL...\n");
        return(NULL);
    }
    mdfr<-getMDFR(path,rsims,verbose);
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
