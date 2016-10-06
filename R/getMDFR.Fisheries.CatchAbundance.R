#'
#'@title Get fishery catch abundance time series from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get fishery catch abundance time series from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM object, rsimTCSAM.resLst object, or named list of the latter
#'@param category - 'captured','discarded','retained', 'discard mortality', or 'index'
#'@param cast - casting formula for excluding x,m,s factor levels from an average-at-size across unspecified factors
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details Extracts the estimated fishery abundance time series.
#'
#'@export
#'
getMDFR.Fisheries.CatchAbundance<-function(rsims,
                                           category=c('captured','discarded','retained','discard mortality','index'),
                                           cast="x",
                                           verbose=FALSE){
    if (verbose) cat("--rsimTCSAM::Getting fishery catch abundance time series.\n");

    category<-category[1];

    if (category=='captured'){
        path<-'mr/F_list/cpN_fyxmsz'; #total captured
    } else if (category=='discarded'){
        path<-'mr/F_list/dsN_fyxmsz'; #total discarded
    } else if (category=='retained'){
        path<-'mr/F_list/rmN_fyxmsz'; #total retained
    } else if (category=='discard mortality'){
        path<-'mr/F_list/dmN_fyxmsz'; #discard MORTALITY
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
