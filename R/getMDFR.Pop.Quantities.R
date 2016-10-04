#'
#'@title Extract population quantities from rsimTCSAM model runs as a dataframe
#'
#'@description Function to extracte population quantities from rsimTCSAM model runs.
#'
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param type - quantity to extract ("R_y","B_yxms","MB_yx","N_yxmsz",or "N_yxms")
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details none.
#'
#'@export
#'
getMDFR.Pop.Quantities<-function(rsims=NULL,
                                 type=c("R_y","B_yxms","MB_yx","N_yxmsz","N_yxms"),
                                 verbose=FALSE){

    if (verbose) cat("rsimTCSAM::getMDFR.Pop.Quantities: Getting population biomass trends\n");

    if (type[1]=="R_y"){
        #recruitment
        if (verbose) cat("Getting recruitment time series\n");
        path<-'mp/R_list/R_y';
        mdfr<-getMDFR(path,rsims,verbose=verbose);
        mdfr$m<-'immature';
        mdfr$s<-'new shell';
    }
    if (type[1]=="B_yxms"){
        #biomass trends
        if (verbose) cat("Getting population biomass trends\n");
        path<-'mr/P_list/B_yxms';
        mdfr<-getMDFR(path,rsims,verbose=verbose);
    }
    if (type[1]=="MB_yx"){
        #mature biomass at mating trends
        if (verbose) cat("Getting population mature biomass-at-mating trends\n");
        path<-'mr/P_list/MB_yx';
        mdfr<-getMDFR(path,rsims,verbose=verbose);
        mdfr$m<-'mature';
    }
    if (substr(type[1],1,3)=="N_yxmsz"){
        #Population abundance-at-size
        if (verbose) cat("Getting population abundance-at-size\n");
        path<-'mr/P_list/N_yxmsz';
        dfr<-getMDFR(path,rsims,verbose=verbose);
        dfr<-removeImmOS(mdfr);
        if (type[1]=="N_yxmsz") mdfr<-dfr;
        if (type[1]=="N_yxms"){
            #abundance trends
            if (verbose) cat("Getting population abundance trends\n");
            mdfr<-reshape2::dcast(dfr,case+y+x+m+s~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[6]<-'val';
        }
    }

    mdfr<-getMDFR.CanonicalFormat(mdfr);
    mdfr$type<-"population";

    if (verbose) cat("--rsimTCSAM::getMDFR.Pop.Quantities: Done. \n");
    return(mdfr);
}
