#'
#'@title Get natural mortality rates from model results from rsimTCSAM model runs as a dataframe
#'
#'@description Function to get natural mortality rates from model results from rsimTCSAM model runs as a dataframe.
#'
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param type - flag indicating which M's to extract ('M_cxm', 'M_yxm' or 'M_yxmsz') 
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details Extracts natural mortality rates
#'
#'@export
#'
getMDFR.Pop.NaturalMortality<-function(rsims,
                                      type=c('M_cxm','M_yxm','M_yxmsz'),
                                      verbose=FALSE){
    if (verbose) cat("--Getting natural mortality info\n");
    if (inherits(rsims,'rsimTCSAM')){
        #rsims is a rsimTCSAM model object
        rsims<-list(rep=rsims);
        class(rsims)<-c("rsimTCSAM.resLst",class(rsims));
    }
    if (inherits(rsims,'rsimTCSAM.resLst')){
        #rsims is a rsimTCSAM.resLst model object
        rsims<-list(rsim=rsims);
    }
    
    if (type[1]=='M_cxm'){
        mdfrp<-getMDFR('mp/M_cxm',rsims,verbose=verbose);
        mdfrp$y<-'';
        ums<-as.character(unique(mdfrp$case))
        for (um in ums){
            idx<-(mdfrp$case==um);
            rsim<-rsims[[um]];
            if (inherits(rsim,"rsimTCSAM.resLst")) rsim<-rsim$rep;
            mdfrp$y[idx]<-reformatTimeBlocks(mdfrp$pc[idx],rsim$mc$dims)
        }
        mdfrp<-mdfrp[,c("case","pc","y","x","m","val")];
    } else if (type[1]=='M_yxm'){
        mdfrp<-getMDFR('mp/M_yxmsz',rsims,verbose=verbose);
        mdfrp<-reshape2::dcast(mdfrp,formula='case+y+x+m~.',fun.aggregate=mean,value.var='val');
        names(mdfrp)<-c("case","y","x","m","val");
    } else if (type[1]=='M_yxmsz'){
        mdfrp<-getMDFR('mp/M_yxmsz',rsims,verbose=verbose);
        mdfrp<-mdfrp[,c("case","y","x","m","s","z","val")];
    }

    mdfr<-getMDFR.CanonicalFormat(mdfrp);

    if (verbose) cat("--Done. \n");
    return(mdfr);
}
