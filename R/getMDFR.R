#'
#'@title Get model arrays as a melted dataframe from rsimTCSAM models
#'
#'@description Function to get model objects as a melted dataframe from rsimTCSAM models.
#'
#'@param path - path in models to requested array (using '/' as separator for list levels)
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param verbose - flag (T/F) to print diagnostics
#'
#'@return Melted dataframe (ala package reshape2).
#'
#'@details Returned dataframe is a melted (ala reshape2) version of the requested array in
#'canonical format. The array value is in column 'val'. Uses \code{reshape2::melt(...)}.
#'
#'@export
#'
getMDFR<-function(path,rsims,verbose=FALSE){
    mdfr<-NULL;
    if (inherits(rsims,'rsimTCSAM')){
        #rsims is a rsimTCSAM model object
        obj<-getObj(path,rsims,verbose=verbose);
        if (!is.null(obj)){
            mdfrp<-reshape2::melt(obj,value.name='val',as.is=TRUE);
            mdfrp$case<-'rsim';
            mdfr<-rbind(mdfr,mdfrp)
        }
    } else if (class(rsims)=='list'){
        #rsims is a list of rsimTCSAM model objects
        nl<-length(rsims);
        nms<-names(rsims);
        for (l in 1:nl){
            rsim1<-rsims[[l]];
            mdfrp<-getMDFR(path,rsims=rsim1,verbose=verbose);
            if (!is.null(mdfrp)){
                if (!is.null(nms[l])) mdfrp$case<-nms[l];
                mdfr<-rbind(mdfr,mdfrp);
            }
        }
    }

    if (!is.null(mdfr)){
        cns<-colnames(mdfr);
        chks<-c('y','z','zp');
        for (chk in chks){
            idx<-which(cns==chk);
            if (length(idx)>0) mdfr[,chk]<-as.numeric(mdfr[,chk]);
        }
    }
    
    mdfr<-getMDFR.CanonicalFormat(mdfr);
    
    return(mdfr);
}
