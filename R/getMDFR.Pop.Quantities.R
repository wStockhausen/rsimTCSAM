#'
#'@title Extract population quantities from rsimTCSAM model runs as a dataframe
#'
#'@description Function to extracte population quantities from rsimTCSAM model runs.
#'
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param type - quantity to extract ("R_y","B_yxms","MB_yx","N_yxmsz","N_yxms","N_yxm","N_yx","iN_xmsz","fN_xmsz")
#'@param verbose - flag (T/F) to print debug info
#'
#'@return dataframe in canonical format
#'
#'@details none.
#'
#'@export
#'
getMDFR.Pop.Quantities<-function(rsims=NULL,
                                 type=c("R_y","MB_yx",
                                        "B_yxms","B_yxm","B_yx",
                                        "N_yxmsz","N_yxmz","N_yxz",
                                        "N_yxms","N_yxm","N_yx",
                                        "iN_xmsz","fN_xmsz"),
                                 verbose=FALSE){

    if (verbose) cat("rsimTCSAM::getMDFR.Pop.Quantities: Getting population trends\n");
    
    types<-c("R_y","MB_yx",
             "B_yxms","B_yxm","B_yx",
             "N_yxmsz","N_yxmz","N_yxz",
             "N_yxms","N_yxm","N_yx",
             "iN_xmsz","fN_xmsz");
    if (!(type[1] %in% types)){
        cat("rsimTCSAM::getMDFR.Pop.Quantities: Unknown type requested: '",type[1],"'.\n",sep='');
        return(NULL);
    }

    if (type[1]=="R_y"){
        #recruitment
        if (verbose) cat("Getting recruitment time series\n");
        path<-'mp/R_list/R_y';
        mdfr<-getMDFR(path,rsims,verbose=verbose);
        mdfr$m<-'immature';
        mdfr$s<-'new shell';
    }
    if (type[1]=="MB_yx"){
        #mature biomass at mating trends
        if (verbose) cat("Getting population mature biomass-at-mating trends\n");
        path<-'mr/P_list/MB_yx';
        mdfr<-getMDFR(path,rsims,verbose=verbose);
        mdfr$m<-'mature';
    }
    if (substr(type[1],1,3)=="B_y"){
        #biomass trends
        if (verbose) cat("Getting population biomass trends\n");
        path<-'mr/P_list/B_yxms';
        dfr<-getMDFR(path,rsims,verbose=verbose);
        if (type[1]=="B_yxms") mdfr<-dfr;
        if (type[1]=="B_yxm"){
            if (verbose) cat("Getting population B_yxm.\n");
            mdfr<-reshape2::dcast(dfr,case+y+x+m~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[5]<-'val';
        } 
        if (type[1]=="B_yx"){
            if (verbose) cat("Getting population B_yxm.\n");
            mdfr<-reshape2::dcast(dfr,case+y+x~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[4]<-'val';
        } 
    }
    if (substr(type[1],1,3)=="N_y"){
        #Population abundance-at-size
        if (verbose) cat("Getting population abundance-at-size\n");
        path<-'mr/P_list/N_yxmsz';
        dfr<-getMDFR(path,rsims,verbose=verbose);
        dfr<-removeImmOS(dfr);
        if (type[1]=="N_yxmsz") mdfr<-dfr;
        if (type[1]=="N_yxmz"){
            #abundance trends
            if (verbose) cat("Getting population N_yxmz.\n");
            mdfr<-reshape2::dcast(dfr,case+y+x+m+z~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[6]<-'val';
        }
        if (type[1]=="N_yxz"){
            #abundance trends
            if (verbose) cat("Getting population N_yxz.\n");
            mdfr<-reshape2::dcast(dfr,case+y+x+z~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[5]<-'val';
        }
        if (type[1]=="N_yxms"){
            #abundance trends
            if (verbose) cat("Getting population abundance trends N_yxms\n");
            mdfr<-reshape2::dcast(dfr,case+y+x+m+s~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[6]<-'val';
        }
        if (type[1]=="N_yxm"){
            #abundance trends
            if (verbose) cat("Getting population abundance trends N_yxm\n");
            mdfr<-reshape2::dcast(dfr,case+y+x+m~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[5]<-'val';
        }
        if (type[1]=="N_yx"){
            #abundance trends
            if (verbose) cat("Getting population abundance trends N_yx\n");
            mdfr<-reshape2::dcast(dfr,case+y+x~.,fun.aggregate=sum,value.var='val');
            names(mdfr)[4]<-'val';
        }
    }
    if (type[1]=="iN_xmsz"){
        #Initial population abundance-at-size
        if (verbose) cat("Getting initial population abundance-at-size\n");
        path<-'mr/iN_xmsz';
        mdfr<-getMDFR(path,rsims,verbose=verbose);
        mdfr<-removeImmOS(mdfr);
    }
    if (type[1]=="fN_xmsz"){
        #Final population abundance-at-size
        if (verbose) cat("Getting final population abundance-at-size\n");
        path<-'mr/P_list/N_yxmsz';
        dfr<-getMDFR(path,rsims,verbose=verbose);
        dfr<-removeImmOS(dfr);
        dfr$y<-as.numeric(dfr$y)
        dfrp<-reshape2::dcast(dfr,case~.,fun.aggregate=max,value.var='y');
        names(dfrp)[2]<-'y';
        mdfr<-NULL;
        for (case in dfrp$case){
            idx<-(dfr$case==case)&(dfr$y==dfrp$y[dfrp$case==case]);
            mdfrp<-dfr[idx,];
            mdfr<-rbind(mdfr,mdfrp);
        }
    }

    mdfr<-getMDFR.CanonicalFormat(mdfr);
    mdfr$type<-"population";

    if (verbose) cat("--rsimTCSAM::getMDFR.Pop.Quantities: Done. \n");
    return(mdfr);
}
