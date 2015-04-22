#'
#'@title Select and read a model configuration file for a TCSAM model.
#'
#'@description Function to select and read a model configuration file for a TCSAM.
#'
#'@param fn - the file to read (NULL brings up a file chooser dialog)
#'@param ext - extension of type of file to select (if fn is NULL)
#'
#'@return model configuration list object
#'
#'@import wtsUtilities 
#'
#'@export
#'
readModelConfiguration.TCSAM<-function(fn=NULL,ext='*'){
    #get file name to read
    if (is.null(fn)){
        fn<-selectFile(ext);
        if (is.null(fn)) return(NULL);
    }
    
    conn<-file(fn,open='r');
    res<-readLines(con=conn);    
    close(conn);
    
    #parse res to remove comments and collect quoted and unquoted text elements
    rsp<-parseText(res,stripComments=TRUE);
    
    cat('Parsed file:\n')
    for (i in 1:length(rsp)){cat(rsp[[i]],'\n');}
    cat('\n\n')
    
    i<-1;
    chk<-rsp[[i]][1];i<-i+1;
    if (chk!='ModelConfiguration') {
        cat("First non-comment line should begin with keyword 'ModelConfiguration'.\n")
        cat('First word is "',chk,'".\n',sep='');
        cat('Aborting...\n');
        stop();
    }    
    
    #model type
    modelType<-rsp[[i]][1]; i<-i+1;
    cat("ModelType = '",modelType,"'\n",sep='');
    
    #parse model dimensions
    lst.dims<-parseMC.Dims(rsp,i); i<-lst.dims$i;
    dims<-lst.dims$dims;    
    
    mny<-dims$y$mny;#start year for simulation
    mxy<-dims$y$mxy;#final year for simulation
    asy<-dims$y$asy;#assessment year for simulation (mxy+1)
    zbs<-dims$z$vls;#size bins
    
    #--parse model parameters    
    params <- list();
    
    #PARAMETERS key word
    cat("reading parameters\n")
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'PARAMETERS');
    
    #mating and fishing times
    params$mate.time<-parseNum(rsp[[i]][1]); i<-i+1;
    params$fish.time<-parseNum(rsp[[i]][1]); i<-i+1;
    cat('--read mate.time, fish.time\n')
        
    #weight at size
    res<-parseMC.WatZ(rsp,i,dims); i<-res$i;
    params$wAtZ<-res$params;
    dims$wAtZ<-list();
    dims$wAtZ$n<-length(params$wAtZ$blocks);
    dims$wAtZ$nms<-names(params$wAtZ$blocks);
    dims$wAtZ$lbls<-names(params$wAtZ$blocks);
    cat("wAtZ = ",addQuotes(dims$wAtZ$lbls),'\n')
    cat('--read weight-at-size parameters\n')
    
    #natural mortality
    res<-parseMC.NaturalMortality(rsp,i,dims); i<-res$i;
    params$natmort<-res$params;
    dims$natmort<-list();
    dims$natmort$n<-length(params$natmort$blocks);
    dims$natmort$nms<-names(params$natmort$blocks);
    dims$natmort$lbls<-names(params$natmort$blocks);
    cat("natmort = ",addQuotes(dims$natmort$lbls),'\n')
    cat('--read natural mortality function parameters\n')
    
    #molting
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Molting');
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        z50_xs<-dimArray(list(dims=dims),'x.s');
        sdv_xs<-dimArray(list(dims=dims),'x.s');
        for (xp in 1:dims$x$n){
            for (sp in 1:dims$s$n){
                x<-rsp[[i]][1]; 
                s<-rsp[[i]][2]; 
                z50_xs[x,s]<-parseNum(rsp[[i]][3]); 
                sdv_xs[x,s]<-parseNum(rsp[[i]][4]); 
                i<-i+1;
            }#sp
        }#xp
        blocks[[t]]<-list(years=years,
                          z50_xs=z50_xs,
                          sdv_xs=sdv_xs
                          );
    }#blocks
    params$molting<-list(blocks=blocks);
    cat('--read molting parameters\n')
    
    #molt to maturity
    res<-parseMC.PrM2M(rsp,i,dims); i<-res$i;
    params$prM2M<-res$params;
    dims$prM2M<-list();
    dims$prM2M$n<-length(params$prM2M$blocks);
    dims$prM2M$nms<-names(params$prM2M$blocks);
    dims$prM2M$lbls<-names(params$prM2M$blocks);
    cat("prM2M = ",addQuotes(dims$prM2M$lbls),'\n')
    cat('--read molt-to-maturity parameters\n')
    
    #growth
    res<-parseMC.Growth(rsp,i,dims); i<-res$i;
    params$growth<-res$params;
    dims$growth<-list();
    dims$growth$n   <-length(params$growth$blocks);
    dims$growth$nms <-names(params$growth$blocks);
    dims$growth$lbls<-names(params$growth$blocks);
    cat("growth = ",addQuotes(dims$growth$lbls),'\n')
    cat('--read growth parameters\n')
    
    #recruitment
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Recruitment');
    #read initialization section
    inits<-list();
    inits$lnR      <- parseNum(rsp[[i]][1]); #ln-scale mean recruitment
    inits$cvR      <- parseNum(rsp[[i]][2]); #ln-scale value for ln-scale recruitment standard deviation
    inits$lgtMnXR  <- parseNum(rsp[[i]][3]); #logit-scale nominal sex ratio
    inits$lgtSdXR  <- parseNum(rsp[[i]][4]); #logit-scale standard deviation for sex ratio deviations
    inits$lnAlphaZ <- parseNum(rsp[[i]][5]); #ln-scale alpha parameter for rec. size distribution
    inits$lnBetaZ  <- parseNum(rsp[[i]][6]); #ln-scale beta parameter for rec. size distribution
    i<-i+1;
    #read time blocks
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        block<-list();
        t<-rsp[[i]][1];
        eval(parse(text=paste('years<-',t)));
        block$years<-years;
        block$lnR      <- parseNum(rsp[[i]][2]); #ln-scale mean recruitment
        block$cvR      <- parseNum(rsp[[i]][3]); #ln-scale value for ln-scale recruitment standard deviation
        block$lgtMnXR  <- parseNum(rsp[[i]][4]); #logit-scale nominal sex ratio
        block$lgtSdXR  <- parseNum(rsp[[i]][5]); #logit-scale standard deviation for sex ratio deviations
        block$lnAlphaZ <- parseNum(rsp[[i]][6]); #ln-scale alpha parameter for rec. size distribution
        block$lnBetaZ  <- parseNum(rsp[[i]][7]); #ln-scale beta parameter for rec. size distribution
        blocks[[t]]<-block;
        i<-i+1;
    }#blocks
    params$rec<-list(inits=inits,
                     blocks=blocks);
    cat('--read recruitment parameters\n')
    
    #selectivity/retention functions
    res<-parseMC.SelFcns(rsp,i,dims); i<-res$i;
    params$selfcns<-res$selfcns;
    dims$selfcns$n<-length(params$selfcns);
    dims$selfcns$nms<-1:dims$selfcns$n;
    dims$selfcns$lbls<-vector(mode='character',length=dims$selfcns$n);
    for (s in 1:dims$selfcns$n){
        dims$selfcns$lbls[s]<-params$selfcns[[s]]$label;
    }
    cat("selfcns = ",addQuotes(dims$selfcns$lbls),'\n')
    cat('--read selectivity/retention function parameters\n')
    
    #fisheries
    res<-parseMC.Fisheries(rsp,i,dims); i<-res$i;
    params$fisheries<-res$fisheries;
    cat('--read fisheries parameters\n')
    
    #surveys
    res<-parseMC.Surveys(rsp,i,dims); i<-res$i;
    params$surveys<-res$surveys;
    cat('--read surveys parameters\n')
    
    #-----model configuration
    mc<-list(type=modelType,dims=dims,params=params)
    return(mc);
}

#mc<-readModelConfiguration.TCSAM();
