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
        fn<-selectFile(ext)
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
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'WatZ');
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        a<-dimArray(list(dims=dims),'x.m');
        b<-dimArray(list(dims=dims),'x.m');
        for (xp in 1:dims$x$n){
            for (mp in 1:dims$m$n){
                x<-rsp[[i]][1]; 
                m<-rsp[[i]][2]; 
                av<-parseNum(rsp[[i]][3]); 
                bv<-parseNum(rsp[[i]][4]);
                a[x,m]<-av;
                b[x,m]<-bv;
                i<-i+1;
            }#mp
        }#xp
        blocks[[t]]<-list(years=years,
                          a_xm=a,
                          b_xm=b
                         );
    }#blocks
    params$wAtZ<-list(blocks=blocks);
    cat('--read weight-at-size parameters\n')
    
    #natural mortality
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'M');
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        M0_xms  <- dimArray(list(dims=dims),'x.m.s');
        cvM_xms <- dimArray(list(dims=dims),'x.m.s');
        for (xp in 1:dims$x$n){
            for (mp in 1:dims$m$n){
                for (sp in 1:dims$s$n){
                    x<-rsp[[i]][1]; 
                    m<-rsp[[i]][2]; 
                    s<-rsp[[i]][3]; 
                    M0_xms[x,m,s]  <- parseNum(rsp[[i]][4]);;
                    cvM_xms[x,m,s] <- parseNum(rsp[[i]][5]);;
                    i<-i+1;
                }
            }
        }
        blocks[[t]]<-list(years=years,
                          M0_xms=M0_xms,
                          cvM_xms=cvM_xms
                         );
    }#blocks
    params$nm<-list(blocks=blocks);
    cat('--read natural mortality parameters\n')
    
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
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'MoltToMaturity');
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
    params$moltToMaturity<-list(blocks=blocks);
    cat('--read molt-to-maturity parameters\n')
    
    #growth
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Growth');
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        a<-dimArray(list(dims=dims),'x');
        b<-dimArray(list(dims=dims),'x');
        s<-dimArray(list(dims=dims),'x');
        for (xp in 1:dims$x$n){
            x<-rsp[[i]][1]; 
            a[x]<-exp(parseNum(rsp[[i]][2])); 
            b[x]<-exp(parseNum(rsp[[i]][3])); 
            s[x]<-exp(parseNum(rsp[[i]][4])); 
            i<-i+1;
        }#xp
        blocks[[t]]<-list(years=years,
                          a_x=a,
                          b_x=b,
                          s_x=s
                         );
    }#blocks
    params$growth<-list(blocks=blocks);
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
    
    #fisheries
    resF<-parseMC.Fisheries(rsp,i,dims); i<-resF$i;
    params$fisheries<-resF$fisheries;
    cat('--read fisheries parameters\n')
    
    #surveys
    resS<-parseMC.Surveys(rsp,i,dims); i<-resS$i;
    params$surveys<-resS$surveys;
    cat('--read surveys parameters\n')
    
    #-----model configuration
    mc<-list(type=modelType,dims=dims,params=params)
    return(mc);
}

#mc<-readModelConfiguration();
