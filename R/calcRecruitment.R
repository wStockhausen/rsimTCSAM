#'
#'@title Calculate recruitment-at-size by year, sex
#'
#'@param mc - model configuration object
#'@param showPlot - flag (T/F) to show plots
#'
#'@return R_yxz: 3d array with numbers  of crab recruiting by year/sex/size
#'
#'@import ggplot2
#'@import reshape2
#'
#'@export
#'
calcRecruitment<-function(mc,showPlot=TRUE){
    d <- mc$dims;      #model dimensions info
    p <- mc$params$rec;#recruitment parameters
    
    #calc total recruitment by year
    devs_y <- dimArray(mc,'y',val=NA);
    R_y    <- dimArray(mc,'y',val=NA);
    dims   <- dim(R_y);
    dmnms  <- dimnames(R_y);
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        ndvs <- length(yrs);
        sdR <- sqrt(log(1+(tb$cvR)^2));
        devs <- rnorm(ndvs,mean=0, sd=sdR);
        devs<-devs-mean(devs);#enforce sum to zero
#        names(devs)<-yrs;
        devs_y[yrs] <- devs;
        r_y <- exp(tb$lnR+devs-(sdR^2)/2);
#        names(r_y) <- yrs;
        R_y[yrs] <- r_y;#changes R_y from array to vector for some reason
    }
    devs_y<-as.array(devs_y,dim=dims,dimnames=dmnms);#change back to array
    dimnames(devs_y)<-dmnms;#make sure names of dimnames are correct
    R_y<-as.array(R_y,dim=dims,dimnames=dmnms);#change back to array
    dimnames(R_y)<-dmnms;#make sure names of dimnames are correct
    
    #calc sex-specific recruitment by year
    R_yx  <- dimArray(mc,'y.x',val=NA);
    if (d$x$n==1){
        sdXR <- 0;
        R_yx[1+(1:ndvs),1] <- 1;
    } else {
        for (t in names(p$blocks)){
            tb<-p$blocks[[t]];
            yrs <-as.character(tb$years);
            ndvs<- length(yrs);
            if (tb$lgtSdXR>0){
                devs<- rnorm(ndvs,mean=0, sd=tb$lgtSdXR);
            } else {
                devs<-0*1:ndvs;
            }
            mXR <- 1/(1+exp(-(tb$lgtMnXR+devs+(tb$lgtSdXR^2)/2)));
            R_yx[yrs,'male']   <- mXR;
            R_yx[yrs,'female'] <- 1-mXR;
        }#t
    }
    
    #calc annual size distributions for all years
    R_yz  <- calcRatZ.All(mc,showPlot=showPlot)
    
    #calc year/sex/size-specific recruimtent
    R_yxz <- dimArray(mc,'y.x.z');
    for (y in d$y$nms){
        #cat('y =',y,'\n')
        for (x in d$x$nms){
            #cat('x =',x,'\n')
            R_yxz[y,x,] <- R_y[y]*R_yx[y,x]*R_yz[y,];
        }
    }
    
    if (showPlot){
        mdfr<-melt(R_y,value.name='n')
        py <- ggplot(mapping=aes(x=y,y=n),data=mdfr)
        py <- py + geom_line();
        py <- py + labs(x='year',y='Total Annual Recruitment',title='Recruitment')
        print(py)
        mdfr<-melt(R_yx,value.name='p')
        px <- ggplot(mapping=aes(x=y,y=p),data=mdfr[mdfr$x=='male',])
        px <- px + geom_line();
        px <- px + labs(x='year',y='fraction male',title='sex ratio')
        px <- px + ylim(c(0,1))
        print(px)
    }
    return(list(devs_y=devs_y,R_y=R_y,R_yx=R_yx,R_yz=R_yz,R_yxz=R_yxz));
}
#---------------------------------------------------------------------
#'
#'@title Calculate proportions recruiting-at-size over all years
#'
#'@param mc - model configuration object
#'@param showPlot - flag (T/F) to show plots
#'
#'@return R_yz: d array with annual proportions of crab recruiting by size
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcRatZ.All<-function(mc,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'TC','calcRatZ()');
    }
    
    d <- mc$dims;      #model dimensions info
    p <- mc$params$rec;#recruitment parameters
    
    #calc size distribution
    R_yz  <- dimArray(mc,'y.z',val=NA);
    mdfr<-NULL;
    for (t in names(p$blocks)){        
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        prs<-calcRatZ(mc,tb,showPlot=FALSE);
	    for (y in yrs) {R_yz[y,] <- prs;}
        mdfrp<-melt(prs,value.name='val');
        mdfrp$t<-t;
        mdfr<-rbind(mdfr,mdfrp);
    }
    
    if (showPlot){
        pz <- ggplot(mapping=aes(x=z,y=val),data=mdfr)
        pz <- pz + geom_bar(stat='identity');
        pz <- pz + labs(x='size (mm)',y='pr(Z)',title='Recruitment Size Distributions')
        pz <- pz + facet_grid(t~.);
        print(pz)
    }
    return(R_yz)
}
#---------------------------------------------------------------------
#'
#'@title Calculate proportions recruiting-at-size for one set of parameters
#'
#'@param mc - model configuration list object
#'@param tb - time block or inits list object
#'@param showPlot - flag (T/F) to show plots
#'
#'@return 1d array with proportions of crab recruiting by size
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcRatZ<-function(mc,tb,showPlot=FALSE){
    
    d<-mc$dims;
    
    #calc size distribution
    alpha <- exp(tb$lnAlphaZ);
    beta  <- exp(tb$lnBetaZ);
    zbs   <- d$z$vls - d$zc$vls[1];#size increment from the lowest size cutpoint
    #print(zbs)
    prs<-dimArray(mc,'z');
    prs[] <- dgamma(zbs,shape=alpha/beta,scale=beta);
    prs   <- prs/sum(prs);#standardized to sum to 1
    
    if (showPlot){
        mdfr<-melt(prs,value.name='val');
        pz <- ggplot(mapping=aes(x=z,y=val),data=mdfr)
        pz <- pz + geom_bar(stat='identity');
        pz <- pz + labs(x='size (mm)',y='pr(Z)',title='Recruitment Size Distribution')
        print(pz)
    }
    return(prs)
}

