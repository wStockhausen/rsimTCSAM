#'
#'@title Calculate recruitment-at-size (in millions) by year, sex
#'
#'@param mc - model configuration object
#'@param showPlot - flag (T/F) to show plots
#'
#'@return list with elements:
#'devs_y: annual ln-scale recruitment deviations
#'R_y: annual total recruitment (millions)
#'Rx_c: (logit-scale) mean recruitment sex ratio (males), by time block
#'R_cz: proportion recruitment-at-size by time block
#'R_yx: proportion of recruitment by sex by year
#'R_yxz: 3d array with numbers of crab (in millions) recruiting by year/sex/size
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
    devs_y <- dimArray(mc,'y',val=NA); atts.devs_y<-attributes(devs_y);
    R_y    <- dimArray(mc,'y',val=NA); atts.R_y   <-attributes(R_y);
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        ndvs <- length(yrs);
        sdR <- sqrt(log(1+(tb$cvR)^2));
        devs <- rnorm(ndvs,mean=0, sd=sdR);
        devs<-devs-mean(devs);#enforce sum to zero
#        names(devs)<-yrs;
        devs_y[yrs] <- devs;
        r_y <- exp(tb$lnR+devs);
#        names(r_y) <- yrs;
        R_y[yrs] <- r_y;#changes R_y from array to vector for some reason
    }
#     devs_y<-as.array(devs_y,dim=dims,dimnames=dmnms);#change back to array
#     dimnames(devs_y)<-dmnms;#make sure names of dimnames are correct
#     R_y<-as.array(R_y,dim=dims,dimnames=dmnms);#change back to array
#     dimnames(R_y)<-dmnms;#make sure names of dimnames are correct
    attributes(devs_y)<-atts.devs_y; #reassign array attributes (which got stripped) to vector
    attributes(R_y)   <-atts.R_y;    #reassign array attributes (which got stripped) to vector
    
    #calc annual size distributions by time block
    R_cz  <- calcRatZ.All(mc,showPlot=showPlot)
    
    #calc sex-specific recruitment by year
    Rx_c  <- dimArray(mc,'pc_rec',val=0);   atts.Rx_c<-attributes(Rx_c);
    R_yx  <- dimArray(mc,'y.x',val=NA);
    R_yxz <- dimArray(mc,'y.x.z');
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs <-as.character(tb$years);
        ndvs<- length(yrs);
        if (tb$lgtSdXR>0){
            devs<- rnorm(ndvs,mean=0, sd=tb$lgtSdXR);
        } else {
            devs<-0*1:ndvs;
        }
        
        if (d$x$n==1){
            #1 sex, so ignore value
            Rx_c[t]          <- 1.0;
            R_yx[yrs,'male'] <- 1.0;
        } else {
            Rx_c[t] <- 1/(1+exp(-(tb$lgtMnXR+(tb$lgtSdXR^2)/2)));
            mXR     <- 1/(1+exp(-(tb$lgtMnXR+devs+(tb$lgtSdXR^2)/2)));
            R_yx[yrs,'male']   <- mXR;
            R_yx[yrs,'female'] <- 1-mXR;
        }
        for (y in yrs){
            for (x in d$x$nms){R_yxz[y,x,] <- R_y[y]*R_yx[y,x]*R_cz[t,];}
        }#y in yrs
    }#t
    attributes(Rx_c)<-atts.Rx_c;#re-assign array attributes to vector

    if (showPlot){
        mdfr<-melt(Rx_c,value.name='p')
        pc <- ggplot(mapping=aes(x=pc,y=p),data=mdfr)
        pc <- pc + geom_bar(stat='identity');
        pc <- pc + labs(x='time block',y='fraction male',title='sex ratio')
        pc <- pc + ylim(c(0,1.2))
        print(pc)
        mdfr<-melt(R_y,value.name='n')
        py <- ggplot(mapping=aes(x=y,y=n),data=mdfr)
        py <- py + geom_line();
        py <- py + labs(x='year',y='Total Annual Recruitment',title='Recruitment (millions)')
        print(py)
        mdfr<-melt(R_yx,value.name='p')
        px <- ggplot(mapping=aes(x=y,y=p),data=mdfr)
        px <- px + geom_line();
        px <- px + labs(x='year',y='fraction male',title='sex ratio')
        px <- px + ylim(c(0,1.2))
        print(px)
    }
    return(list(devs_y=devs_y,R_y=R_y,Rx_c=Rx_c,R_cz=R_cz,R_yx=R_yx,R_yxz=R_yxz));
}
#---------------------------------------------------------------------
#'
#'@title Calculate proportions recruiting-at-size over all time blocks
#'
#'@param mc - model configuration object
#'@param showPlot - flag (T/F) to show plots
#'
#'@return R_cz: d array with annual proportions (by time block) of crab recruiting by size
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
    R_cz  <- dimArray(mc,'pc_rec.z');
    mdfr<-NULL;
    for (t in names(p$blocks)){        
        tb<-p$blocks[[t]];
        prs<-calcRatZ(mc,tb,showPlot=FALSE);
	    R_cz[t,] <- prs;
    }
    
    if (showPlot){
        mdfr<-melt(R_cz,value.name='val')
        pz <- ggplot(mapping=aes(x=z,y=val),data=mdfr)
        pz <- pz + geom_bar(stat='identity');
        pz <- pz + labs(x='size (mm CW)',y='pr(Z)',title='Recruitment Size Distributions')
        pz <- pz + facet_grid(pc~.);
        print(pz)
    }
    return(R_cz)
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

