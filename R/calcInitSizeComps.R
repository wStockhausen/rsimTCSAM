#'
#'@title Calculate initial size compositions.
#'
#'@description Function to calculate initial size compositions.
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'
#'@return n_xmsz: 4-d array of initial population abundance by sex/maturity/shell condition/size
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcInitSizeComps<-function(mc,mp,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'TC','calcInitSizeDist');
    }
    
    d<-mc$dims;
    p<-mc$params$rec;
    
    #get initial recruitment information
    R0   <- NA;              #total equilibrium recruitment
    r_x  <- dimArray(mc,'x');#equilibrium fraction of recruitment by sex
    r_z  <- dimArray(mc,'z');#equilibrium fraction of recruitment by size
    t<-names(p$blocks)[1];   #use equilibrium values from FIRST recruitment time block
    tb<-p$blocks[[t]];
    if (!(d$y$mny %in% tb$years)){
        cat('First recruitment time block should contain model start year.\n');
        cat('Aborting...\n');
        stop();
    }
    sdR <- sqrt(log(1+(tb$cvR)^2));
    R0  <- exp(tb$lnR-(sdR^2)/2);
    if (d$x$n==1){
        r_x['male'] <- 1;
    } else {
        mXR <- 1/(1+exp(-tb$lgtMnXR));
        r_x['male']   <- mXR;
        r_x['female'] <- 1-mXR;
    }
    r_z[]<-mp$R_list$R_yz[1,];#assumed same as 1st year of model
    
    #calculate equilibrium size distributions
    n_xmsz<-dimArray(mc,'x.m.s.z');
    for (x in d$x$nms){
        R_z    <- R0*r_x[x]*r_z;             #recruits-at-size
        S1_msz <- mp$S1_yxmsz[1,x,,,];       #survival to mating/molting
        P_sz   <- mp$prMolt_yxsz[1,x,,];     #pr(molt)
        Th_sz  <- mp$prMolt2Mat_yxsz[1,x,,];#pr(molt to maturity)
        T_szz  <- mp$T_yxszz[1,x,,,];        #size transition matrix
        S2_msz <- mp$S2_yxmsz[1,x,,,];       #survival after molting/mating
        n_xmsz[x,,,]<-calcEquilibriumSizeComp.TM(mc,R_z,S1_msz,P_sz,Th_sz,T_szz,S2_msz);
    }#x
    
    if (showPlot){
        mdfr<-melt(n_xmsz,value.name='val')
        pl <- ggplot(aes(x=z,y=val,fill=s),data=mdfr);
        pl <- pl + geom_bar(alpha=0.5,stat='identity',position='dodge');
        pl <- pl + labs(x='size (mm)',y='initial abundance')
        pl <- pl + guides(fill=guide_legend('shell condition'))
        pl <- pl + facet_wrap(~x+m,ncol=2);
        print(pl)
    }
    return(n_xmsz);
}

