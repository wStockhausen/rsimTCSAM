#'
#'@title Calculate probability-at-size (pre-molt) for immature crab that a molt is to maturity, by sex and shell condition.
#'
#'@title Function to calculate probability-at-size (pre-molt) for immature crab that a molt is to maturity, by sex and shell condition.
#'
#'@param mc - model configuration object
#'@param showPlot - flag to show plots
#'
#'@return prMoltToMat_yxsz
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcPrMoltToMaturity<-function(mc,showPlot=TRUE){
    d<-mc$dims;
    p<-mc$params$moltToMaturity;
    
    if (mc$type=='TC'){
        prMoltToMat_yxsz <- dimArray(mc,'y.x.s.z',val=0);    
        mdfr<-NULL;
        for (t in names(p$blocks)){
            tb<-p$blocks[[t]];
            yrs<-as.character(tb$years);
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms) {
                        z50 <- tb$z50_xms[x,m,s];
                        sdv <- tb$sdv_xms[x,m,s];
                        mp_z<-dimArray(mc,'z');
                        mp_z[] <- plogis(d$z$vls,z50,sdv);
                        mdfrp<-melt(mp_z,value.name='val');
                        mdfrp$fac<-paste(x,s,sep=', ');
                        mdfrp$t<-t;
                        mdfr<-rbind(mdfr,mdfrp);
                        for (y in yrs) prMoltToMat_yxsz[y,x,s,]<-mp_z;
                    }#s
                }#m      
            }#x
        }#t
    } else {
        throwModelTypeError(mc$type,'calcPrMoltToMaturity()');
    }
    
    if (showPlot){
        p <- ggplot(aes(x=z,y=val,color=fac),data=mdfr)
        p <- p + geom_line()
        p <- p + labs(x='size (mm)',y='pr(molt-to-maturity|size, molt)')
        p <- p + guides(color=guide_legend('sex, shell condition'));
        p <- p + facet_wrap(~t,ncol=1);
        print(p);
    }
    
    return(prMoltToMat_yxsz);
}