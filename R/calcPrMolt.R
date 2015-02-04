#'
#'@title Calculate probability at size of molting for immature crab by sex, shell condition.
#'
#'@title Function to calculate probability at size of molting for immature crab by sex, shell condition.
#'
#'@param mc - model configuration object
#'@param showPlot - flag to show plots
#'
#'@return prMolt_yxsz
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcPrMolt<-function(mc,showPlot=TRUE){
    d<-mc$dims;
    p<-mc$params$molting;
    
    if (mc$type=='TC'){
        prMolt_yxsz <- dimArray(mc,'y.x.s.z',val=0);    
        mdfr<-NULL;
        for (t in names(p$blocks)){
            tb<-p$blocks[[t]];
            yrs<-as.character(tb$years);
            for (x in d$x$nms){
                for (s in d$s$nms) {
                    z50 <- tb$z50_xms[x,m,s];
                    sdv <- tb$sdv_xms[x,m,s];
                    mp_z<-dimArray(mc,'z');
                    mp_z[] <- 1.0 - plogis(d$z$vls,z50,sdv);
                    mdfrp<-melt(mp_z,value.name='val');
                    mdfrp$fac<-paste(x,s,sep=', ');
                    mdfrp$t<-t;
                    mdfr<-rbind(mdfr,mdfrp);
                    for (y in yrs) prMolt_yxsz[y,x,s,]<-mp_z;
                }#s
            }#x
        }#t
    } else {
        throwModelTypeError(mc$type,'calcPrMolt()');
    }
    
    if (showPlot){
        p <- ggplot(aes(x=z,y=val,color=fac),data=mdfr)
        p <- p + geom_line()
        p <- p + labs(x='size (mm)',y='pr(molt|size)')
        p <- p + guides(color=guide_legend('sex, shell condition'));
        p <- p + facet_wrap(~t,ncol=1);
        print(p);
    }
    
    return(prMolt_yxsz);
}