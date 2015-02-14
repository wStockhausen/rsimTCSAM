#'
#'@title Calculate natural mortality.
#'
#'@param mc - model configuration object
#'
#'@return array M_yxmsz
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcNaturalMortality<-function(mc,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'TC','calcNaturalMortality');
    }
    
    d<-mc$dims;
    p<-mc$params$nm;
    
    M_yxmsz <- dimArray(mc,'y.x.m.s.z');
    mdfr<-NULL;
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        for (y in yrs) {
            for (x in d$x$nms) {
                for (m in d$m$nms) {
                    for (s in d$s$nms) {
                        M_yxmsz[y,x,m,s,] <- tb$M0_xms[x,m,s];
                    }
                }
            }
        }
        mdfrp<-melt(tb$M0_xms,value.name='val');
        mdfrp$fac<-paste(mdfrp$x,mdfrp$m,mdfrp$s,sep=' ');
        mdfrp$tb<-t;
        mdfr<-rbind(mdfr,mdfrp);
    }
    
    if (showPlot){
        pl <- ggplot(aes(x=x,y=val,fill=m),data=mdfr)
        pl <- pl + geom_bar(stat='identity',position='dodge')
        pl <- pl + labs(x='sex',y='M0')
        pl <- pl + guides(fill=guide_legend('maturity'))
        pl <- pl + facet_grid(tb~s)
        print(pl);
    }
    return(M_yxmsz);
}
