#'
#'@title Calculate natural mortality.
#'
#'@param mc - model configuration object
#'
#'@return list  with array M_yxmsz and M_cxm
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
    p<-mc$params$natmort;
    
    M_cxm   <- dimArray(mc,'pc_natmort.x.m');
    M_yxmsz <- dimArray(mc,'y.x.m.s.z');
    mdfr<-NULL;
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        for (y in yrs) {
            for (x in d$x$nms) {
                for (m in d$m$nms) {
                    M_cxm[t,x,m]   <- tb$M0_xm[x,m];
                    for (s in d$s$nms) {
                        M_yxmsz[y,x,m,s,] <- M_cxm[t,x,m];
                    }
                }
            }
        }
    }
    
    if (showPlot){
        mdfr<-melt(M_cxm,value.name='val');
        pl <- ggplot(aes(x=x,y=val,fill=m),data=mdfr)
        pl <- pl + geom_bar(stat='identity',position='dodge')
        pl <- pl + labs(x='sex',y='M0')
        pl <- pl + guides(fill=guide_legend('maturity'))
        pl <- pl + facet_grid(pc~.)
        print(pl);
    }
    return(list(M_cxm=M_cxm,M_yxmsz=M_yxmsz));
}
