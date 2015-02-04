#'
#'@title Calculate natural mortality.
#'
#'@param mc - model configuration object
#'
#'@return array M_yxmsz
#'
#'@import ggplot2
#'
#'@export
#'
calcNaturalMortality<-function(mc,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'calcNaturalMortality');
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
        p <- ggplot(aes(x=tb,y=val,fill=fac),data=mdfr)
        p <- p + geom_bar(stat='identity',position='dodge')
        p <- p + labs(x='time block',y='M0')
        p <- p + guides(fill=guide_legend('sex, maturity, \nshell condition'))
        print(p);
    }
    return(M_yxmsz);
}
