#'
#'@title Calculate weight-at-size by year, sex, maturity state and shell condition.
#'
#'@param mc - model configuration object
#'
#'@return W_yxmsz: 5d array with weight-at-size in KG by year/sex/maturity state/shell condition
#'
#'@details Input parameters are for weight-at-size in grams, but converted to kg so biomass
#'is in 1000s t.
#'
#'@import ggplot2
#'@import reshape2
#'
#'@export
#'
calcWatZ<-function(mc,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'TC','calcWatZ()');
    }
    
    d<-mc$dims;
    p<-mc$params$wAtZ;
    
    W_yxmsz <- dimArray(mc,'y.x.m.s.z',val=0);
    mdfr<-NULL;
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        W_xmz<-dimArray(mc,'x.m.z');
        for (x in d$x$nms){
            for (m in d$m$nms){
                lnA<-log(tb$a_xm[x,m]);#no dependence on maturity/shell condition
                B  <-tb$b_xm[x,m];
                W_xmz[x,m,] <- exp(lnA+B*log(d$z$vls));#weight in g
                for (s in d$s$nms){
                    for (y in yrs) {W_yxmsz[y,x,m,s,]<-W_xmz[x,m,]/1000;} #in kg
                }
            }
        }
        mdfrp<-melt(W_xmz,value.name='val');
        mdfrp$t<-t;
        mdfr<-rbind(mdfr,mdfrp);
    }
    if (showPlot){
        pz <- ggplot(mapping=aes(x=z,y=val,color=x,linetype=m),data=mdfr)
        pz <- pz + geom_line();
        pz <- pz + labs(x='size (mm)',y='weight (g)',title='')
        pz <- pz + guides(color=guide_legend('sex',order=1),linetype=guide_legend('maturity',order=2));
        pz <- pz + facet_wrap(~t,ncol=1);
        print(pz)
    }
    
    return(W_yxmsz);
}
