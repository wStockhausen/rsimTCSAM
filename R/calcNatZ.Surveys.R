#'
#'@title Calculate survey numbers at size through the model time interval
#'
#'@description Function to calculate survey numbers at size through the model time interval
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param N_yxmsz - initial numbers at size array
#'
#'@return N_vyxmsz: 6-d array of survey numbers by year/sex/maturity/shell condition/size
#'
#'@import reshape2
#'@import ggplot2
#'@import wtsUtilities
#'
#'@export
#'
calcNatZ.Surveys<-function(mc,mp,N_yxmsz,showPlot=TRUE){
    #calculate time series of survey abundance
    d<-mc$dims;
    Q_vyxmsz <- mp$S_list$Q_vyxmsz;
    N_vyxmsz <- dimArray(mc,'v.y.x.m.s.z');
    for (v in d$v$nms){
        N_vyxmsz[v,,,,,]<-Q_vyxmsz[v,,,,,]*N_yxmsz;
    }#v
    if (showPlot){
        mdfr<-melt(N_vyxmsz,value.name='val');
        ddfr<-dcast(mdfr,v+x+y~.,fun.aggregate=sum,value.var='val');
        p <- ggplot(aes(x=y,y=`.`,color=x,linetype=x,shape=x),data=ddfr);
        p <- p + geom_line(alpha=0.8,width=2);
        p <- p + geom_point(alpha=0.8);
        p <- p + labs(x='year',y='Survey Abundance (millions)');
        p <- p + guides(color=guide_legend('',order=1,alpha=1),
                        linetype=guide_legend('',order=1),
                        shape=guide_legend('',order=1));
        p <- p + facet_grid(v ~ .);
        print(p);
        
        #size comps
        for (vp in d$v$nms){
            p <- ggplot(aes(x=y,y=z,color=val,size=val),data=mdfr[mdfr$v==vp,]);
            p <- p + geom_point(alpha=0.4);
            p <- p + scale_size_area(max_size=6);
            p <- p + scale_color_gradientn(colours=createColorPalette('jet',100,alpha=0.4))
            p <- p + labs(x='year',y='size (mm)',title=paste(vp,'Survey Abundance'));
            p <- p + guides(size=guide_legend('millions',order=1),color=guide_colorbar(''));
            p <- p + facet_grid(m + s ~ x);
            print(p);
        }#vp
    }
    
    return(N_vyxmsz)
}