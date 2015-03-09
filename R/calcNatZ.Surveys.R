#'
#'@title Calculate survey numbers at size through the model time interval
#'
#'@description Function to calculate survey numbers at size through the model time interval
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param N_yxmsz - initial numbers at size array
#'
#'@return list with the following elements:
#'N_vyxms: 5-d array of survey numbers by year/sex/maturity/shell condition
#'B_vyxms: 5-d array of survey biomass by year/sex/maturity/shell condition
#'N_vyxmsz: 6-d array of survey numbers by year/sex/maturity/shell condition/size
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
    
    #--calculate survey abundance/biomass aggregated over size
    N_vyxms<-dimArray(mc,'v.y.x.m.s'); #survey numbers by y, x, m, s
    B_vyxms<-dimArray(mc,'v.y.x.m.s'); #survey biomass by y, x, m, s
    W_yxmsz<-mp$W_yxmsz;               #weight-at-size
    for (v in d$v$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        N_vyxms[v,y,x,m,s]<-N_vyxms[v,y,x,m,s]+sum(N_vyxmsz[v,y,x,m,s,]);
                        B_vyxms[v,y,x,m,s]<-B_vyxms[v,y,x,m,s]+sum(W_yxmsz[y,x,m,s,]*N_vyxmsz[v,y,x,m,s,]);
                    }#s
                }#m
            }#x
        }#y
    }#v
    
    if (showPlot){
        #abundance
        mdfr<-melt(N_vyxms,value.name='val');
        ddfr<-dcast(mdfr,v+y+x~.,fun.aggregate=sum,value.var='val')
        p <- ggplot(aes(x=y,y=`.`,color=x,shape=x),data=ddfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Survey Abundance (millions)")
        p <- p + guides(color=guide_legend('sex',order=1),
                        shape=guide_legend('sex'))
        p <- p + facet_grid(v~.)
        print(p);
        
        #biomass
        mdfr<-melt(B_vyxms,value.name='val');
        ddfr<-dcast(mdfr,v+y+x~.,fun.aggregate=sum,value.var='val')
        p <- ggplot(aes(x=y,y=`.`,color=x,shape=x),data=ddfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Survey Biomass (1000s t)")
        p <- p + guides(color=guide_legend('sex',order=1),
                        shape=guide_legend('sex'))
        p <- p + facet_grid(v~.)
        print(p);
        
        #size comps
        mdfr<-melt(N_vyxmsz,value.name='val');
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
    
    return(invisible(list(N_vyxms=N_vyxms,B_vyxms=B_vyxms,N_vyxmsz=N_vyxmsz)));
}