#'
#'@title Calculate size-specific survey catchabilities.
#'
#'@description Function to calculate size-specific survey catchabilities.
#'
#'@param mc - model configuration object
#'@param showPlot - flag to show plots
#'
#'@return list with the following elements:
#'sel_vyxmsz - size-specific selectivity
#'Q_vyxms - fully-selected survey catchability
#'Q_vyxmsz - size-specific survey catchability
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcSurveyCatchabilities<-function(mc,showPlot=TRUE){
    d<-mc$dims;
    vs<-mc$params$surveys
    
    Q_vyx     <-dimArray(mc,'v.y.x',val=NA);  #sex-specific capture rates by year for survey v
    devs_vy   <-dimArray(mc,'v.y',val=NA);
    Q_vyxms   <-dimArray(mc,'v.y.x.m.s',val=NA);
    sel_vyxmsz<-dimArray(mc,'v.y.x.m.s.z',val=NA)
    Q_vyxmsz  <-dimArray(mc,'v.y.x.m.s.z',val=NA);
    
    for (v in names(vs)){
        sel_yxz<-dimArray(mc,'y.x.z',val=NA);#sex-specific capture selectivity by year for survey v
        blocks<-vs[[v]]$blocks;
        mdfr<-list();
        for (t in names(blocks)){
            b<-blocks[[t]];
            yrs<-as.character(b$years);
            devs<-rnorm(length(yrs),mean=0,sd=b$sdQ);
            devs<-devs-mean(devs);#enforce sum to zero
            devs_vy[v,yrs]<-devs;
            Q_b<-b$mnQ*exp(-(b$sdQ^2)/2+devs);#annual Q's in time block for males
            sel_xz<-dimArray(mc,'x.z');#sex-specific capture selectivity for time block
            for (x in d$x$nms){
                #set catchabilities
                fac<-1;
                if (x=='female') fac<-b$offQX;
                Q_vyx[v,yrs,x]<-fac*Q_b;
                
                #calc selectivity/retention curves
                si<-b$sel[[x]];#selectivity info
                sel_xz[x,]<-calcSelectivity(si$type,d$z$vls,si$params);
                for (y in yrs) {sel_yxz[y,x,]<-sel_xz[x,];}
            }#x
            mdfrp<-melt(sel_xz,value.name='val');
            mdfrp$t<-t;
            mdfr<-rbind(mdfr,mdfrp);
            for (y in yrs){
                for (x in d$x$nms){
                    for (m in d$m$nms){
                        for (s in d$s$nms) {
                            #no dependence on m,s assumed
                            Q_vyxms[v,y,x,m,s]     <- Q_vyx[v,y,x];
                            sel_vyxmsz[v,y,x,m,s,] <- sel_yxz[y,x,];
                            Q_vyxmsz[v,y,x,m,s,]   <- Q_vyx[v,y,x]*sel_yxz[y,x,];
                        }#s
                    }#m
                }#x
            }#y
        }#t
        if (showPlot){
            p <- ggplot(aes(x=z,y=val,color=x,shape=x),data=mdfr);
            p <- p + geom_point(size=5,alpha=0.5);
            p <- p + geom_line();
            p <- p + labs(x='size (mm)',y='survey selectivity',title=v)
            p <- p + guides(color=guide_legend(''),shape=guide_legend(''))
            p <- p + facet_grid(t~.)
            print(p)
        }
    }#v
    if (showPlot){
        mdfr<-melt(Q_vyx,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=x),data=mdfr);
        p <- p + geom_line();
        p <- p + labs(x='year',y='Catchability',title=v)
        p <- p + guides(color=guide_legend(''))
        p <- p + facet_grid(v~.)
        print(p)
    }
    
    return(list(sel_vyxmsz=sel_vyxmsz,devs_vy=devs_vy,
                Q_vyxms=Q_vyxms,Q_vyxmsz=Q_vyxmsz))
}