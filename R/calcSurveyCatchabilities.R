#'
#'@title Calculate size-specific survey catchabilities.
#'
#'@description Function to calculate size-specific survey catchabilities.
#'
#'@param mc - model configuration object
#'@param sel_cz - selectivity functions
#'@param showPlot - flag to show plots
#'
#'@return list with the following elements:\cr
#'sel_vyxmsz - size-specific selectivity\cr
#'devs_vy  - ln-scale annual deviations in Q\cr
#'Q_vxy    - fully-selected survey catchability\cr
#'Q_vyxms  - fully-selected survey catchability\cr
#'Q_vyxmsz - size-specific survey catchability
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcSurveyCatchabilities<-function(mc,sel_cz,showPlot=TRUE){
    d<-mc$dims;
    vs<-mc$params$surveys
    
    devs_vy   <-dimArray(mc,'v.y',val=NA);    #ln-scale capture rate devs by year for survey v
    Q_vxy     <-dimArray(mc,'v.x.y',val=NA);  #sex-specific capture rates by year for survey v
    Q_vyxms   <-dimArray(mc,'v.y.x.m.s',val=NA);
    sel_vyxmsz<-dimArray(mc,'v.y.x.m.s.z',val=NA)
    Q_vyxmsz  <-dimArray(mc,'v.y.x.m.s.z',val=NA);
    
    for (v in names(vs)){
        sel_xyz<-dimArray(mc,'x.y.z',val=NA);#sex-specific capture selectivity by year for survey v
        blocks<-vs[[v]]$blocks;
        mdfr<-list();
        for (t in names(blocks)){
            b<-blocks[[t]];
            yrs<-as.character(b$years);
            ps<-b$params;
            sel_xz<-dimArray(mc,'x.z',val=NA)
            for (x in d$x$nms){
                if (x=='male'){
                    devs<-rnorm(length(yrs),mean=0,sd=ps[[x]]$sdQ);
                    devs<-devs-mean(devs);#enforce sum to zero
                    devs_vy[v,yrs]<-devs;
                }
                Q_vxy[v,x,yrs]<-exp(ps[[x]]$lnQ+devs_vy[v,yrs]);#annual Q's in time block
                
                #get selectivity curve
                sel_xz[x,]<-sel_cz[ps[[x]]$idSel,];
                for (y in yrs) {sel_xyz[x,y,]<-sel_xz[x,];}
            }#x
            mdfrp<-melt(sel_xz,value.name='val');
            mdfrp$t<-t;
            mdfr<-rbind(mdfr,mdfrp);
            for (y in yrs){
                for (x in d$x$nms){
                    for (m in d$m$nms){
                        for (s in d$s$nms) {
                            #no dependence on m,s assumed
                            Q_vyxms[v,y,x,m,s]     <- Q_vxy[v,x,y];
                            sel_vyxmsz[v,y,x,m,s,] <- sel_xyz[x,y,];
                            Q_vyxmsz[v,y,x,m,s,]   <- Q_vxy[v,x,y]*sel_xyz[x,y,];
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
        mdfr<-melt(Q_vxy,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=x),data=mdfr);
        p <- p + geom_line();
        p <- p + labs(x='year',y='Catchability',title=v)
        p <- p + guides(color=guide_legend(''))
        p <- p + facet_grid(v~.)
        print(p)
    }
    
    return(list(sel_vyxmsz=sel_vyxmsz,devs_vy=devs_vy,
                Q_vxy=Q_vxy,Q_vyxms=Q_vyxms,Q_vyxmsz=Q_vyxmsz))
}