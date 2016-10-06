#'
#'@title Calculate size-specific fishery capture and retention rates
#'
#'@description Function to calculate size-specific fishery capture and retention rates.
#'
#'@param mc - model configuration object
#'@param sel_cz - selectivity functions
#'@param showPlot - flag to show plots
#'
#'@return list with the following elements
#'sel_fyxmsz - size-specific selectivity
#'ret_fyxmsz - size-specific retention
#'devs_fy    - ln-scale deviations in annual capture rates
#'hm_fxy     - handling mortality rate
#'cpF_fyxms  - fully-selected fishing capture rate
#'cpF_fyxmsz - size-specific fishing capture rate
#'tmF_fyxmsz - size-specific fishing mortality rate
#'rmF_fyxmsz - size-specific retention mortality rate
#'dmF_fyxmsz - size-specific discard mortality rate
#'tmF_yxmsz  - total fishing mortality across all fisheries
#'
#'@details None.
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcFishingMortalities<-function(mc,sel_cz,showPlot=TRUE){
    d<-mc$dims;
    fs<-mc$params$fisheries;
    
    hm_fxy    <-dimArray(mc,'f.x.y',val=0);
    devs_fy   <-dimArray(mc,'f.y',val=0);
    cpF_fyxms <-dimArray(mc,'f.y.x.m.s',val=0);
    sel_fyxmsz<-dimArray(mc,'f.y.x.m.s.z',val=0)
    ret_fyxmsz<-dimArray(mc,'f.y.x.m.s.z',val=0);
    cpF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    tmF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    rmF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    dmF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    
    cpF_fxy   <-dimArray(mc,'f.x.y',val=NA);  #sex-specific capture rates by year for fishery f
    for (f in names(fs)){
        sel_xyz<-dimArray(mc,'x.y.z',val=NA);#sex-specific capture selectivity by year for fishery f
        ret_xyz<-dimArray(mc,'x.y.z',val=NA);#sex-specific retention by year for fishery f
        blocks<-fs[[f]]$blocks;
        mdfr<-list();
        for (t in names(blocks)){
            b<-blocks[[t]];
            yrs<-as.character(b$years);
            ps<-b$params;
            sel_xz<-dimArray(mc,'x.z');#sex-specific capture selectivity for time block
            ret_xz<-dimArray(mc,'x.z');#sex-specific retention for time block
            for (x in d$x$nms){
                hm_fxy[f,x,yrs]<-ps[[x]]$hm;
                if (x=='male'){
                    devs<-rnorm(length(yrs),mean=0,sd=ps[[x]]$sdF);
                    devs<-devs-mean(devs);#enforce sum to zero
                    devs_fy[f,yrs]<-devs;
                }
                cpF_fxy[f,x,yrs]<-exp(ps[[x]]$lnF+devs_fy[f,yrs]);
                
                #calc selectivity/retention curves
                sel_xz[x,]<-sel_cz[ps[[x]]$idSel,];
                ret_xz[x,]<-0*sel_xz[x,];
                if (ps[[x]]$idRet>0) ret_xz[x,]<-sel_cz[ps[[x]]$idRet,];
                for (y in yrs) {
                    sel_xyz[x,y,]<-sel_xz[x,];
                    ret_xyz[x,y,]<-ret_xz[x,];
                }
            }#x
            sdfr<-melt(sel_xz,value.name='val');
            sdfr$type<-'selectivity';
            rdfr<-melt(ret_xz,value.name='val');
            rdfr$type<-'retention';
            mdfrp<-rbind(sdfr,rdfr)
            mdfrp$t<-t;
            mdfr<-rbind(mdfr,mdfrp);
            #print(dimnames(cpF_fyxms))
            for (y in yrs){
                for (x in d$x$nms){
                    for (m in d$m$nms){
                        for (s in d$s$nms) {
                            #no dependence on m,s assumed
                            #cat(f,y,x,m,s,'\n')
                            cpF_fyxms[f,y,x,m,s]   <- cpF_fxy[f,x,y];
                            sel_fyxmsz[f,y,x,m,s,] <- sel_xyz[x,y,];
                            ret_fyxmsz[f,y,x,m,s,] <- ret_xyz[x,y,];
                            #fishing capture rates
                            cpF_fyxmsz[f,y,x,m,s,]  <- cpF_fyxms[f,y,x,m,s]*sel_fyxmsz[f,y,x,m,s,];
                            #retention mortality rates
                            rmF_fyxmsz[f,y,x,m,s,]  <- cpF_fyxmsz[f,y,x,m,s,]*(ret_fyxmsz[f,y,x,m,s,]);
                            #discard mortality rates
                            dmF_fyxmsz[f,y,x,m,s,]  <- cpF_fyxmsz[f,y,x,m,s,]*((1-ret_fyxmsz[f,y,x,m,s,])*hm_fxy[f,x,y]);
                            #total fishing mortality rates
                            tmF_fyxmsz[f,y,x,m,s,]  <- rmF_fyxmsz[f,y,x,m,s,]+dmF_fyxmsz[f,y,x,m,s,];
                        }#s
                    }#m
                }#x
            }#y
        }#b
        if (showPlot){
            p <- ggplot(aes(x=z,y=val,color=x,shape=x,linetype=type),data=mdfr);
            p <- p + geom_point(size=5,alpha=0.5);
            p <- p + geom_line();
            p <- p + labs(x='size (mm)',y='selectivity/retention',title=f)
            p <- p + guides(color=guide_legend('',order=1),
                            shape=guide_legend('',order=1),
                            linetype=guide_legend('',order=2));
            p <- p + facet_grid(t~.)
            print(p)
        }
    }#f
    if (showPlot){
        mdfr<-melt(cpF_fxy,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=x),data=mdfr);
        p <- p + geom_line();
        p <- p + labs(x='year',y='Capture Rate')
        p <- p + guides(color=guide_legend(''))
        p <- p + facet_grid(f~.)
        print(p)
    }
    
    #calc total fishing mortality
    tmF_yxmsz<-dimArray(mc, 'y.x.m.s.z');
    for (f in d$f$nms){
        tmF_yxmsz <- tmF_yxmsz + tmF_fyxmsz[f,,,,,];
    }
    
    return(list(sel_fyxmsz=sel_fyxmsz,ret_fyxmsz=ret_fyxmsz,devs_fy=devs_fy,hm_fxy=hm_fxy,
                cpF_fxy=cpF_fxy,cpF_fyxms=cpF_fyxms,cpF_fyxmsz=cpF_fyxmsz,
                tmF_yxmsz=tmF_yxmsz,tmF_fyxmsz=tmF_fyxmsz,
                rmF_fyxmsz=rmF_fyxmsz,dmF_fyxmsz=dmF_fyxmsz))
}
