#'
#'@title Calculate size-specific fishery capture and retention rates.
#'
#'@description Function to calculate size-specific fishery capture and retention rates.
#'
#'@param mc - model configuration object
#'@param showPlot - flag to show plots
#'
#'@return list with the following elements
#'sel_fyxmsz - size-specific selectivity
#'ret_fyxmsz - size-specific retention
#'devs_fy    - ln-scale deviation in annual capture rates
#'hm_fy      - handling mortality rate
#'cpF_fyxms  - fully-selected fishing capture rate
#'cpF_fyxmsz - size-specific fishing capture rate
#'tmF_fyxmsz - size-specific fishing mortality rate
#'rmF_fyxmsz - size-specific retention mortality rate
#'dmF_fyxmsz - size-specific discard mortality rate
#'tmF_yxmsz  - total fishing mortality across all fisheries
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcFishingMortalities<-function(mc,showPlot=TRUE){
    d<-mc$dims;
    fs<-mc$params$fisheries;
    
    hm_fy     <-dimArray(mc,'f.y',val=0);
    devs_fy   <-dimArray(mc,'f.y',val=0);
    cpF_fyxms   <-dimArray(mc,'f.y.x.m.s',val=0);
    sel_fyxmsz<-dimArray(mc,'f.y.x.m.s.z',val=0)
    ret_fyxmsz<-dimArray(mc,'f.y.x.m.s.z',val=0);
    cpF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    tmF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    rmF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    dmF_fyxmsz <-dimArray(mc,'f.y.x.m.s.z',val=0);
    
    cpF_fyx   <-dimArray(mc,'f.y.x',val=NA);  #sex-specific capture rates by year for fishery f
    for (f in names(fs)){
        sel_yxz<-dimArray(mc,'y.x.z',val=NA);#sex-specific capture selectivity by year for fishery f
        ret_yxz<-dimArray(mc,'y.x.z',val=NA);#sex-specific retention by year for fishery f
        blocks<-fs[[f]]$blocks;
        mdfr<-list();
        for (t in names(blocks)){
            b<-blocks[[t]];
            yrs<-as.character(b$years);
            hm_fy[f,yrs]<-b$hm;
            devs<-rnorm(length(yrs),mean=0,sd=b$sdF);
            devs<-devs-mean(devs);#enforce sum to zero
            devs_fy[f,yrs]<-devs;
            F_b<-exp(b$lnF+devs);#annual F's in time block for males
            sel_xz<-dimArray(mc,'x.z');#sex-specific capture selectivity for time block
            ret_xz<-dimArray(mc,'x.z');#sex-specific retention for time block
            for (x in d$x$nms){
                #set capture rates
                fac<-1;
                if (x=='female') fac<-exp(b$lnFX);
                cpF_fyx[f,yrs,x]<-fac*F_b;
                
                #calc selectivity/retention curves
                si<-b$sel[[x]];#selectivity info
                ri<-b$ret[[x]];#retention info
                sel_xz[x,]<-calcSelectivity(si$type,d$z$vls,si$params);
                ret_xz[x,]<-0*sel_xz[x,];
                if (!is.null(ri)) ret_xz[x,]<-calcSelectivity(ri$type,d$z$vls,ri$params);
                for (y in yrs) {
                    sel_yxz[y,x,]<-sel_xz[x,];
                    ret_yxz[y,x,]<-ret_xz[x,];
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
                            cpF_fyxms[f,y,x,m,s]     <- cpF_fyx[f,y,x];
                            sel_fyxmsz[f,y,x,m,s,] <- sel_yxz[y,x,];
                            ret_fyxmsz[f,y,x,m,s,] <- ret_yxz[y,x,];
                            #fishing capture rates
                            cpF_fyxmsz[f,y,x,m,s,]  <- cpF_fyxms[f,y,x,m,s]*sel_fyxmsz[f,y,x,m,s,];
                            #retention mortality rates
                            rmF_fyxmsz[f,y,x,m,s,]  <- cpF_fyxmsz[f,y,x,m,s,]*(ret_fyxmsz[f,y,x,m,s,]);
                            #discard mortality rates
                            dmF_fyxmsz[f,y,x,m,s,]  <- cpF_fyxmsz[f,y,x,m,s,]*((1-ret_fyxmsz[f,y,x,m,s,])*b$hm);
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
        mdfr<-melt(cpF_fyx,value.name='val');
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
    
    return(list(sel_fyxmsz=sel_fyxmsz,ret_fyxmsz=ret_fyxmsz,devs_fy=devs_fy,hm_fy=hm_fy,
                cpF_fyx=cpF_fyx,cpF_fyxms=cpF_fyxms,cpF_fyxmsz=cpF_fyxmsz,
                tmF_yxmsz=tmF_yxmsz,tmF_fyxmsz=tmF_fyxmsz,
                rmF_fyxmsz=rmF_fyxmsz,dmF_fyxmsz=dmF_fyxmsz))
}
