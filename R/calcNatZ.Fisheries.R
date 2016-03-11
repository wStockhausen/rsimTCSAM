#'
#'@title Calculate fishery catches, mortality at size through the model time interval
#'
#'@description Function to calculate fishery catches, mortality at size through the model time interval
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param N_yxmsz - numbers at size array
#'
#'@return list with the following elements:
#'cpN_fyxms:  5-d array of fishery catches (NOT MORTALITY) by fishery/year/sex/maturity/shell condition
#'cpB_fyxms:  5-d array of fishery catches (NOT MORTALITY) by fishery/year/sex/maturity/shell condition
#'cpN_fyxmsz: 6-d array of fishery catches (NOT MORTALITY) by fishery/year/sex/maturity/shell condition/size
#'dsN_fyxms:  5-d array of fishery discard numbers (NOT MORTALITY) by fishery/year/sex/maturity/shell condition
#'dsB_fyxms:  5-d array of fishery discard biomass (NOT MORTALITY) by fishery/year/sex/maturity/shell condition
#'dsN_fyxmsz: 6-d array of fishery discard numbers (NOT MORTALITY) by fishery/year/sex/maturity/shell condition/size
#'tmN_fyxms:  5-d array of total fishery mortality as numbers by fishery/year/sex/maturity/shell condition
#'tmB_fyxms:  5-d array of total fishery mortality as biomass by fishery/year/sex/maturity/shell condition
#'tmN_fyxmsz: 6-d array of total fishery mortality as numbers by fishery/year/sex/maturity/shell condition/size
#'rmN_fyxms:  5-d array of retention mortality as numbers by fishery/year/sex/maturity/shell condition
#'rmB_fyxms:  5-d array of retention mortality as biomass by fishery/year/sex/maturity/shell condition
#'rmN_fyxmsz: 6-d array of retention mortality by year/sex/maturity/shell condition/size
#'dmN_fyxms:  5-d array of discard mortality as numbers by fishery/year/sex/maturity/shell condition
#'dmB_fyxms:  5-d array of discard mortality as biomass by fishery/year/sex/maturity/shell condition
#'dmN_fyxmsz: 6-d array of discard mortality by fishery/year/sex/maturity/shell condition/size
#'
#'@import reshape2
#'@import ggplot2
#'@import wtsUtilities
#'
#'@export
#'
calcNatZ.Fisheries<-function(mc,mp,N_yxmsz,showPlot=TRUE){
    d<-mc$dims;
    W_yxmsz <- mp$W_yxmsz;           #weight-at-size retained
    M_yxmsz <- mp$M_yxmsz;           #natural mortality rate
    
    #calculate time series of fisheries catches
    cpF_fyxmsz <- mp$F_list$cpF_fyxmsz;#capture rates
    tmF_fyxmsz <- mp$F_list$tmF_fyxmsz;#mortality rates
    rmF_fyxmsz <- mp$F_list$rmF_fyxmsz;#retention mortality rates
    dmF_fyxmsz <- mp$F_list$dmF_fyxmsz;#discard mortality rates
    tmF_yxmsz  <- mp$F_list$tmF_yxmsz; #total mortality rates across all fisheries
    
    #calc survival to midpoint of fisheries
    S_yxmsz <- exp(-mc$params$fish.time*mp$M_yxmsz);
    
    fac_yxmsz  <- dimArray(mc, 'y.x.m.s.z');
    fac_yxmsz[,,,,]<-(1/tmF_yxmsz)*(1-exp(-tmF_yxmsz))*(S_yxmsz*N_yxmsz);
    
    cpN_fyxmsz <- dimArray(mc,'f.y.x.m.s.z');#numbers captured
    dsN_fyxmsz <- dimArray(mc,'f.y.x.m.s.z');#numbers discarded
    rmN_fyxmsz <- dimArray(mc,'f.y.x.m.s.z');#numbers retained
    dmN_fyxmsz <- dimArray(mc,'f.y.x.m.s.z');#numbers discarded and killed
    tmN_fyxmsz <- dimArray(mc,'f.y.x.m.s.z');#numbers killed
    for (f in d$f$nms){
        cpN_fyxmsz[f,,,,,]<-cpF_fyxmsz[f,,,,,]*fac_yxmsz[,,,,];#numbers captured-at-size
        tmN_fyxmsz[f,,,,,]<-tmF_fyxmsz[f,,,,,]*fac_yxmsz[,,,,];#total mortality-at-size
        rmN_fyxmsz[f,,,,,]<-rmF_fyxmsz[f,,,,,]*fac_yxmsz[,,,,];#retained mortality-at-size
        dmN_fyxmsz[f,,,,,]<-dmF_fyxmsz[f,,,,,]*fac_yxmsz[,,,,];#discard mortality-at-size
        dsN_fyxmsz[f,,,,,]<-cpN_fyxmsz[f,,,,,]-rmN_fyxmsz[f,,,,,];#numbers discarded-at-size
    }#f
    
    #--Captured abundance/biomass (1000s mt) aggregated over size [NOT mortality]
    cpN_fyxms<-dimArray(mc,'f.y.x.m.s');     #captured abundance by f,y,x,m,s
    cpB_fyxms<-dimArray(mc,'f.y.x.m.s');     #captured biomass by f,y,x,m,s
    for (f in d$f$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        cpN_fyxms[f,y,x,m,s]<-cpN_fyxms[f,y,x,m,s]+sum(cpN_fyxmsz[f,y,x,m,s,]);
                        cpB_fyxms[f,y,x,m,s]<-cpB_fyxms[f,y,x,m,s]+sum(W_yxmsz[y,x,m,s,]*cpN_fyxmsz[f,y,x,m,s,]);
                    }#s
                }#m
            }#x
        }#y
    }#f
    
    #--Retained catch abundance/biomass (1000s mt) aggregated over size
    rmN_fyxms<-dimArray(mc,'f.y.x.m.s');     #retained abundance by f,y,x,m,s
    rmB_fyxms<-dimArray(mc,'f.y.x.m.s');     #retained biomass by f,y,x,m,s
    for (f in d$f$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        rmN_fyxms[f,y,x,m,s]<-rmN_fyxms[f,y,x,m,s]+sum(rmN_fyxmsz[f,y,x,m,s,]);
                        rmB_fyxms[f,y,x,m,s]<-rmB_fyxms[f,y,x,m,s]+sum(W_yxmsz[y,x,m,s,]*rmN_fyxmsz[f,y,x,m,s,]);
                    }#s
                }#m
            }#x
        }#y
    }#f
    
    #--Discarded catch abundance/biomass (millions/1000s mt) [NOT mortality] aggregated over size
    dsN_fyxms<-dimArray(mc,'f.y.x.m.s');  #discard numbers by f,y,x,m,s (NOT mortality)
    dsB_fyxms<-dimArray(mc,'f.y.x.m.s');  #discard biomass by f,y,x,m,s (NOT mortality)
    for (f in d$f$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        dsN_fyxms[f,y,x,m,s]<-dsN_fyxms[f,y,x,m,s]+sum(cpN_fyxmsz[f,y,x,m,s,]-rmN_fyxmsz[f,y,x,m,s,]);
                        dsB_fyxms[f,y,x,m,s]<-dsB_fyxms[f,y,x,m,s]+sum(W_yxmsz[y,x,m,s,]*(cpN_fyxmsz[f,y,x,m,s,]-rmN_fyxmsz[f,y,x,m,s,]));
                    }#s
                }#m
            }#x
        }#y
    }#f
    
    #--Discarded catch MORTALITY in abundance/biomass (millions/1000s mt) aggregated over size
    dmN_fyxms<-dimArray(mc,'f.y.x.m.s');  #discard mortality in numbers by f,y,x,m,s
    dmB_fyxms<-dimArray(mc,'f.y.x.m.s');  #discard mortality in biomass by f,y,x,m,s
    for (f in d$f$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        dmN_fyxms[f,y,x,m,s]<-dmN_fyxms[f,y,x,m,s]+sum(dmN_fyxmsz[f,y,x,m,s,]);
                        dmB_fyxms[f,y,x,m,s]<-dmB_fyxms[f,y,x,m,s]+sum(W_yxmsz[y,x,m,s,]*dmN_fyxmsz[f,y,x,m,s,]);
                    }#s
                }#m
            }#x
        }#y
    }#f
    
    #--total mortality numbers/biomass (millions, 1000's t) aggregated over size
    tmN_fyxms<-dimArray(mc,'f.y.x.m.s');  #total numbers killed by f,x,y,m,s
    tmB_fyxms<-dimArray(mc,'f.y.x.m.s');  #total biomass killed by f,x,y,m,s
    tmN_fyxms <- rmN_fyxms+dmN_fyxms;
    tmB_fyxms <- rmB_fyxms+dmB_fyxms;
    
    if (showPlot){
        #retained mortality
        mdfr<-melt(rmN_fyxms,value.name='val');
        ddfr<-dcast(mdfr,f+y~`.`,fun.aggregate=sum,value.var='val')
        p <- ggplot(aes(x=y,y=`.`),data=ddfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Retained Catch Abundance (millions)")
        p <- p + facet_grid(f~.)
        print(p);
        mdfr<-melt(rmB_fyxms,value.name='val');
        ddfr<-dcast(mdfr,f+y~`.`,fun.aggregate=sum,value.var='val')
        p <- ggplot(aes(x=y,y=`.`),data=ddfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Retained Catch Biomass (1000s t)")
        p <- p + facet_grid(f~.)
        print(p);
        
        #discard catch
        mdfr<-melt(dsN_fyxms,value.name='val');
        ddfr<-dcast(mdfr,f+y+x~`.`,fun.aggregate=sum,value.var='val')
        p <- ggplot(aes(x=y,y=`.`,color=x,shape=x),data=ddfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Discard Catch Abundance (millions)")
        p <- p + guides(color=guide_legend('sex',order=1),
                        shape=guide_legend('sex'))
        p <- p + facet_grid(f~.)
        print(p);
        mdfr<-melt(dsB_fyxms,value.name='val');
        ddfr<-dcast(mdfr,f+y+x~`.`,fun.aggregate=sum,value.var='val')
        p <- ggplot(aes(x=y,y=`.`,color=x,shape=x),data=ddfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Discard Catch Biomass (1000s t)")
        p <- p + guides(color=guide_legend('sex',order=1),
                        shape=guide_legend('sex'))
        p <- p + facet_grid(f~.)
        print(p);
        
        #more stuff
        ncdfr<-melt(cpN_fyxmsz,value.name='val'); ncdfr$type<-'captured';
        nmdfr<-melt(tmN_fyxmsz,value.name='val'); nmdfr$type<-'total mortality';
        rmdfr<-melt(rmN_fyxmsz,value.name='val'); rmdfr$type<-'retained mortality';
        dmdfr<-melt(dmN_fyxmsz,value.name='val'); dmdfr$type<-'discard mortality';
        mdfr<-rbind(ncdfr,nmdfr,rmdfr,dmdfr);
        ddfr<-dcast(mdfr,f+type+x+y~.,fun.aggregate=sum,value.var='val');
        p <- ggplot(aes(x=y,y=`.`,color=type,linetype=type,shape=type),data=ddfr);
        p <- p + geom_line(alpha=0.8,size=2);
        p <- p + geom_point(alpha=0.8);
        p <- p + labs(x='year',y='Fishery Catch/Mortality (millions)');
        p <- p + facet_grid(f ~ x);
        print(p);
        
        #size comps
        ddfr<-dcast(mdfr,f+type+x+y+z~.,fun.aggregate=sum,value.var='val');
        for (fp in d$f$nms){
            p <- ggplot(aes(x=y,y=z,color=`.`,size=`.`),data=ddfr[ddfr$f==fp,]);
            p <- p + geom_point(alpha=0.4);
            p <- p + scale_size_area(max_size=6);
            p <- p + scale_color_gradientn(colours=createColorPalette('jet',100,alpha=0.4))
            p <- p + labs(x='year',y='size (mm)',title=paste(fp,'Catch/Mortality (millions)'));
            p <- p + guides(size=guide_legend('millions',order=1),color=guide_colorbar(''));
            p <- p + facet_grid(type ~ x);
            print(p);
        }
    }
    
    return(list(cpN_fyxms=cpN_fyxms,cpB_fyxms=cpB_fyxms,cpN_fyxmsz=cpN_fyxmsz,
                dsN_fyxms=dsN_fyxms,dsB_fyxms=dsB_fyxms,dsN_fyxmsz=dsN_fyxmsz,
                tmN_fyxms=tmN_fyxms,tmB_fyxms=tmB_fyxms,tmN_fyxmsz=tmN_fyxmsz,
                rmN_fyxms=rmN_fyxms,rmB_fyxms=rmB_fyxms,rmN_fyxmsz=rmN_fyxmsz,
                dmN_fyxms=dmN_fyxms,dmB_fyxms=dmB_fyxms,dmN_fyxmsz=dmN_fyxmsz));
}