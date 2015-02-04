# Notes on model configuration object
# 
# mc$dims
#     y
#         mny - min year
#         mxy - max year
#         n   - number of model years to simulate
#         nms - vector of years as names
#         vls - vector of years as values
#     x
#         n   - number of sexes
#         nms - names of sexes
#     m
#         n   - number of maturity states
#         nms - names of maturity states
#     s
#         n   - number of shell conditions
#         nms - names of shell conditions
#     z
#         n - number of size bins
#         nms - size bins (midpoints) as names
#         vls - size bins as values
#     zp   (same as z, but used for post-molt bins in growth matrices)
#         n - number of size bins
#         nms - size bins (midpoints) as names
#         vls - size bins as values
#     zc
#         n - number of size cutpoints
#         nms - sizebin cutpoints as names
#         vls - sizebin cutpoints as values
#     f
#         n - number of fisheries
#         nms - names of fisheries
#     v
#         n - number of surveys
#         nms - names of surveys
# mc$params
#     recruitment
#         lnR  - ln-scale mean recruitment
#         sigR - ln-scale recruitment standard deviation
#         lnXR - nominal sex ratio
#         sigXR - ln-scale standard deviation for sex ratio deviations
#         aZ    - alpha parameter for rec. size distribution
#         bZ    - beta parameter for rec. size distribution
#-----------------------------------------------------------------------------------
#'
#'@title Model configuration for Tanner crab.
#'
#'@export
#'
ModelConfiguration.TCSAM<-function(){
    #-----dimensions
    mny=1950;mxy=2014;yrs<-mny:mxy;
    zcs<-seq(from=25,to=185,by=5);
    zbs<-midpoints(zcs);
    dims<-list(y=list(n=length(yrs),nms=as.character(yrs),vls=yrs,mny=mny,mxy=mxy),
               x=list(n=2,nms=c('male','female')),
               m=list(n=2,nms=c('immature','mature')),
               s=list(n=2,nms=c('new shell','old shell')),
               z=list(n=length(zbs),nms=as.character(zbs),vls=zbs),
               zp=list(n=length(zbs),nms=as.character(zbs),vls=zbs),
               zc=list(n=length(zcs),nms=as.character(zcs),vls=zcs),
               f=list(n=1,nms=c('Tanner crab directed fishery')),
               v=list(n=1,nms=c('NMFS trawl survey'))
               );
    #---------------------------------------------------------------
    #-----parameters
    params <- list();
    
    #weight-at-size
    nt<-1;#number of time blocks
    blocks<-list();
    #start of blocks
        #--time block 1
        years<-mny:mxy;
        a<-dimArray(list(dims=dims),'x.m')
        b<-dimArray(list(dims=dims),'x.m')
        a	b
        a['female','immature']<-0.000637;	b['female','immature']<-2.794;
        a['female',  'mature']<-0.000344	b['female',  'mature']<-2.956
        a[  'male','immature']<-0.000163	b[  'male','immature']<-3.136
        a[  'male',  'mature']<-0.000163	b[  'male',  'mature']<-3.136    
        blocks[[1]]<-list(years=years,
                          a_xm=a,
                          b_xm=b
                         );
    #end of blocks
    params$wAtZ<-list(blocks=blocks);

    #natural mortality
    nt<-2; #number of time blocks
    blocks<-list();
    #start of blocks
        #--time block 1
        years<-c(mny:1989,1995:mxy);
        M0_xms<-dimArraylist(dims=dims),'x.m.s');
        cvM_xms<-dimArraylist(dims=dims),'x.m.s');
        M0<-0.23;
        M0_xms[  'male',  'mature','new shell']<-M0;
        M0_xms[  'male',  'mature','old shell']<-M0;
        M0_xms[  'male','immature','new shell']<-M0;
        M0_xms[  'male','immature','old shell']<-M0;
        M0_xms['female',  'mature','new shell']<-M0;
        M0_xms['female',  'mature','old shell']<-M0;
        M0_xms['female','immature','new shell']<-M0;
        M0_xms['female','immature','old shell']<-M0;
        cvM_xms[,,] <- 0.1;
        blocks[[1]]<-list(years=years,
                          M0_xms=M0_xms,
                          cvM_xms=cvM_xms
                         );
        #--time block 2
        years<-c(1990:1994);
        M0_xms<-dimArraylist(dims=dims),'x.m.s');
        cvM_xms<-dimArraylist(dims=dims),'x.m.s');
        M0<-0.3;
        M0_xms[  'male',  'mature','new shell']<-M0;
        M0_xms[  'male',  'mature','old shell']<-M0;
        M0_xms[  'male','immature','new shell']<-M0;
        M0_xms[  'male','immature','old shell']<-M0;
        M0_xms['female',  'mature','new shell']<-M0;
        M0_xms['female',  'mature','old shell']<-M0;
        M0_xms['female','immature','new shell']<-M0;
        M0_xms['female','immature','old shell']<-M0;
        cvM_xms[,,] <- 0.1;
        blocks[[2]]<-list(years=years,
                          M0_xms=M0_xms,
                          cvM_xms=cvM_xms
                         );
    #end of blocks
    params$nm<-list(blocks=blocks);
    
    #molting: probability at size of molting, parameterized as a descending logistic
    nt<-1; #number of time blocks
    blocks<-list();
    #start of blocks
        #--time block 1
        years<-c(mny:mxy);
        z50_xms<-dimArray(list(dims=dims),'x.m.s');
        cv_xms<-dimArray(list(dims=dims),'x.m.s');
        z50_xms[,'immature',] <- 1000;#all immature crab molt
        z50_xms[,  'mature',] <--1000;#no mature crab molt
        cv_xms[,,]            <- 1;  #pretty steep
        blocks[[1]]<-list(years=years,
                          z50_xms=z50_xms,
                          cv_xms=cv_xms
                          );
    #end of blocks
    params$molting<-list(blocks=blocks);
    
    #molt to maturity: probability that molt is to maturity (terminal molt), at size.
    #parameterized as an increasing logistic function of size
    nt<-1; #number of time blocks
    blocks<-list();
    #start of blocks
        #--time block 1
        years<-c(mny:mxy);
        z50_xms<-dimArray(list(dims=dims),'x.m.s');
        cv_xms<-dimArray(list(dims=dims),'x.m.s');
        z50_xms[,'immature',] <- 1000;#all immature crab molt
        z50_xms[,  'mature',] <--1000;#no mature crab molt
        cv_xms[,,]            <- 1;  #pretty steep
        blocks[[2]]<-list(years=years,
                          z50_xms=z50_xms,
                          cv_xms=cv_xms
                          );
    #end of blocks
    params$moltToMaturity<-list(blocks=blocks);
    
    #growth
    nt<-1; #number of time blocks
    blocks<-list();
    #start of blocks
        #--time block 1
        years<-c(mny:mxy);
        a<-dimArray(list(dims=dims),'x');
        a[c('male','female')] <- exp(c(-0.798507696,-0.597837001));
        b<-dimArray(list(dims=dims),'x');
        b[c('male','female')] <- exp(c(-0.051293294,-0.105360516));
        scale<-dimArray(list(dims=dims),'x');
        scl[c('male','female')] <- exp(-0.287682072);
        blocks[[1]]<-list(a_x=a,
                          b_x=b,
                          s_x=scl
                          );
    #end of blocks
    params$growth<-list(blocks=blocks);

    #recruitment
    nt<-1; #number of time blocks
    blocks<-list();
    #start of blocks
        #--time block 1
        years<-c(mny:mxy);
        blocks[[1]]<-list(lnR     = 4.3,        # ln-scale mean recruitment
                          lnRCV   =-0.43275213, #ln-scale value for recruitment cv
                          lnXR    = 0,          #ln-scale nominal sex ratio
                          lnSigXR = 0,          #ln-scale standard deviation for sex ratio deviations
                          lnAlphaZ= 2.442347,   #ln-scale alpha parameter for rec. size distribution
                          lnBetaZ = 1.386294    #ln-scale beta parameter for rec. size distribution
                          )
    #end of blocks
    params$rec<-list(blocks=blocks);
    
    #fisheries
    f1<-list(name='Tanner crab directed fishery',
             blocks=list(
                        list(years=1965:2014,hm=0.3,mnF=0.3,sdF=0.4,offFX=0.1,
                             sel=list(male  =list(type='logistic',params=list(mu=100,sd=20)),
                                      female=list(type='logistic',params=list(mu= 60,sd=20))),
                             ret=list(male  =list(type='logistic',params=list(mu=140,sd=5)))
                        )
                    )
            );
    params$fisheries<-list(`Tanner crab directed fishery`=f1);
    
    #surveys
    s1<-list(name='NMFS trawl survey',
             blocks=list(
                        list(years=1975:2014,mnQ=0.7,sdQ=0.1,offQX=0.8,
                             sel=list(male  =list(type='logistic',params=list(mu=50,sd=20)),
                                      female=list(type='logistic',params=list(mu=40,sd=20)))
                        )
                    )
            );
    params$surveys<-list(`NMFS trawl survey`=s1);

    #-----model configuration
    mc<-list(type='TC',dims=dims,params=params)
    return(mc)
}
#mc<-ModelConfiguration.TCSAM();
