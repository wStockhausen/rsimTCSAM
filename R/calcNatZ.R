#'
#'@title Project population numbers at size forward through the model time interval
#'
#'@description Function to project population numbers at size forward through the model time interval
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param iN_xmsz - initial numbers at size array
##'
#'@return list with two elements:
#'MB_yx   - 2d array with mature biomass at mating time by year/sex
#'N_yxms  - 4-d array of population numbers by year/sex/maturity/shell condition
#'B_yxms  - 4-d array of population biomass by year/sex/maturity/shell condition
#'N_yxmsz - 5-d array of population numbers by year/sex/maturity/shell condition/size
#'
#'@details None.
#'
#'@import reshape2
#'@import ggplot2
#'@import wtsUtilities
#'
#'@export
#'
calcNatZ<-function(mc,mp,iN_xmsz,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'TC','calcNatZ');
    }
    
    d<-mc$dims;
    
    #calculate time series of population abundance
    N_yxmsz <- dimArray(mc,'y.x.m.s.z');#numbers-at-size
    MB_yx   <- dimArray(mc,'y.x',NA);   #mature biomass at time of mating
    N_yxmsz[1,,,,] <- iN_xmsz;
    for (y in (1:(d$y$n-1))){#note that y index is an integer, not a string
        MB_yx[y,]<-0;
        for (x in d$x$nms){
            N_msz <- dimArray(mc,'m.s.z');
            N_msz[,,] <- N_yxmsz[y,x,,,];#sex-specifc abundance at start of year y
            #calculate mature biomass at time of mating for year y
            for (s in d$s$nms){
                m<-'mature';
                MB_yx[y,x] <- MB_yx[y,x] + sum(mp$W_yxmsz[y,x,m,s,] * mp$S1_yxmsz[y,x,m,s,] * N_msz[m,s,]);
            }#s
            #project population one year forward
            R_z    <- mp$R_list$R_yxz[y,x,];     #recruitment
            S1_msz <- mp$S1_yxmsz[y,x,,,];       #survival to mating/molting
            P_sz   <- mp$prMolt_yxsz[y,x,,];     #pr(molt)
            Th_sz  <- mp$prM2M_yxsz[y,x,,];      #pr(molt to maturity)
            T_szz  <- mp$T_list$T_yxszz[y,x,,,]; #size transition matrix
            S2_msz <- mp$S2_yxmsz[y,x,,,];       #survival after molting/mating
            N_yxmsz[y+1,x,,,]<-runOneYear.TM(mc,R_z,S1_msz,P_sz,Th_sz,T_szz,S2_msz,N_msz);
        }#x
    }#y
    
    #calculate aggregate numbers/biomss (millions, 1000s t)
    N_yxms <- dimArray(mc,'y.x.m.s');#numbers
    B_yxms <- dimArray(mc,'y.x.m.s');#biomass
    for (y in d$y$nms){
        for (x in d$x$nms){
            for (m in d$m$nms){
                for (s in d$s$nms){
                    N_yxms[y,x,m,s]<-sum(N_yxmsz[y,x,m,s,]);
                    B_yxms[y,x,m,s]<-sum(mp$W_yxmsz[y,x,m,s,] * N_yxmsz[y,x,m,s,]);
                }#s
            }#m
        }#x
    }#y
    
    if (showPlot){
        mdfr<-melt(N_yxmsz,value.name='val');
        ddfr<-dcast(mdfr,x+y~.,fun.aggregate=sum,value.var='val');
        p <- ggplot(aes(x=y,y=`.`,color=x,shape=x),data=ddfr);
        p <- p + geom_line(alpha=0.8,size=2);
        p <- p + geom_point(alpha=0.8);
        p <- p + labs(x='year',y='Population Abundance');
        p <- p + guides(color=guide_legend('',order=1,alpha=1),
                        shape=guide_legend('',order=1));
        print(p);
        
        mdfrp<-melt(MB_yx,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=x,shape=x),data=mdfrp);
        p <- p + geom_line(alpha=0.8,size=1);
        p <- p + geom_point(alpha=0.8);
        p <- p + labs(x='year',y='Mature Biomass (at time of mating)');
        p <- p + guides(color=guide_legend('',order=1,alpha=1),
                        shape=guide_legend('',order=1));
        print(p);
        
        #size comps
        p <- ggplot(aes(x=y,y=z,color=val,size=val),data=mdfr);
        p <- p + geom_point(alpha=0.4);
        p <- p + scale_size_area(max_size=6);
        p <- p + scale_color_gradientn(colours=createColorPalette('jet',100,alpha=0.4))
        p <- p + labs(x='year',y='size (mm)',title='Population Abundance');
        p <- p + guides(size=guide_legend('millions',order=1),
                        color=guide_colorbar('',order=2,alpha=1));
        p <- p + facet_grid(m + s ~  x);
        print(p);
    }
    
    return(list(MB_yx=MB_yx,N_yxms=N_yxms,B_yxms=B_yxms,N_yxmsz=N_yxmsz));
}

#'
#'@title Project one sex of the population one time step forward.
#'
#'@description Function to project one sex of the population one time step forward.
#'
#'@param mc     - model configuration list object
#'@param R_z    - vector of recruits-at-size
#'@param S1_msz - 3d array of pr(survival) from start of year to mating/molting by maturity state, shell condition, size class
#'@param P_sz   - 2d array of the probability by size class of molting for immature crab, by shell condition
#'@param Th_sz  - 2d array with pr(molt to maturity|size, molt) for immature crab by shell condition
#'@param T_szz  - 3d array with size transition matrix for growth by immature crab by shell condition
#'@param S2_msz - 3d array of pr(survival) from mating/molting to end of year by maturity state, shell condition, size class
#'@param n_msz  - 3d array with initial nubers at size
#'
#'@return np_msz, a 3d array with projected numbers-at-size
#'
#'@export
#'
runOneYear.TM<-function(mc,R_z,S1_msz,P_sz,Th_sz,T_szz,S2_msz,n_msz){
    #create an identity matrix
    I <- diag(mc$dims$z$n);
    
    #calc the state transition matrices
    l<-calcStateTransitionMatrices(mc,S1_msz,P_sz,Th_sz,T_szz,S2_msz);
    
    #get initial numbers-at-size
    imm.ns <- n_msz['immature','new shell',];#immature, new shell
    imm.os <- n_msz['immature','old shell',];#immature, old shell
    mat.ns <- n_msz[  'mature','new shell',];#  mature, new shell
    mat.os <- n_msz[  'mature','old shell',];#  mature, old shell
    
    #project forward one time step
    np_msz<-dimArray(mc,'m.s.z');
    np_msz['immature','new shell',] <- l$A %*% imm.ns + l$B %*% imm.os + R_z;#immature, new shell
    np_msz['immature','old shell',] <- l$C %*% imm.ns + l$D %*% imm.os;      #immature, old shell
    np_msz[  'mature','new shell',] <- l$E %*% imm.ns + l$F %*% imm.os;      #  mature, new shell
    np_msz[  'mature','old shell',] <- l$G %*% mat.ns + l$H %*% mat.os;      #  mature, old shell
    
    return(np_msz);
}

