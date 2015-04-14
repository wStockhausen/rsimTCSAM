#'
#'@title Calculate values for a selectivity curve
#'
#'@param type - the type of selectivity function to calculate
#'@param z - vector of values at which to calculate the function
#'@param params - the selectivity function parameters, as a vector
#'
#'@return vector matching size of z, with names given by elements of z
#'
#'@export
#'
calcSelectivity<-function(type,z,params){
    fsz <- 0; #fully-selected size for re-scaling selectivity function
    if (tolower(type)=='asclogistic\n'){
        cat('sel function=asclogistic\n')
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic(z,params[1],params[2],fsz);
    } else if (tolower(type)=='asclogistic5095'){
        cat('sel function=asclogistic5095')
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic5095(z,params[1],params[2],fsz);
    } else if (tolower(type)=='asclogistic50D95'){
        cat('sel function=asclogistic5095')
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic50D95(z,params[1],params[2],fsz);
    } else if (tolower(type)=='asclogistic50Ln95'){
        if (length(params)>2) fsz<-params[3];
        res<-asclogistic50LnD95(z,params[1],params[2],fsz);
    } else if (tolower(type)=='asclogisticLn50LnD95'){
        if (length(params)>2) fsz<-params[3];
        res<-asclogisticLn50LnD95(z,params[1],params[2],fsz);
    } else if (tolower(type)=='dbllogistic'){
        if (length(params)>4) fsz<-params[5];
        res<-dbllogistic(z,params[1],params[2],params[3],params[4],fsz);
    } else if (tolower(type)=='dbllogisticLnD50'){
        if (length(params)>4) fsz<-params[5];
        res<-dbllogisticLnD50(z,params[1],params[2],params[3],params[4],fsz);
    } else if (tolower(type)=='dbllogistic50Ln95'){
        if (length(params)>4) fsz<-params[5];
        res<-dbllogistic50LnD95(z,params[1],params[2],params[3],params[4],fsz);
    } else if (tolower(type)=='dbllogisticLn50LnD95'){
        if (length(params)>4) fsz<-params[5];
        res<-dbllogisticLn50LnD95(z,params[1],params[2],params[3],params[4],fsz);
    } else {
        cat('Selectivity/retention function type "',type,'" not recognnized.\n',sep='');
        cat('Aborting...\n');
        stop();
    }
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate the logistic function
#'
#'@description Function to calculate the logistic function
#'
#'@param z    - vector of sizes at which to compute selectivities
#'@param z50 - size at which selectivity  = 0.5 (logit-scale mean)
#'@param sd  - standard deviation in selectivity (logit-scale standard deviation)
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
plogis<-function(z,z50,sd,fsz=0){
    #cat(z,'\n')
    #cat('z50, sd = ',z50,sd,'\n')
    res<-1.0/(1.0+exp(-(z-z50)/sd));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-(fsz-z50)/sd));
    } else if (fsz<0){
        scl<-1/max(res);
    }
    res<-scl*res;
    #print(res);
    names(res)<-as.character(z);
    #print(res)
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and slope
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and slope
#'
#'@param z     - vector of sizes at which to compute selectivities
#'@param z50   - size at which selectivity  = 0.5 (logit-scale mean)
#'@param slope - slope at z50
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
asclogistic<-function(z,z50,slope,fsz=0){
    #cat(z,'\n')
    #cat('z50, lnD = ',z50,lnD,'\n')
    res <- 1.0/(1.0+exp(-slope*(z-z50)));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-slope*(fsz-z50)));
    } else if (fsz<0){
        scl<-1.0/max(res);
    }
    res<-scl*res;
    #print(res);
    names(res)<-as.character(z);
    #print(res)
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and z95
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and z95
#'
#'@param z    - vector of sizes at which to compute selectivities
#'@param z50 - size at which selectivity  = 0.5 (logit-scale mean)
#'@param z95 - size at which selectivity  = 0.95
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
asclogistic5095<-function(z,z50,z95,fsz=0){
    slope<-log(19.0)/(z95-z50);
    return(asclogistic(z,z50,slope,fsz=fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and D95 (=z95-z50)
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and D95
#'
#'@param z    - vector of sizes at which to compute selectivities
#'@param z50 - size at which selectivity  = 0.5 (logit-scale mean)
#'@param D95 - z95-z50
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
asclogistic50D95<-function(z,z50,D95,fsz=0){
    slope<-log(19.0)/D50;
    return(asclogistic(z,z50,slope,fsz=fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and ln(z95-z50)
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and ln(z95-z50)
#'
#'@param z    - vector of sizes at which to compute selectivities
#'@param z50 - size at which selectivity  = 0.5 (logit-scale mean)
#'@param lnD95 - ln-scale difference beteen z50 and z95
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
asclogistic50LnD95<-function(z,z50,lnD95,fsz=0){
    slope<-log(19.0)/exp(lnD95);
    return(asclogistic(z,z50,slope,fsz=fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by ln(z50) and ln(z95-z50)
#'
#'@description Function to calculate an ascending logistic function parameterized by ln(z50) and ln(z95-z50)
#'
#'@param z    - vector of sizes at which to compute selectivities
#'@param lnZ50 - ln-scale size at which selectivity  = 0.5 (logit-scale mean)
#'@param lnD95 - ln-scale difference beteen z50 and z95
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
asclogisticLn50LnD95<-function(z,lnZ50,lnD95,fsz=0){
    z50<-exp(lnZ50);
    slope<-log(19.0)/exp(lnD95);
    return(asclogistic(z,z50,slope,fsz=fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic function parameterized by z50 and slope for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by z50 and slope for ascending/descending limbs
#'
#'@param z     - vector of sizes at which to compute selectivities
#'@param ascZ50   - ascending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param ascSlope - ascending logistic slope at z50
#'@param dscZ50   - descending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param dscSlope - descending logistic slope at z50
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
dbllogistic<-function(z,ascZ50,ascSlope,dscZ50,dscSlope,fsz=0){
    #cat(z,'\n')
    #cat('z50, lnD = ',z50,lnD,'\n')
    res <- 1.0/(1.0+exp(-ascSlope*(z-ascZ50)))*1.0/(1.0+exp(-dscSlope*(z-dscZ50)));
    scl <-1;
    if (fsz>0){
        scl<-(1.0+exp(-ascSlope*(fsz-ascZ50)))*(1.0+exp(-dscSlope*(fsz-dscZ50)));
    } else if (fsz<0){
        scl<-1.0/max(res);
    }
    res<-scl*res;
    #print(res);
    names(res)<-as.character(z);
    #print(res)
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic function parameterized by z50 and slope for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by z50 and slope for ascending/descending limbs
#'
#'@param z     - vector of sizes at which to compute selectivities
#'@param ascZ50   - ascending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param ascSlope - ascending logistic slope at z50
#'@param lnD50    - ln-scale difference between ascending, descending limb z50s
#'@param dscSlope - descending logistic slope at z50
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
dbllogisticLnD50<-function(z,ascZ50,ascSlope,lnD50,dscSlope,fsz=0){
    dscZ50<-ascZ50+exp(lnD50);
    return(dbllogistic(ascZ50,ascSlope,dscZ50,dscSlope,fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic function parameterized by z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@param z      - vector of sizes at which to compute selectivities
#'@param ascZ50 - ascending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param ascD95 - increment from z50 to z95 on ascending limb
#'@param dscZ50 - descending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param dscD95 - increment from z50 to z95 on descending limb
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
dbllogistic50D95<-function(z,ascZ50,ascD95,dscZ50,dscD95,fsz=0){
    ascSlope<-log(19.0)/ascD95;
    dscSlope<-log(19.0)/dscD95;
    return(dbllogistic(ascZ50,ascSlope,dscZ50,dscSlope,fsz));
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate a double logistic function parameterized by ln-scale versions of z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@description Function to calculate a double logistic function parameterized by ln-scale z50 and D95=z95-z50 for ascending/descending limbs
#'
#'@param z      - vector of sizes at which to compute selectivities
#'@param ascLnZ50 - ln-scale ascending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param ascLnD95 - ln-scale increment from z50 to z95 on ascending limb
#'@param dscLnZ50 - ln-scale descending logistic size at which selectivity  = 0.5 (logit-scale mean)
#'@param dscLnD95 - ln-scale increment from z50 to z95 on descending limb
#'@param fsz   - if fsz>0, fsz=fully-selected size. if fsz<0, function is normalized to max. if fsz=0, no re-scaling is done
#'
#'@return vector with selectivity values at the elements of z
#'
dbllogisticLn50LnD95<-function(z,ascLnZ50,ascLnD95,dscLnZ50,dscLnD95,fsz=0){
    ascZ50<-exp(ascLnZ50);
    ascD95<-exp(ascLnD95);
    dscZ50<-exp(dscLnZ50);
    dscD95<-exp(dscLnD95);
    return(dbllogistic5095(ascZ50,ascD95,dscZ50,dscD95,fsz));
}
