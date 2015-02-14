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
    if (tolower(type)=='logistic'){
        res<-plogis(z,params[1],params[2]);
        res<-res/max(res);
    } else if (tolower(type)=='asclogistic50ln95'){
        res<-asclogistic50Ln95(z,params[1],params[2]);
        res<-res/max(res);
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
#'
#'@return vector with selectivity values at the elements of z
#'
plogis<-function(z,z50,sd){
    #cat(z,'\n')
    #cat('z50, sd = ',z50,sd,'\n')
    res<-1.0/(1.0+exp(-(z-z50)/sd));
    #print(res);
    names(res)<-as.character(z);
    #print(res)
    return(res)
}
#-----------------------------------------------------------------------------------
#'
#'@title Calculate an ascending logistic function parameterized by z50 and ln(z95-z50)
#'
#'@description Function to calculate an ascending logistic function parameterized by z50 and ln(z95-z50)
#'
#'@param z    - vector of sizes at which to compute selectivities
#'@param z50 - size at which selectivity  = 0.5 (logit-scale mean)
#'@param lnD - ln-scale difference beteen z50 and z95
#'
#'@return vector with selectivity values at the elements of z
#'
asclogistic50Ln95<-function(z,z50,lnD){
    #cat(z,'\n')
    #cat('z50, lnD = ',z50,lnD,'\n')
    res <- 1.0/(1.0+exp(-log(19.0)*(z-z50)/exp(lnD)));
    #print(res);
    names(res)<-as.character(z);
    #print(res)
    return(res)
}
