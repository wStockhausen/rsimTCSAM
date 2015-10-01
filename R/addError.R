#'
#'@title Add error to "observations.
#'
#'@description Function to add error to "observations".
#'
#'@param x - value(s) to add errors to
#'@param cv - error cv
#'@param sdv - error standard deviation
#'@param ss - sample size for multinomial distributions
#'@param type - error distribution (NORMAL, LOGNORMAL, MULTINOMIAL)
#'
#'@return number (or, if x is a vector, vector) with error(s) added
#'
#'@export
#'
addError<-function(x,cv=NULL,ss=NULL,sdv=NULL,type='NORMAL'){
    y<-x;
    if (toupper(type)=='NORMAL'){
        if (!is.null(cv)) sdv <- x*cv;
        y <- rnorm(length(x),mean=x,sd=sdv);
    } else
    if (toupper(type)=='LOGNORMAL'){
        if (!is.null(cv)) sdv <- sqrt(log(1+cv^2));
        y <- rlnorm(length(x),meanlog=log(x),sd=sdv);
    } else
    if (toupper(type)=='MULTINOMIAL'){
        xs<-sum(x,na.rm=TRUE);
        xp<-x/xs;
        y <- (xs/ss)*as.vector(rmultinom(1,ss,xp));
    }
    return(y);
}