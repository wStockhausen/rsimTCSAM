#'
#'@title Parse model dimensions information from a character string vector.
#'
#'@param rz - character string vector from which to extract model dimensions info
#'@param i - index into character string vector at which to start parsing
#'
#'@return a list with elements
#'dims : model dimensions list object
#'i    : final index + 1 into rz at end of parsing dims info
#'
#'@export
#'
parseMC.Dims<-function(rz,i){
    
    #create dims list object
    dims<-list(y=list(n=0,nms=NULL,vls=NULL,mny=0,mxy=0,asy=0),
               x=list(n=0,nms=NULL),
               m=list(n=0,nms=NULL),
               s=list(n=0,nms=NULL),
               z=list(n=0,nms=NULL,vls=NULL),
               zp=list(n=0,nms=NULL,vls=NULL),
               zc=list(n=0,nms=NULL,vls=NULL),
               f=list(n=0,nms=NULL),
               v=list(n=0,nms=NULL),
               selfcns=list(n=0,nms=NULL,lbls=NULL)
               );
    
    #set text row counter
    j<-0;
    
    #DIMESIONS keyword
    chk<-rz[[i+j]][1]; j<-j+1;
    checkKeyword(chk,'DIMENSIONS');
    
    #years
    mny<-parseNum(rz[[i+j]][1]); j<-j+1;#start model year
    asy<-parseNum(rz[[i+j]][1]); j<-j+1;#assessment year (=mxy+1)
    mxy<-asy-1;
    cat('model years = ',mny,":",mxy,'. Assessment year = ',asy,'\n',sep='')
    dims$y$mny<-mny;
    dims$y$mxy<-mxy
    dims$y$asy<-asy;
    dims$y$n<-asy-mny+1;
    dims$y$vls<-mny:asy
    dims$y$nms<-as.character(mny:asy);
    
    #size bins and cutpoints
    mnZC<-parseNum(rz[[i+j]][1]); j<-j+1;
    mxZC<-parseNum(rz[[i+j]][1]); j<-j+1;
    delZ<-parseNum(rz[[i+j]][1]); j<-j+1;
    dims$zc$vls<-seq(from=mnZC,to=mxZC,by=delZ);
    dims$zc$nms<-as.character(dims$zc$vls);
    dims$zc$n<-length(dims$zc$vls);
    dims$z$vls<-dims$zc$vls[(2:dims$zc$n)-1]+0.5*first_difference(dims$zc$vls);
    dims$z$nms<-as.character(dims$z$vls);
    dims$z$n<-length(dims$z$vls);
    dims$zp$vls<-dims$z$vls
    dims$zp$nms<-dims$z$nms;
    dims$zp$n<-dims$z$n;
    cat('size bins =',dims$z$nms,'\n');
    
    #sexes
    dims$x$n<-parseNum(rz[[i+j]][1]); j<-j+1;
    dims$x$nms<-rz[[i+j]][1:dims$x$n]; j<-j+1;
    cat("sexes = ",addQuotes(dims$x$nms),'\n')
    
    #maturity states
    dims$m$n<-parseNum(rz[[i+j]][1]); j<-j+1;
    dims$m$nms<-rz[[i+j]][1:dims$m$n]; j<-j+1;
    cat("maturity states = ",addQuotes(dims$m$nms),'\n')
    
    #shell condition
    dims$s$n<-parseNum(rz[[i+j]][1]); j<-j+1;
    dims$s$nms<-rz[[i+j]][1:dims$s$n]; j<-j+1;
    cat("shell condition = ",addQuotes(dims$s$nms),'\n')
    
    #fisheries
    dims$f$n<-parseNum(rz[[i+j]][1]); j<-j+1;
    dims$f$nms<-rz[[i+j]][1:dims$f$n]; j<-j+1;
    cat("fisheries = ",addQuotes(dims$f$nms),'\n')
    
    #surveys
    dims$v$n<-parseNum(rz[[i+j]][1]); j<-j+1;
    dims$v$nms<-rz[[i+j]][1:dims$v$n]; j<-j+1;
    cat("surveys = ",addQuotes(dims$v$nms),'\n')
    
    #selectivity/retention functions
    dims$selfcns$n<-parseNum(rz[[i+j]][1]); j<-j+1;
    dims$selfcns$nms<-as.character(1:dims$selfcns$n);
    dims$selfcns$lbls<-rz[[i+j]][1:dims$selfcns$n]; j<-j+1;
    cat("selfcns = ",addQuotes(dims$selfcns$lbls),'\n')
    
    cat("--finished reading dims list object\n\n")
    
    return(list(dims=dims,i=i+j));
}
