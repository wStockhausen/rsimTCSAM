#'
#'@title Parse molt-to-maturity section of a model configuration file.
#'
#'@description Function to parse the 'molt-to-maturity' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of molt-to-maturity section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'params - list object with molt-to-maturity time blocks
#'
#'@export
#'
parseMC.PrM2M<-function(rsp,i,dims){
    #get y dims for parsing time blocks
    mny<-dims$y$mny;
    mxy<-dims$y$mxy;
    asy<-dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'MoltToMaturity');
    blocks<-list();
    nt<-wtsUtilities::parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        z50_x<-dimArray(list(dims=dims),'x');
        sdv_x<-dimArray(list(dims=dims),'x');
        for (xp in 1:dims$x$n){
            x<-rsp[[i]][1]; 
            z50_x[x]<-wtsUtilities::parseNum(rsp[[i]][2]); 
            sdv_x[x]<-wtsUtilities::parseNum(rsp[[i]][3]); 
            i<-i+1;
        }#xp
        blocks[[t]]<-list(years=years,
                          z50_x=z50_x,
                          sdv_x=sdv_x
                          );
    }#blocks
    return(list(i=i,params=list(blocks=blocks)))
}
