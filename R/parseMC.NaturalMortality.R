#'
#'@title Parse natural mortality section of a model configuration file.
#'
#'@description Function to parse the 'natural mortality' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of natural mortality section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'params - list object with natural mortality time blocks
#'
#'@import wtsUtilities
#'
#'@export
#'
parseMC.NaturalMortality<-function(rsp,i,dims){
    cat('parsing Natural Mortality\n')
    #get y dims for parsing time blocks
    mny<-dims$y$mny;
    mxy<-dims$y$mxy;
    asy<-dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'NaturalMortality');
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        M0_xm  <- dimArray(list(dims=dims),'x.m');
        cvM_xm <- dimArray(list(dims=dims),'x.m');
        for (xp in 1:dims$x$n){
            for (mp in 1:dims$m$n){
                x<-rsp[[i]][1]; 
                m<-rsp[[i]][2]; 
                M0_xm[x,m]  <- parseNum(rsp[[i]][3]);;
                cvM_xm[x,m] <- parseNum(rsp[[i]][4]);;
                i<-i+1;
            }
        }
        blocks[[t]]<-list(years=years,
                          M0_xm=M0_xm,
                          cvM_xm=cvM_xm
                         );
    }#blocks
    cat('finished parsing Natural Mortality\n')
    return(list(i=i,params=list(blocks=blocks)))
}
