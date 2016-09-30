#'
#'@title Parse growth section of a model configuration file.
#'
#'@description Function to parse the 'growth' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of growth section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'params - list object with growth time blocks
#'
#'@export
#'
parseMC.Growth<-function(rsp,i,dims){
    #get y dims for parsing time blocks
    mny<-dims$y$mny;
    mxy<-dims$y$mxy;
    asy<-dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Growth');
    blocks<-list();
    nt<-wtsUtilities::parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        a<-dimArray(list(dims=dims),'x');
        b<-dimArray(list(dims=dims),'x');
        s<-dimArray(list(dims=dims),'x');
        for (xp in 1:dims$x$n){
            x<-rsp[[i]][1]; 
            a[x]<-exp(wtsUtilities::parseNum(rsp[[i]][2])); 
            b[x]<-exp(wtsUtilities::parseNum(rsp[[i]][3])); 
            s[x]<-exp(wtsUtilities::parseNum(rsp[[i]][4])); 
            i<-i+1;
        }#xp
        blocks[[t]]<-list(years=years,
                          a_x=a,
                          b_x=b,
                          s_x=s
                         );
    }#blocks
    return(list(i=i,params=list(blocks=blocks)))
}
