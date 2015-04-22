#'
#'@title Parse weight-at-size section of a model configuration file.
#'
#'@description Function to parse the 'weight-at-size' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of weight-at-size section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'params - list object with weight-at-size time blocks
#'
#'@import wtsUtilities
#'
#'@export
#'
parseMC.WatZ<-function(rsp,i,dims){
    #get y dims for parsing time blocks
    mny<-dims$y$mny;
    mxy<-dims$y$mxy;
    asy<-dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'WatZ');
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        t<-rsp[[i]][1]; i<-i+1;
        eval(parse(text=paste('years<-',t)));
        a<-dimArray(list(dims=dims),'x.m');
        b<-dimArray(list(dims=dims),'x.m');
        for (xp in 1:dims$x$n){
            for (mp in 1:dims$m$n){
                x<-rsp[[i]][1]; 
                m<-rsp[[i]][2]; 
                av<-parseNum(rsp[[i]][3]); 
                bv<-parseNum(rsp[[i]][4]);
                a[x,m]<-av;
                b[x,m]<-bv;
                i<-i+1;
            }#mp
        }#xp
        blocks[[t]]<-list(years=years,
                          a_xm=a,
                          b_xm=b
                         );
    }#blocks
    return(list(i=i,params=list(blocks=blocks)))
}
