#'
#'@title Parse surveys section of a model configuration file.
#'
#'@description Function to parse the 'surveys' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of surveys section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'surveys - list object with surveys info
#'
#'@export
#'
parseMC.Surveys<-function(rsp,i,dims){
    n <- dims$v$n;
    mny <- dims$y$mny;
    mxy <- dims$y$mxy;
    asy <- dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Surveys');
    surveys<-list();
    for (vp in 1:n){
        v<-rsp[[i]][1]; i<-i+1;#survey name
        resOI<-parseMC.OutputInfo(rsp,i); i<-resOI$i;
        blocks<-list();
        nt<-wtsUtilities::parseNum(rsp[[i]][1]); i<-i+1;
        for (tp in 1:nt){
            t<-rsp[[i]][1]; i <- i+1;
            eval(parse(text=paste('years<-',t)));
            params<-list();
            for (ix in 1:dims$x$n){
                x    <-rsp[[i]][1];           #sex
                lnQ  <-wtsUtilities::parseNum(rsp[[i]][2]); #log(Q)
                sdQ  <-wtsUtilities::parseNum(rsp[[i]][3]); #sd(Q)
                idSel<-wtsUtilities::parseNum(rsp[[i]][4]); #selectivity function index
                params[[x]]<-list(lnQ=lnQ,sdQ=sdQ,idSel=idSel);
                i<-i+1;
            }#ix
            blocks[[t]]<-list(years=years,params=params);
        }#t
        surveys[[v]]<-list(name=v,output=resOI$output,blocks=blocks);
    }#vp
    return(list(i=i,surveys=surveys));
}