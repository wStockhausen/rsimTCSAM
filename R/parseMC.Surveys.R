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
#'@import wtsUtilities
#'
#'@export
#'
parseMC.Surveys<-function(rsp,i,dims){
    n <- dims$f$n;
    mny <- dims$y$mny;
    mxy <- dims$y$mxy;
    asy <- dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Surveys');
    surveys<-list();
    for (vp in 1:n){
        v<-rsp[[i]][1]; i<-i+1;#survey name
#         flags<-as.logical(rsp[[i]][1:3]); i<-i+1;
#         error<-as.numeric(rsp[[i]][1:3]); i<-i+1;
        resOI<-parseMC.OutputInfo(rsp,i); i<-resOI$i;
        blocks<-list();
        nt<-parseNum(rsp[[i]][1]); i<-i+1;
        for (tp in 1:nt){
            t<-rsp[[i]][1]; i <- i+1;
            eval(parse(text=paste('years<-',t)));
            lnQ  <-parseNum(rsp[[i]][1]);
            sdQ  <-parseNum(rsp[[i]][2]);
            lnQX <-parseNum(rsp[[i]][3]);
            i<-i+1;
            sel<-list();
            nc<-parseNum(rsp[[i]][1]); i<-i+1;
            for (ic in 1:nc){
                x <-rsp[[i]][1];               #sex
                ct<-rsp[[i]][2];               #curve type
                ft<-rsp[[i]][3];               #function type
                np<-parseNum(rsp[[i]][4]);     #number of parameters for function
                ps<-parseNum(rsp[[i]][4+1:np]);#parameter values
                if (ct=='selectivity'){
                    sel[[x]]<-list(type=ft,params=ps);
                } else {
                    cat('Unrecognized curve type "',ct,'"\n');
                    cat('Should be "selectivity"\n');
                    cat('Aborting...');
                    stop();
                }
                i<-i+1;
            }#ic
            block<-list(years=years,
                        lnQ=lnQ,sdQ=sdQ,lnQX=lnQX,
                        sel=sel);
            blocks[[t]]<-block;
        }#t
        surveys[[v]]<-list(name=v,output=resOI$output,blocks=blocks);
    }#vp
    return(list(i=i,surveys=surveys));
}