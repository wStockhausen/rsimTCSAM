#'
#'@title Parse fisheries section of a model configuration file.
#'
#'@description Function to parse the 'fisheries' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of fisheries section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'fisheries - list object with fisheries info
#'
#'@import wtsUtilities
#'
#'@export
#'
parseMC.Fisheries<-function(rsp,i,dims){
    n <- dims$f$n;
    mny <- dims$y$mny;
    mxy <- dims$y$mxy;
    asy <- dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Fisheries');
    fisheries<-list();
    for (fp in 1:n){
        f<-rsp[[i]][1]; i<-i+1;
#         retFlags<-as.logical(rsp[[i]][1:3]); i<-i+1;
#         dscFlags<-as.logical(rsp[[i]][1:3]); i<-i+1;
#         totFlags<-as.logical(rsp[[i]][1:3]); i<-i+1;
#         retErr<-as.numeric(rsp[[i]][1:3]); i<-i+1;
#         dscErr<-as.numeric(rsp[[i]][1:3]); i<-i+1;
#         totErr<-as.numeric(rsp[[i]][1:3]); i<-i+1;
        retOI<-parseMC.OutputInfo(rsp,i); i<-retOI$i;#retained catch output info
        dscOI<-parseMC.OutputInfo(rsp,i); i<-dscOI$i;#discard catch output info
        totOI<-parseMC.OutputInfo(rsp,i); i<-totOI$i;#total catch output info
        blocks<-list();
        nt<-parseNum(rsp[[i]][1]); i<-i+1;
        for (tp in 1:nt){
            t<-rsp[[i]][1]; i<-i+1;
            eval(parse(text=paste('years<-',t)));
            hm   <-parseNum(rsp[[i]][1]);
            mnF  <-parseNum(rsp[[i]][2]);
            sdF  <-parseNum(rsp[[i]][3]);
            offFX<-parseNum(rsp[[i]][4]);
            i<-i+1;
            sel<-list();
            ret<-list();
            nc<-parseNum(rsp[[i]][1]); i<-i+1;
            for (ic in 1:nc){
                x <-rsp[[i]][1];               #sex
                ct<-rsp[[i]][2];               #curve type (sel or ret)
                ft<-rsp[[i]][3];               #function type
                np<-parseNum(rsp[[i]][4]);     #number of parameters for function
                ps<-parseNum(rsp[[i]][4+1:np]);#parameter values
                if (ct=='selectivity'){
                    sel[[x]]<-list(type=ft,params=ps);
                } else if (ct=='retention'){
                    ret[[x]]<-list(type=ft,params=ps);
                } else {
                    cat('Unrecognized curve type "',ct,'"\n');
                    cat('Should be "selectivity" or "retention" \n');
                    cat('Aborting...');
                    stop();
                }
                i<-i+1;
            }#ic
            block<-list(years=years,
                        hm=hm,mnF=mnF,sdF=sdF,offFX=offFX,
                        sel=sel,ret=ret);
            blocks[[t]]<-block;
        }#tp
        fisheries[[f]]<-list(name=f,
                             output=list(ret=retOI$output,
                                         dsc=dscOI$output,
                                         tot=totOI$output),
                             blocks=blocks);
    }#fp
    
    return(list(i=i,fisheries=fisheries));
}