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
        retOI<-parseMC.OutputInfo(rsp,i); i<-retOI$i;#retained catch output info
        dscOI<-parseMC.OutputInfo(rsp,i); i<-dscOI$i;#discard catch output info
        totOI<-parseMC.OutputInfo(rsp,i); i<-totOI$i;#total catch output info
        blocks<-list();
        nt<-wtsUtilities::parseNum(rsp[[i]][1]); i<-i+1;
        for (tp in 1:nt){
            t<-rsp[[i]][1]; i<-i+1;
            eval(parse(text=paste('years<-',t)));
            params<-list();
            for (ix in 1:dims$x$n){
                x    <-rsp[[i]][1];          #sex
                hm   <-wtsUtilities::parseNum(rsp[[i]][2]);#handling mortality
                lnF  <-wtsUtilities::parseNum(rsp[[i]][3]);#log(F)
                sdF  <-wtsUtilities::parseNum(rsp[[i]][4]);#sd(F)
                idSel<-wtsUtilities::parseNum(rsp[[i]][5]);#selectivityfunction index
                idRet<-wtsUtilities::parseNum(rsp[[i]][6]);#retention function index
                params[[x]]<-list(hm=hm,lnF=lnF,sdF=sdF,idSel=idSel,idRet=idRet);
                i<-i+1;
            }#ix
            blocks[[t]]<-list(years=years,params=params);
        }#tp
        fisheries[[f]]<-list(name=f,
                             output=list(ret=retOI$output,
                                         dsc=dscOI$output,
                                         tot=totOI$output),
                             blocks=blocks);
    }#fp
    
    return(list(i=i,fisheries=fisheries));
}