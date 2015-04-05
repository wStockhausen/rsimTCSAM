#'
#'@title Parse selectivity section of a model configuration file.
#'
#'@description Function to parse the 'selectivity' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of selectivity section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'surveys - list object with selectivity info
#'
#'@import wtsUtilities
#'
#'@export
#'
parseMC.SelFcns<-function(rsp,i,dims){
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Selectivity');
    
    selfcns<-list();
    ns<-parseNum(rsp[[i]][1]); i<-i+1;
    for (s in 1:ns){
        idx<-parseNum(rsp[[i]][1]);
        ft<-rsp[[i]][2];               #function type
        np<-parseNum(rsp[[i]][3]);     #number of parameters for function
        ps<-parseNum(rsp[[i]][3+1:np]);#parameter values
        selfcns[[s]]<-list(type=ft,params=ps);
        i<-i+1;
    }
    
    return(list(i=i,selfcns=selfcns));
}
