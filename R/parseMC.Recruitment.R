#'
#'@title Parse recruitment section of a model configuration file.
#'
#'@description Function to parse the 'recruitment' section of a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of recruitment section in rsp
#'@param dims - model configuration dims list
#'
#'@return list with
#'i - index starting next section
#'params - list object with recruitment time blocks
#'
#'@import wtsUtilities
#'
#'@export
#'
parseMC.Recruitment<-function(rsp,i,dims){
    #get y dims for parsing time blocks
    mny<-dims$y$mny;
    mxy<-dims$y$mxy;
    asy<-dims$y$asy;
    
    chk<-rsp[[i]][1]; i<-i+1;
    checkKeyword(chk,'Recruitment');
    #read initialization section
    inits<-list();
    inits$lnR      <- parseNum(rsp[[i]][1]); #ln-scale mean recruitment
    inits$cvR      <- parseNum(rsp[[i]][2]); #ln-scale value for ln-scale recruitment standard deviation
    inits$lgtMnXR  <- parseNum(rsp[[i]][3]); #logit-scale nominal sex ratio
    inits$lgtSdXR  <- parseNum(rsp[[i]][4]); #logit-scale standard deviation for sex ratio deviations
    inits$lnAlphaZ <- parseNum(rsp[[i]][5]); #ln-scale alpha parameter for rec. size distribution
    inits$lnBetaZ  <- parseNum(rsp[[i]][6]); #ln-scale beta parameter for rec. size distribution
    i<-i+1;
    #read time blocks
    blocks<-list();
    nt<-parseNum(rsp[[i]][1]); i<-i+1;
    for (tp in 1:nt){
        block<-list();
        t<-rsp[[i]][1];
        eval(parse(text=paste('years<-',t)));
        block$years<-years;
        block$lnR      <- parseNum(rsp[[i]][2]); #ln-scale mean recruitment
        block$cvR      <- parseNum(rsp[[i]][3]); #ln-scale value for ln-scale recruitment standard deviation
        block$lgtMnXR  <- parseNum(rsp[[i]][4]); #logit-scale nominal sex ratio
        block$lgtSdXR  <- parseNum(rsp[[i]][5]); #logit-scale standard deviation for sex ratio deviations
        block$lnAlphaZ <- parseNum(rsp[[i]][6]); #ln-scale alpha parameter for rec. size distribution
        block$lnBetaZ  <- parseNum(rsp[[i]][7]); #ln-scale beta parameter for rec. size distribution
        blocks[[t]]<-block;
        i<-i+1;
    }#blocks
    return(list(i=i,params=list(inits=inits,blocks=blocks)))
}
