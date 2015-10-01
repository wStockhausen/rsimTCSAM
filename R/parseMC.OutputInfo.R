#'
#'@title Parse output info for a fishery catch type or survey from a model configuration file.
#'
#'@description Function to parse output info for a fishery catch type or survey from a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of output info section in rsp
#'
#'@return list with elements
#'i - index starting next section
#'output - list object with output info for abundance, biomass and sizecomps data types
#'
#'@export
#'
parseMC.OutputInfo<-function(rsp,i){
    output<-list();
    output$abundance<-parseMC.OILine(rsp,i); i<-i+1;
    output$biomass  <-parseMC.OILine(rsp,i); i<-i+1;
    output$sizecomps<-parseMC.OILine(rsp,i); i<-i+1;
    return(list(i=i,output=output));
}
#'
#'@title Parse output info for a fishery/survey data type from a model configuration file.
#'
#'@description Function to read output info for a fishery/survey data type from a model configuration file.
#'
#'@param rsp - parsed text list from model configuration file
#'@param i - index to start of output info section in rsp
#'
#'@details Note that i is not incremented here.
#'
#'@return list with elements
#'flag - flag (TRUE/FALSE) to output this data type
#'aggType - TCSAM aggregation type (e.g., BY_TOTAL, BY_SEX)
#'errType - TCSAM error structure type (e.g., LOGNORMAL, MULTINOMIAL)
#'err - output cv or sample size
#'addErr - flag to actually add error to output
#'wgt - likelihood weight
#'
#'@importFrom wtsUtilities parseNum
#'
#'@export
#'
parseMC.OILine<-function(rsp,i){
    lst<-list();
    lst$flag   <-as.logical(rsp[[i]][1]);
    lst$aggType<-rsp[[i]][2];
    lst$errType<-rsp[[i]][3];
    lst$err    <-wtsUtilities::parseNum(rsp[[i]][4]);
    lst$addErr <-as.logical(rsp[[i]][5]);
    lst$wgt    <-wtsUtilities::parseNum(rsp[[i]][6]);
    return(lst);
}
