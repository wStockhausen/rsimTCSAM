#'
#'@title Throw an unrecognized model type error.
#'
#'@param imt - input model type
#'@param xmt - expected model types
#'@param str - string to print (name of calling function)
#'
#'@details Stops the R run.
#'
throwModelTypeError<-function(imt,xmt,str){
    cat('Model type not recognized in ',str,'.\n',sep='')
    cat("Input type was '",imt,"'.\n",sep='')
    cat("Expected types were ",paste("'",xmt,"'",sep='',collapse=", "),".\n",sep='')
    cat("Aborting...\n");
    stop();
}