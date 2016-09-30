#'
#'@title Create a rsimTCSAM.resLst object from a model run
#'
#'@description Function to create a rsimTCSAM.resLst object from a model run.
#'
#'@param path - path to .RData object containing model run
#'
#'@return a rsimTCSAM.resLst object.
#'
#'@details Uses \code{tcltk::tk_chose.file} to open a file dialog to select model .RData file
#'if path is NULL. A rsimTCSAM.resLst object is a list with elements
#'\itemize{
#'  \item{rep - a rsimTCSAM object, or NULL}
#'}
#'
#'@export
#'
getResLst<-function(path=NULL){
    if (is.null(path)){
        path<-tcltk::tk_choose.files(default="rsim.RData",
                                     caption="Select file with model results",
                                     filters=matrix(data=c("*.RData",".RData"),nrow=1,ncol=2));
        if (length(path)==0) {
            cat("User canceled file selection!! Returning NULL as model results.\n")
            return(NULL);#user canceled file selection
        }
    }

    if (!file.exists(path)) {
        cat("Warning from getResLst(...).\n");
        cat("--The following file does not exist:\n\t'",path,"'\n",sep='');
        cat("--Returning NULL.\n")
        return(NULL);
    }


    cat("Reading model report from file:\n",path,"\n")
    source(path,local=TRUE);
    if(!any(names(res)=='mc')){
            cat("The file '",path,"'\n",
                "\tdoes not appear to be a rsimTCSAM model results file.\n",
                "\trsimTCSAM results files are R lists, with 'mc' as the first element.\n",
                "\tReturning NULL.\n",sep="");
            return(NULL);
    }
    class(res)<-c('rsimTCSAM',class(res));#set class attribute to 'rsimTCSAM' for identification
    
    resLst<-list(rep=res);
    class(resLst)<-c('rsimTCSAM.resLst',class(resLst));

    return(resLst);
}
