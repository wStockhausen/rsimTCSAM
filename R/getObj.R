#'
#'@title Get an object from a list object using the specified path into the list.
#'
#'@description Function to get an object from a list object using the specified path into the list.
#'
#'@param path - path in the list to the desired object (using '/' as separator for list levels)
#'@param lst - the list to extract the object from
#'@param verbose - flag (T/F) to print diagnostics
#'
#'@return the specified object, or NULL if not found
#'
#'@export
#'
getObj<-function(path,lst,verbose=FALSE){
    if (verbose) cat("Getting object at '",path,"'\n",sep='')
    lvls<-strsplit(path,"/",fixed=TRUE);
    obj<-tryCatch(lst[[lvls[[1]]]],error=function(e){NULL;});
    return(obj);
}

