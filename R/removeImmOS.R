#'
#'@title Remove rows marked as immature, old shell from a dataframe
#'
#'@description Function to remove rows marked as immature, old shell from a dataframe.
#'
#'@param mdfr - datarame from which to remove the immature, old shell factor combination
#'@param maturity - name of column specifying maturity state
#'@param shell_condition - name of column specifying shell condition
#'
#'@return dataframe with immature, old shell factor combination removed
#'
#'@export
#'
removeImmOS<-function(mdfr,maturity='m',shell_condition='s'){
    idx<-!((mdfr[[maturity]]=='immature')&(mdfr[[shell_condition]]=='old shell'))
    mdfr<-mdfr[idx,];
    return(mdfr);
}