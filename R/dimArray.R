#'
#'@title Get an array
#'
#'@param mc - model configurqtion object
#'@param str - string describing array indices
#'
#'@details To match arrays in TCSAM2015, indices reflecting parameter combinations 
#'(such as for selectivity functions) should be indicated as "pc_xxx", 
#'where xxx is the rsim dimension name (e.g., 'selfcns'). The resulting array will use
#''pc' as the corresponding dimension name
#'
#'@return array of appropriate dimensions filled with default value(s)
#'
#'@export
#'
dimArray<-function(mc,str,val=0){
    d<-mc$dims;
    strp<-strsplit(str,'.',fixed=TRUE);
    ns<-length(strp[[1]]);#number of dimensions
    dms<-0*(1:ns);        #dimensions sizes
    dmnms<-list();        #list of dimensions (names) with index values
    #cat("indices = ")
    for (i in 1:ns){
        idx <- strp[[1]][i];#index name
        cdx <- strsplit(idx,"_")[[1]];#check if this has form "pc_xxx"
        if (cdx[1]=='pc'){
            #cat("'",cdx[2],"' ",sep='')
            dms[i]      <-d[[cdx[2]]]$n;
            dmnms[['pc']]<-d[[cdx[2]]]$nms;
        } else {
            #cat("'",idx,"' ",sep='')
            dms[i]      <-d[[idx]]$n;
            dmnms[[idx]]<-d[[idx]]$nms;
        }
    }
    A<-array(data=val,dim=dms,dimnames=dmnms);
    return(A);
}