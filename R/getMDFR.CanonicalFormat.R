#'
#'@title Convert dataframe to canonical format
#'
#'@description Function to convert dataframe to canonical format.
#'
#'@param mdfr - dataframe to convert
#'
#'@return dataframe in canonical format
#'
#'@details returns a dataframe in canonical format
#'
#'@export
#'
getMDFR.CanonicalFormat<-function(mdfr){
    #check existing columns and add missing ones
    nms<-names(mdfr);
    if (!('case' %in% nms))     mdfr[['case']]    <-"";
    if (!('type' %in% nms))     mdfr[['type']]    <-"";
    if (!('category' %in% nms)) mdfr[['category']]<-"";
    if (!('fleet' %in% nms))    mdfr[['fleet']]   <-"";
    if (!('pc' %in% nms))       mdfr[['pc']]      <-"";
    if (!('y' %in% nms)) mdfr[['y']]<-"all";
    if (!('x' %in% nms)) mdfr[['x']]<-"all";
    if (!('m' %in% nms)) mdfr[['m']]<-"all";
    if (!('s' %in% nms)) mdfr[['s']]<-"all";
    if (!('z' %in% nms)) mdfr[['z']]<-"all";
    if (!('val' %in% nms)) mdfr[['val']]<-NA;
    if (!('zp' %in% nms)) {
        if (!('lci' %in% nms)) mdfr[['lci']]<-NA;
        if (!('uci' %in% nms)) mdfr[['uci']]<-NA;
    }
    
    #check for fisheries and surveys info
    if ('f' %in% nms) {
        mdfr$type<-'fishery';
        mdfr$fleet<-mdfr$f;
    } else if ('v' %in% nms) {
        mdfr$type<-'survey';
        mdfr$fleet<-mdfr$v;
    } else {
        mdfr$type<-'population';
        mdfr$fleet<-'';
    }
    
    #re-order to canconical format
    if ('zp' %in% nms){
        mdfr<-mdfr[,c('case','type','fleet','category','pc','y','x','m','s','z','zp','val')]
    } else {
        mdfr<-mdfr[,c('case','type','fleet','category','pc','y','x','m','s','z','val','lci','uci')]
    }

    return(mdfr);
}
