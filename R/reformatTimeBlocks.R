#'
#'@title Reformat rsim, tcsam notation for time blocks to a standard format
#'
#'@description Reformat rsim, tcsam notation for time blocks to a standard format.
#'
#'@param tb0 - time block string to reformat
#'@param dims - model 'dims' list from tscam or rsim object
#'
#'@return time block string in "standardized" format
#'
#'@details For 2 time periods and a single year, the standardized format 
#'looks like 'yyy1-yyy2;yyy3-yyy4;yyy5'.
#'
#'@export
#'
reformatTimeBlocks<-function(tb0,dims){
    mny<-dims$y$mny;
    mxy<-dims$y$mxy;
    asy<-dims$y$asy;
    tb1<-gsub("c(",  "",tb0,fixed=TRUE);#remove "c("s
    tb2<-gsub( ")",  "",tb1,fixed=TRUE);#remove ")"s
    tb3<-gsub( ",", ";",tb2,fixed=TRUE);#replace commas with semi-colons
    tb4<-gsub( "[",  "",tb3,fixed=TRUE);#remove "["s
    tb5<-gsub( "]",  "",tb4,fixed=TRUE);#remove "]"s

    tb6<-gsub("mny",mny,tb5,fixed=TRUE);#replace 'mny' with value
    tb7<-gsub("mxy",mxy,tb6,fixed=TRUE);#replace 'mxy' with value
    tb8<-gsub("asy",asy,tb7,fixed=TRUE);#replace 'asy' with value
    
    tb9<-gsub( ":",  "-",tb8,fixed=TRUE);#replace ":" with "-"
    return(tb9);
}