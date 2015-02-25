#'
#'@title Compare multiple size comps.
#'
#'@description Function to compare multiple size comps.
#'
#'@param n_xmsz - array (or list of arrays) dimensioned xmsz
#'@param title - title for plot
#'@param showPlot - flag to show plot immediately
#'
#'@return ggplot2 object
#'
#'@import ggplot2
#'@import reshape2
#'
#'@export
#'
compareSizeCompsGG<-function(n_xmsz=NULL,
                             title='',
                             ylab='Abundance (millions)',  
                             showPlot=TRUE){
    
    oneModel<-FALSE;
    #size comps come in as array(s)
    if (is.array(n_xmsz)){
        mdfr<-melt(n_xmsz,value.name='val');
        mdfr$model<-'';
        oneModel<-TRUE;
    } else if (is.list(n_xmsz)) {
        #n_xmsz is a list of array objects with models as names 
        mdls<-names(n_xmsz);
        mdfr<-NULL;
        for (mdl in mdls){
            mdfrp<-melt(n_xmsz[[mdl]],value.name='val');
            mdfrp$model<-mdl;
            mdfr<-rbind(mdfr,mdfrp);
        }
    } else {
        cat('Invalid specification for n_xmsz in compareSizeCompsGG\n');
        cat('n_xmsz must be either an array or a list of arrays\n');
        cat("Aborting...\n")
        stop();
    }
            
    if (oneModel){
        #plotting one model
        pl <- ggplot(aes(x=z,y=val,fill=s),data=mdfr);
        pl <- pl + facet_grid(m~x);
    } else {
        #plotting multiple models
        pl <- ggplot(aes(x=z,y=val,fill=model),data=mdfr);
        pl <- pl + facet_grid(m+s~x);
    }
    pl <- pl + geom_bar(alpha=1,stat='identity',position='dodge');
    pl <- pl + labs(x='size (mm)',y=ylab);
    pl <- pl + guides(fill=guide_legend(''));
    if (title!='') pl <- pl + ggtitle(title);
    if (showPlot) print(pl);
    
    return(invisible(pl))
}