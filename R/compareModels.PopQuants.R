#'
#'@title Compare population quantities from TCSAM2015 and rsimTCSAM model runs.
#'
#'@description Function to compare population quantities from TCSAM2015 and rsimTCSAM model runs.
#'
#'@param tcsams - single TCSAM2015 model report object, or named list of such
#'@param rsims - single rsimTCSAM results object, or named list of such
#'@param showPlot - flag to show/print plots immediately
#'@param pdf - name of pdf file to record plot output to
#'@param width - pdf page width (in inches)
#'@param height - pdf page width (in inches)
#'@param verbose - flag (T/F) to print debug info
#'
#'@return list of ggplot2 objects
#'
#'@details none.
#'
#'@export
#'
compareModels.PopQuants<-function(tcsams=NULL,
                                      rsims=NULL,
                                      showPlot=TRUE,
                                      pdf=NULL,
                                      width=8,
                                      height=6,
                                      verbose=FALSE){
    #set up pdf device, if requested
    if (!is.null(pdf)){
        pdf(file=pdf,width=width,height=height);
        on.exit(dev.close())
    }
    
    if (inherits(tcsams,'tcsam2015.rep')){
        tcsams<-list(tcsam=tcsams);#wrap in list
    }
    if (inherits(rsims,'rsimTCSAM')){
        rsims<-list(rsim=rsims);#wrap in list
    }
    
    plots<-list();
    
    #abundance trends
    if (verbose) cat("Plotting population abundance trends\n");
    mdfr<-getMDFR('mr/P_list/N_yxmsz',tcsams,rsims);
    dfr<-reshape2::dcast(mdfr,modeltype+model+y+x+m+s~.,fun.aggregate=sum,value.var='val');
    p1<-plotMDFR.XY(dfr,x='y',value.var='.',faceting='x~m',
                   plotABline=TRUE,plotPoints=FALSE,
                   xlab='year',ylab='Abundance (millions)',units="",
                   linetype='s',guideTitleLineType='',
                   colour='model',guideTitleColour='');
    if (showPlot||!is.null(pdf)) print(p1);
    p2<-plotMDFR.XY(dfr[dfr$y>=1980,],x='y',value.var='.',faceting='x~m',
                   plotABline=TRUE,plotPoints=FALSE,
                   xlab='year',ylab='Abundance (millions)',units="",
                   linetype='s',guideTitleLineType='',
                   colour='model',guideTitleColour='');
    if (showPlot||!is.null(pdf)) print(p2);
    plots$N_yxms<-list(p1,p2);
    
    #biomass trends
    if (verbose) cat("Plotting population biomass trends\n");
    mdfr<-getMDFR('mr/P_list/B_yxms',tcsams,rsims);
    p1<-plotMDFR.XY(mdfr,x='y',value.var='val',faceting='x~m',
                   plotABline=TRUE,plotPoints=FALSE,
                   xlab='year',ylab='Biomass (1000s t)',units="",
                   linetype='s',guideTitleLineType='',
                   colour='model',guideTitleColour='');
    if (showPlot||!is.null(pdf)) print(p1);
    p2<-plotMDFR.XY(mdfr[mdfr$y>=1980,],x='y',value.var='val',faceting='x~m',
                   plotABline=TRUE,plotPoints=FALSE,
                   xlab='year',ylab='Biomass (1000s t)',units="",
                   linetype='s',guideTitleLineType='',
                   colour='model',guideTitleColour='');
    if (showPlot||!is.null(pdf)) print(p2);
    plots$B_yxms<-list(p1,p2);
    
    #mature biomass at mating trends
    if (verbose) cat("Plotting population mature biomass-at-mating trends\n");
    mdfr<-getMDFR('mr/P_list/MB_yx',tcsams,rsims);
    p1<-plotMDFR.XY(mdfr,x='y',value.var='val',faceting='x~.',
                   plotABline=TRUE,plotPoints=FALSE,
                   xlab='year',ylab='Mating Biomass (1000s t)',units="",
                   colour='model',guideTitleColour='');
    if (showPlot||!is.null(pdf)) print(p1);
    p2<-plotMDFR.XY(mdfr[mdfr$y>=1980,],x='y',value.var='val',faceting='x~.',
                   plotABline=TRUE,plotPoints=FALSE,
                   xlab='year',ylab='Mating Biomass (1000s t)',units="",
                   colour='model',guideTitleColour='');
    if (showPlot||!is.null(pdf)) print(p2);
    plots$MB_yx<-list(p1,p2);
    

    #recruitment
    if (verbose) cat("Plotting recruitment time series\n");
    path<-'mp/R_list/R_y';
    mdfr<-getMDFR(path,tcsams,rsims);
    p<-plotMDFR.XY(mdfr,x='y',agg.formula=NULL,faceting=NULL,
                   xlab='year',ylab='Recruitment',units='millions',lnscale=FALSE,
                   colour='model',guideTitleColor='',
                   shape='model',guideTitleShape='');
    if (showPlot||!is.null(pdf)) print(p);
    plots$R_y<-p;
    p<-plotMDFR.XY(mdfr,x='y',agg.formula=NULL,faceting=NULL,
                   xlab='year',ylab='Recruitment',units='millions',lnscale=TRUE,
                   colour='model',guideTitleColor='',
                   shape='model',guideTitleShape='');
    if (showPlot||!is.null(pdf)) print(p);
    plots$lnR_y<-p;
    
    #Population abundance-at-size
    if (verbose) cat("Plotting poulation abundance-at-size\n");
    path<-'mr/P_list/N_yxmsz';
    mdfr<-getMDFR(path,tcsams,rsims);
    mdfr<-removeImmOS(mdfr);
    p<-plotMDFR.Bubbles(mdfr,x='y',y='z',
                        agg.formula='model+y+x+z',faceting='model~x',
                        xlab='year',ylab='size (mm CW)',units="millions",
                        colour='.',guideTitleColour='',useColourGradient=TRUE,alpha=0.5);
    if (showPlot||!is.null(pdf)) print(p);
    plots$N_yxmsz<-p;
        
    #Population abundance-at-size (bubble plots)
    if (verbose) cat("Plotting poulation abundance-at-size (bubble plots)\n");
    path<-'mr/P_list/N_yxmsz';
    mdfr<-getMDFR(path,tcsams,rsims);
    mdfr<-removeImmOS(mdfr);
    p<-plotMDFR.Bubbles(mdfr,x='y',y='z',
                        agg.formula='model+y+x+z',faceting='model~x',
                        xlab='year',ylab='size (mm CW)',units="millions",
                        colour='.',guideTitleColour='',useColourGradient=TRUE,alpha=0.5);
    if (showPlot||!is.null(pdf)) print(p);
    plots$N_yxmsz<-p;
        
    #Population abundance-at-size (line plots)
    cat("TODO: In compareModels.PopQuants(...), implement plots for population abundance-at-size as line plots.\n");
    #p<-compareModels.SizeComps(mdfr);
    
    return(invisible(plots))
}