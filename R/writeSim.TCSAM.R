#'
#'@title Write simulation results to TCSAM2015 input file.
#'
#'@description Function to write simulation results to TCSAM2015 input file.
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param mr - model results list object
#'@param fn - output file name
#'@param showPlot - flag to show plots
#'
#'@export
#'
writeSim.TCSAM<-function(mc,mp,mr,fn='TCAM2015.Data.dat',showPlot=TRUE){
    d<-mc$dims;
    mxy<-as.character(d$y$mxy);
    
    conn<-file(fn,open="w");
    on.exit(close(conn));
    
    cat("#####################################################################\n",file=conn);
    cat("#TCSAM2015 Model Configuration File                                 #\n",file=conn);
    cat("#####################################################################\n",file=conn);
    cat('rsimTest',"  # model scenario","\n",file=conn);
    cat(d$y$mny,"  # model start year (pop. model start year)","\n",file=conn);
    cat(d$y$mxy,"  # assessment year (final pop model numbers at size are given for July 1, assessment year)","\n",file=conn);
    cat(d$z$n,  "  # Number of size bins in the model","\n",file=conn);
    cat('#size bin cut pts\n',file=conn);
    cat(d$zc$vls,"\n",sep=' ',file=conn);
    cat(d$f$n,"  #number of fisheries\n",file=conn);
    ctr<-1;
    for (f in d$f$nms) {
        str<-gsub("[[:blank:]]","_",f);
        cat(str,"  #fishery ",ctr,': ',f,"\n",sep='',file=conn);
        ctr<-ctr+1;
    }
    cat(d$v$n,"  #number of surveys\n",file=conn);
    ctr<-1;
    for (v in d$v$nms) {
        str<-gsub("[[:blank:]]","_",v);
        cat(str,"  #survey ",ctr,': ',v,"\n",sep='',file=conn);
        ctr<-ctr+1;
    }
    cat("#--THE FOLLOWING ARE DEFAULT INPUTS: PLEASE CHANGE!\n",file=conn)
    cat("FALSE     #run operating model only\n",file=conn);
    cat("TRUE                    #fit to priors\n",file=conn);
    cat("Model.ParametersInfo.dat       #model parameters info file\n",file=conn)
    cat("Model.Datasets.dat             #model datasets file\n",file=conn);
    cat("Model.Options.dat              #model options file\n",file=conn);
    cat("OFF                     #jitter resampling option\n",file=conn);
    cat("0.2                     #jitter range\n",file=conn);
    cat("OFF                     #prior resampling option\n",file=conn);
    cat("1                       #prior variance inflation factor\n",file=conn);
    
    cat("\n\n",file=conn)
    cat("#####################################################################\n",file=conn);
    cat("#TCSAM2015 Model Datasets File                                      #\n",file=conn);
    cat("#####################################################################\n",file=conn);
    cat("Data.BioInfo.dat       #biological info file\n",file=conn)
    cat(d$f$n,"  #number of fishery data files\n",file=conn);
    for (f in d$f$nms) {
        str<-gsub("[[:blank:]]","_",f);
        cat(paste('Data.Fishery.',str,'.dat',sep=''),"  #data file for",f,"\n",file=conn);
    }
    cat(d$v$n,"  #number of survey data files\n",file=conn);
    for (v in d$v$nms) {
        str<-gsub("[[:blank:]]","_",v);
        cat(paste('Data.Survey.',str,'.dat',sep=''),"  #data file for",v,"\n",file=conn);
    }

    cat("\n\n",file=conn)
    cat("#####################################################################\n",file=conn);
    cat("#TCSAM2015 Model Biological Info File                               #\n",file=conn);
    cat("#####################################################################\n",file=conn);
    cat("BIO_DATA     #required keyword\n",file=conn);
    cat("#------------------\n",file=conn);
    cat("#recruitment lag (in years)\n",file=conn);
    cat("5  #recLag\n",file=conn);
    cat("#------------------\n",file=conn);
    cat("#length bins\n",file=conn);
    cat("32 #NUMBER OF SIZE BINS\n",file=conn);
    cat("#SIZE BINS\n",file=conn);
    cat(d$z$nms,"\n",file=conn);
    cat("#------------------\n",file=conn);
    cat("#WEIGHT-AT-SIZE\n",file=conn);
    cat("KG   #units\n",file=conn);
    cat(d$x$n*d$m$n,"    #number of factor combinations (sex x maturity state)\n",file=conn);
    for (x in d$x$nms){
        for (m in d$m$nms){
            cat(toupper(x),toupper(m),"\n",file=conn);
            cat(mp$W_yxmsz[mxy,x,m,1,],"\n",file=conn);
        }
    }
    cat("#------------------\n",file=conn);
    cat("#Timing of fisheries and mating\n",file=conn);            
    cat(mc$params$fish.time,"   #typical midPtFisheries\n",file=conn);
    cat(mc$params$mate.time,"   #typical matingTime\n",file=conn);
    cat("0   #number of atypical years\n",file=conn);
    cat("#data for atypical years\n",file=conn);
    cat("#year   midPtFisheries  matingTime\n",file=conn);
        
    #Fisheries Data
    writeSim.TCSAM.Fisheries(mc,mp,mr,conn,showPlot=showPlot);
    
    #Surveys Data
    writeSim.TCSAM.Surveys(mc,mp,mr,conn,showPlot=showPlot);
    
    return();
}