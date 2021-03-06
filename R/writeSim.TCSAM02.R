#'
#'@title Write simulation results to TCSAM02 input file2.
#'
#'@description Function to write simulation results to TCSAM02 input files.
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param mr - model results list object
#'@param out.dir - folder for output files
#'@param showPlot - flag to show plots
#'
#'@return list with elements F_list and S_list
#'
#'@export
#'
writeSim.TCSAM02<-function(mc,mp,mr,out.dir='.',showPlot=TRUE){
    d<-mc$dims;
    mxy<-as.character(d$y$mxy);
    
    #create data file names
    Model.Config        <-'Model.Config.dat';
    Model.ParametersInfo<-'Model.ParametersInfo.dat';
    Model.Datasets      <-'Model.Datasets.dat';
    Model.Options       <-'Model.Options.dat';
    Model.Data.BioInfo  <-'Model.Data.BioInfo.dat';
    
    fnFshs<-list();
    for (f in d$f$nms) {
        str<-gsub("[[:blank:]]","_",f);
        fnFshs[[f]]<-paste('Model.Data.Fishery.',str,'.dat',sep='');
    }
    
    fnSrvs<-list();
    for (v in d$v$nms) {
        str<-gsub("[[:blank:]]","_",v);
        fnSrvs[[v]]<-paste('Model.Data.Survey.',str,'.dat',sep='');
    }
    
    #write out Model Configuration file
    conn<-file(file.path(out.dir,Model.Config),open="w");    
    cat("#####################################################################\n",file=conn);
    cat("#TCSAM02 Model Configuration File                                 #\n",file=conn);
    cat("#####################################################################\n",file=conn);
    cat('rsimTest',"  # model scenario","\n",file=conn);
    cat(d$y$mny,"  # model start year (pop. model start year)","\n",file=conn);
    cat(d$y$asy,"  # assessment year (final pop model numbers at size are given for July 1, assessment year)","\n",file=conn);
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
    cat("#--THE FOLLOWING ARE DEFAULT INPUTS: PLEASE CHANGE AS NECESSARY!\n",file=conn)
    cat("FALSE     #run operating model only\n",file=conn);
    cat("TRUE                    #fit to priors\n",file=conn);
    cat(Model.ParametersInfo,"\t\t\t#model parameters info file\n",file=conn);
    cat(Model.Datasets,"\t\t\t#model datasets file\n",file=conn);
    cat(Model.Options,"\t\t\t#model options file\n",file=conn);
    cat("ON                      #jitter resampling option\n",file=conn);
    cat("0.2                     #jitter range\n",file=conn);
    cat("OFF                     #prior resampling option\n",file=conn);
    cat("1                       #prior variance inflation factor\n",file=conn);
    close(conn);
    
    #write out Model Datasets file
    conn<-file(file.path(out.dir,Model.Datasets),open="w");    
    cat("#####################################################################\n",file=conn);
    cat("#TCSAM02 Model Datasets File                                      #\n",file=conn);
    cat("#####################################################################\n",file=conn);
    cat(Model.Data.BioInfo,"\t\t\t#biological info file\n",file=conn)
    cat(d$f$n,"  #number of fishery data files\n",file=conn);
    for (f in d$f$nms) {
        cat(fnFshs[[f]],"\t\t\t#data file for",f,"\n",file=conn);
    }
    cat(d$v$n,"  #number of survey data files\n",file=conn);
    for (v in d$v$nms) {
        cat(fnSrvs[[v]],"\t\t\t#data file for",v,"\n",file=conn);
    }
    close(conn);
    
    #write out biological info file
    conn<-file(file.path(out.dir,Model.Data.BioInfo),open="w");    
    cat("\n\n",file=conn)
    cat("#####################################################################\n",file=conn);
    cat("#TCSAM02 Model Biological Info File                               #\n",file=conn);
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
    close(conn);
        
    #Fisheries Data
    fshs<-writeSim.TCSAM02.Fisheries(mc,mp,mr,fnFshs,out.dir=out.dir,showPlot=showPlot);
    
    #Surveys Data
    srvs<-writeSim.TCSAM02.Surveys(mc,mp,mr,fnSrvs,out.dir=out.dir,showPlot=showPlot);
    
    #Parameters info
    fn<-file.path(out.dir,'rsim.ParametersInfo.dat')
    conn<-file(fn,open="w");    
    cat("#---PARAMETERS INFO from rsimTCSAM----------\n",file=conn);
    cat("#--recruitment\n",file=conn);
    blocks<-mc$params$rec$blocks;
    cat("blocks:",names(blocks),'\n',file=conn)
    for (t in names(blocks)){
        tb<-blocks[[t]];
        cat(t,':\n',sep='',file=conn);
        cat('years:',tb$years,'\n',sep='\t',file=conn);
        cat('devs:',mp$R_list$devs_y[as.character(tb$years)],'\n',sep='\t',file=conn);
    }#t
    cat("#--molt-to-maturity\n",file=conn);
    cat("#sex","shell condition","years",mc$dims$z$nms,'\n',sep='\t',file=conn)
    blocks<-mc$params$moltToMaturity$blocks;
    cat("blocks:",names(blocks),'\n',file=conn)
    for (t in names(blocks)){
        tb<-blocks[[t]];
        yr<-as.character(tb$years[1]);#get first year of block
        cat(t,'\n',file=conn)
        for (x in mc$dims$x$nms){
            for (s in mc$dims$s$nms){
                lgt<-log(mp$prM2M_yxsz[yr,x,s,]/(1-mp$prM2M_yxsz[yr,x,s,]));
                cat(x,s,':',lgt,'\n',sep='\t',file=conn);
            }#s
        }#x        
    }#t
    cat("#--fisheries\n",file=conn);
    fs<-names(mc$params$fisheries);
    for (f in fs){
        cat("fishery:",f,'\n',file=conn)
        blocks<-mc$params$fisheries[[f]]$blocks;
        cat("blocks:",names(blocks),'\n',file=conn)
        for (t in names(blocks)){
            tb<-blocks[[t]];
            cat(t,':\n',sep='',file=conn);
            cat('years:',tb$years,'\n',sep='\t',file=conn);
            cat('devs:',mp$F_list$devs_fy[f,as.character(tb$years)],'\n',sep='\t',file=conn);
        }#t
    }#f
    close(conn);
    
    #OFL info
    fn<-file.path(out.dir,'rsim.OFLInfo.dat')
    conn<-file(fn,open="w");    
    cat("#---OFL INFO from rsimTCSAM----------\n",file=conn);
    cat("#--recruitment\n",file=conn);
    mny<-as.character(mc$dims$y$mny);
    mxy<-as.character(mc$dims$y$mxy);
    asy<-as.character(mc$dims$y$asy);
    yrs<-as.character(1982:mc$dims$y$mxy);
    avgRM<-mean(mp$R_list$R_y[yrs]*mp$R_list$R_yx[yrs,  'male'],na.rm=TRUE);
    avgRF<-mean(mp$R_list$R_y[yrs]*mp$R_list$R_yx[yrs,'female'],na.rm=TRUE);
    cat("avgRM = ",avgRM,"\n",file=conn);
    cat("avgRF = ",avgRF,"\n",file=conn);
    cat("MMB100 = ",mr$P_list$MB_yx[mny,'male'],"\n",file=conn);
    cat("MMBcur = ",mr$P_list$MB_yx[mxy,'male'],"\n",file=conn);
    cat("MMBrat = ",mr$P_list$MB_yx[mxy,'male']/mr$P_list$MB_yx[mny,'male'],"\n",file=conn);
    cat("init n-at-z\n",file=conn);
    for (x in mc$dims$x$nms) {
        for (m in mc$dims$m$nms) {
            for (s in mc$dims$s$nms) {
                cat(x,m,s,paste(mr$iN_xmsz[x,m,s,],collapse=" "),"\n",file=conn);
            }
        }
    }
    cat("final n-at-z\n",file=conn);
    for (x in mc$dims$x$nms) {
        for (m in mc$dims$m$nms) {
            for (s in mc$dims$s$nms) {
                cat(x,m,s,paste(mr$P_list$N_yxmsz[asy,x,m,s,],collapse=" "),"\n",file=conn);
            }
        }
    }
    close(conn);
    
    return(list(F_list=fshs,S_list=srvs));
}