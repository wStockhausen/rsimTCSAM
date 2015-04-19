#'
#'@title Write model fisheries results to output file.
#'
#'@description Function to write model fisheries results to output file.
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param mr - model results list object
#'@param fnFshs - files to write fishery data to
#'@param out.dir - folder to write survey data to
#'@param showPlot - flag to show plots
#'
#'@return list with retained, discarded, captured numbers and biomass by f,y,x,s
#'
#'@import ggplot2
#'@import reshape2
#'
#'@export
#'
writeSim.TCSAM.Fisheries<-function(mc,mp,mr,fnFshs,out.dir='.',showPlot=TRUE){
    #model dimensions
    d <- mc$dims;
        
    #--total catch abundance/biomass (millions/1000s mt) [NOT mortality]
    cpN_fyxms <-mr$F_list$cpN_fyxms; #captured abundance by f, y, x, m, s
    cpB_fyxms <-mr$F_list$cpB_fyxms; #captured biomass by f, y, x, m, s
    cpN_fyxmsz<-mr$F_list$cpN_fyxmsz;#size-specific numbers captured
    NC_fyxs<-dimArray(mc,'f.y.x.s');  #captured numbers caught by f, y, x, s 
    BC_fyxs<-dimArray(mc,'f.y.x.s');  #captured biomass caught by f, y, x, s     
    ncy_f<-dimArray(mc,'f',val=0);  #years that will be output for total catch data
    ncc_f<-dimArray(mc,'f',val=0);  #number of factor combinations that will be output for total catch data
    #--Retained catch abundance/biomass (millions/1000s mt)
    rmN_fyxms <-mr$F_list$rmN_fyxms; #retained abundance by f, y, x. m, s
    rmB_fyxms <-mr$F_list$rmB_fyxms; #retained biomass by f, y, x, m, s
    rmN_fyxmsz<-mr$F_list$rmN_fyxmsz;#size-specific numbers retained
    NR_fyxs<-dimArray(mc,'f.y.x.s');  #retained numbers caught by f, y, x, s 
    BR_fyxs<-dimArray(mc,'f.y.x.s');  #retained biomass caught by f, y, x, s 
    nry_f<-dimArray(mc,'f',val=0);  #years that will be output for retained catch data
    nrc_f<-dimArray(mc,'f',val=0);  #number of factor combinations that will be output for retained catch data
    #--Discarded catch abundance/biomass (millions/1000s mt) [NOT mortality]
    dsN_fyxms <-mr$F_list$dsN_fyxms; #discarded abundance by f, y, x. m, s
    dsB_fyxms <-mr$F_list$dsB_fyxms; #discarded biomass by f, y, x, m, s
    dsN_fyxmsz<-mr$F_list$dsN_fyxmsz;#size-specific numbers discarded
    ND_fyxs<-dimArray(mc,'f.y.x.s');  #discard numbers caught by f, y, x, s (NOT mortality)
    BD_fyxs<-dimArray(mc,'f.y.x.s');  #discard biomass caught by f, y, x, s (NOT mortality)
    ndy_f<-dimArray(mc,'f',val=0);  #years that will be output for discarded catch data
    ndc_f<-dimArray(mc,'f',val=0);  #number of factor combinations that will be output for discarded catch data
    for (f in d$f$nms){
        for (y in d$y$nms) {
            if (sum(cpN_fyxms[f,y,,,],na.rm=TRUE)>0){ncy_f[f]<-ncy_f[f]+1;}
            if (sum(rmN_fyxms[f,y,,,],na.rm=TRUE)>0){nry_f[f]<-nry_f[f]+1;}
            if (sum(cpN_fyxms[f,y,,,],na.rm=TRUE)>0){ndy_f[f]<-ndy_f[f]+1;}
            for (x in d$x$nms) {
                for (s in d$s$nms) {
                    NC_fyxs[f,y,x,s]<-sum(cpN_fyxms[f,y,x,,s],na.rm=TRUE);
                    BC_fyxs[f,y,x,s]<-sum(cpB_fyxms[f,y,x,,s],na.rm=TRUE);
                    NR_fyxs[f,y,x,s]<-sum(rmN_fyxms[f,y,x,,s],na.rm=TRUE);
                    BR_fyxs[f,y,x,s]<-sum(rmB_fyxms[f,y,x,,s],na.rm=TRUE);
                    ND_fyxs[f,y,x,s]<-sum(dsN_fyxms[f,y,x,,s],na.rm=TRUE);
                    BD_fyxs[f,y,x,s]<-sum(dsB_fyxms[f,y,x,,s],na.rm=TRUE);
                }#s
            }#x
        }#y
        for (x in d$x$nms) {
            for (s in d$s$nms) {
                if (any(NC_fyxs[f,,x,s]>0,na.rm=TRUE)){ncc_f[f]<-ncc_f[f]+1;}
                if (any(NR_fyxs[f,,x,s]>0,na.rm=TRUE)){nrc_f[f]<-nrc_f[f]+1;}
                if (any(ND_fyxs[f,,x,s]>0,na.rm=TRUE)){ndc_f[f]<-ndc_f[f]+1;}
            }#s
        }#x
    }
    
    #write results to data file
    for (f in d$f$nms){
        cat("writing fishery data to '",file.path(out.dir,fnFshs[[f]]),"'\n",sep='');
        conn<-file(file.path(out.dir,fnFshs[[f]]),open="w");
        fsh<-mc$params$fisheries[[f]];
        cat("#####################################################################\n",file=conn);
        cat("#TCSAM2015 Model File for",f,"\n",file=conn);
        cat("#####################################################################\n",file=conn);
        cat("FISHERY_DATA     #required keyword\n",file=conn);
        cat(gsub('[[:blank:]]',"_",f),"    #fishery name\n",file=conn);
        cat("FALSE","   #has effort data?\n",file=conn);
        anyRet<-any(fsh$output$ret$abundance$flag,fsh$output$ret$biomass$flag,fsh$output$ret$sizecomps$flag);
        anyDsc<-any(fsh$output$dsc$abundance$flag,fsh$output$dsc$biomass$flag,fsh$output$dsc$sizecomps$flag);
        anyTot<-any(fsh$output$tot$abundance$flag,fsh$output$tot$biomass$flag,fsh$output$tot$sizecomps$flag);
        cat(anyRet,"   #has retained catch?\n",file=conn);
        cat(anyDsc,"   #has observed discard catch\n",file=conn);
        cat(anyTot,"   #has observed total catch\n",file=conn);
        cat("#------------EFFORT DATA-----------#\n",file=conn);
        cat("#-----no effort data\n",file=conn);
        cat("#------------RETAINED CATCH DATA------------#\n",file=conn);
        if (anyRet){
            #retained catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(fsh$output$ret$abundance$flag,"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(fsh$output$ret$biomass$flag,  "   #has aggregate catch biomass (weight)\n",file=conn);
            cat(fsh$output$ret$sizecomps$flag,"   #has size frequency data\n",file=conn);
            if (fsh$output$ret$abundance$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat(toupper(fsh$output$ret$abundance$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$ret$abundance$errType),"\t\t#likelihood type\n",file=conn);
                cat(nry_f[f],"    	#number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(nrc_f[f],"		#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        if (sum(NR_fyxs[f,,x,s],na.rm=TRUE)>0){
                            cat(toupper(x),"ALL_MATURITY",toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year    value    cv\n",file=conn);
                            for (y in d$y$nms){
                                if (sum(NR_fyxs[f,y,,],na.rm=TRUE)>0) cat(y,NR_fyxs[f,y,x,s],fsh$output$ret$abundance$err,'\n',file=conn);
                            }#y
                        }
                    }#s
                }#x
            }#ret$abundance$flag
            if (fsh$output$ret$biomass$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat(toupper(fsh$output$ret$biomass$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$ret$biomass$errType),"\t\t#likelihood type\n",file=conn);
                cat(nry_f[f],"    	#number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(nrc_f[f],"    	#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        if (sum(NR_fyxs[f,,x,s],na.rm=TRUE)>0){
                            cat(toupper(x),"ALL_MATURITY",toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year    value    cv\n",file=conn);
                            for (y in d$y$nms){
                                if (sum(NR_fyxs[f,y,,],na.rm=TRUE)>0) cat(y,BR_fyxs[f,y,x,s],fsh$output$ret$biomass$err,'\n',file=conn);
                            }#y
                        }
                    }#s
                }#x
            }#ret$biomass$flag
            if (fsh$output$ret$sizecomps$flag){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat(toupper(fsh$output$ret$sizecomps$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$ret$sizecomps$errType),"\t\t#likelihood type\n",file=conn);
                cat(nry_f[f],"     #number of years of data\n",file=conn);
                cat("MILLIONS    #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(nrc_f[f],"   #number of factor combinations\n",file=conn);
                mdfr<-melt(rmN_fyxmsz[f,,,,,],value.name='var');
                ddfr<-dcast(mdfr,x+s+y~z,fun.aggregate=sum,na.rm=TRUE,value.var='var');
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        idx<-(ddfr$x==x)&(ddfr$s==s);
                        dp<-ddfr[idx,];
                        if (sum(dp[,3+(1:d$z$n)],na.rm=TRUE)>0){
                            cat(toupper(x),'ALL_MATURITY',toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year  ss  ",d$z$nms,'\n',file=conn);
                            for (r in 1:nrow(dp)){
                                if (sum(dp[r,3+(1:d$z$n)],na.rm=TRUE)>0){
                                    cat(dp[r,3],fsh$output$ret$sizecomps$err,file=conn);
                                    for (j in 3+(1:d$z$n)) cat(' ',dp[r,j],file=conn);
                                    cat('\n',file=conn)
                                }
                            }#r
                        }
                    }#s
                }#x
            }#ret$sizecomps$flag
        } else {
            cat("#-----no retained catch data\n",file=conn);
        }
        cat("#------------DISCARDED CATCH DATA------------#\n",file=conn);
        if (anyDsc){
            #discarded catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(fsh$output$dsc$abundance$flag,"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(fsh$output$dsc$biomass$flag,  "   #has aggregate catch biomass (weight)\n",file=conn);
            cat(fsh$output$dsc$sizecomps$flag,"   #has size frequency data\n",file=conn);
            if (fsh$output$dsc$abundance$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat(toupper(fsh$output$dsc$abundance$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$dsc$abundance$errType),"\t\t#likelihood type\n",file=conn);
                cat(ndy_f[f],"        #number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(ndc_f[f],"		#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        if (sum(ND_fyxs[f,,x,s],na.rm=TRUE)>0){
                            cat(toupper(x),"ALL_MATURITY",toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year    value    cv\n",file=conn);
                            for (y in d$y$nms){
                                if (sum(ND_fyxs[f,y,,],na.rm=TRUE)>0) cat(y,ND_fyxs[f,y,x,s],fsh$output$dsc$abundance$err,'\n',file=conn);
                            }#y
                        }
                    }#s
                }#x
            }#dsc$abundance$flag
            if (fsh$output$dsc$biomass$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat(toupper(fsh$output$dsc$biomass$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$dsc$biomass$errType),"\t\t#likelihood type\n",file=conn);
                cat(ndy_f[f],"    	#number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(ndc_f[f],"    	#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        if (sum(ND_fyxs[f,,x,s],na.rm=TRUE)>0){
                            cat(toupper(x),"ALL_MATURITY",toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year    value    cv\n",file=conn);
                            for (y in d$y$nms){
                                if (sum(ND_fyxs[f,y,,],na.rm=TRUE)>0) cat(y,BD_fyxs[f,y,x,s],fsh$output$dsc$biomass$err,'\n',file=conn);
                            }#y
                        }
                    }#s
                }#x
            }#dsc$biomass$flag
            if (fsh$output$dsc$sizecomps$flag){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat(toupper(fsh$output$dsc$sizecomps$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$dsc$sizecomps$errType),"\t\t#likelihood type\n",file=conn);
                cat(ndy_f[f],"     #number of years of data\n",file=conn);
                cat("MILLIONS       #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(ndc_f[f],"   #number of factor combinations\n",file=conn);
                mdfr<-melt(dsN_fyxmsz[f,,,,,],value.name='var');
                ddfr<-dcast(mdfr,x+s+y~z,fun.aggregate=sum,na.rm=TRUE,value.var='var');
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        idx<-(ddfr$x==x)&(ddfr$s==s);
                        dp<-ddfr[idx,];
                        if (sum(dp[,3+(1:d$z$n)],na.rm=TRUE)>0){
                            cat(toupper(x),'ALL_MATURITY',toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year  ss  ",d$z$nms,'\n',file=conn);
                            for (r in 1:nrow(dp)){
                                if (sum(dp[r,3+(1:d$z$n)],na.rm=TRUE)>0){
                                    cat(dp[r,3],fsh$output$dsc$sizecomps$err,file=conn);
                                    for (j in 3+(1:d$z$n)) cat(' ',dp[r,j],file=conn);
                                    cat('\n',file=conn)
                                }
                            }#r
                        }
                    }#s
                }#x
            }#dsc$sizecomps$flag
        } else {
            cat("#-----no discarded catch data\n",file=conn);
        }
        cat("#------------TOTAL CATCH DATA------------#\n",file=conn);
        if (anyTot){
            #total catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(fsh$output$tot$abundance$flag,"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(fsh$output$tot$biomass$flag,  "   #has aggregate catch biomass (weight)\n",file=conn);
            cat(fsh$output$tot$sizecomps$flag,"   #has size frequency data\n",file=conn);
            if (fsh$output$tot$abundance$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat(toupper(fsh$output$tot$abundance$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$tot$abundance$errType),"\t\t#likelihood type\n",file=conn);
                cat(ncy_f[f],"        #number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(ncc_f[f],"		#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        if (sum(NC_fyxs[f,,x,s],na.rm=TRUE)>0){
                            cat(toupper(x),"ALL_MATURITY",toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year    value    cv\n",file=conn);
                            for (y in d$y$nms){
                                if (sum(NC_fyxs[f,y,,],na.rm=TRUE)>0) cat(y,NC_fyxs[f,y,x,s],fsh$output$tot$abundance$err,'\n',file=conn);
                            }#y
                        }
                    }#s
                }#x
            }#tot$abundance$flag
            if (fsh$output$tot$biomass$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat(toupper(fsh$output$tot$biomass$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$tot$biomass$errType),"\t\t#likelihood type\n",file=conn);
                cat(ncy_f[f],"        #number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(ncc_f[f],"    	#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        if (sum(NC_fyxs[f,,x,s],na.rm=TRUE)>0){
                            cat(toupper(x),"ALL_MATURITY",toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year    value    cv\n",file=conn);
                            for (y in d$y$nms){
                                if (sum(NR_fyxs[f,y,,],na.rm=TRUE)>0) cat(y,BC_fyxs[f,y,x,s],fsh$output$tot$biomass$err,'\n',file=conn);
                            }#y
                        }
                    }#s
                }#x
            }#tot$biomass$flag
            if (fsh$output$tot$sizecomps$flag){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat(toupper(fsh$output$tot$sizecomps$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$tot$sizecomps$errType),"\t\t#likelihood type\n",file=conn);
                cat(ncy_f[f],"     #number of years of data\n",file=conn);
                cat("MILLIONS     #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(ncc_f[f],"   #number of factor combinations\n",file=conn);
                mdfr<-melt(cpN_fyxmsz[f,,,,,],value.name='var');
                ddfr<-dcast(mdfr,x+s+y~z,fun.aggregate=sum,na.rm=TRUE,value.var='var');
                for (x in d$x$nms){
                    for (s in d$s$nms){
                        idx<-(ddfr$x==x)&(ddfr$s==s);
                        dp<-ddfr[idx,];
                        if (sum(dp[,3+(1:d$z$n)],na.rm=TRUE)>0){
                            cat(toupper(x),'ALL_MATURITY',toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                            cat("#year  ss  ",d$z$nms,'\n',file=conn);
                            for (r in 1:nrow(dp)){
                                if (sum(dp[r,3+(1:d$z$n)],na.rm=TRUE)>0){
                                    cat(dp[r,3],fsh$output$tot$sizecomps$err,file=conn);
                                    for (j in 3+(1:d$z$n)) cat(' ',dp[r,j],file=conn);
                                    cat('\n',file=conn)
                                }
                            }#r
                        }
                    }#s
                }#x
            }
        } else {
            cat("#-----no total catch data\n",file=conn);
        }
        close(conn);
    }#f
    return(invisible(list(NR_fyxs=NR_fyxs,BR_fyxs=BR_fyxs,
                          ND_fyxs=ND_fyxs,BD_fyxs=BD_fyxs,
                          NC_fyxs=NC_fyxs,BC_fyxs=BC_fyxs)));
}
    
    