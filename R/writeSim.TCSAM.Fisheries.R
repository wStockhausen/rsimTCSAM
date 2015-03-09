#'
#'@title Write model fisheries results to output file.
#'
#'@description Function to write model fisheries results to output file.
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param mr - model results list object
#'@param fnSrvs - files to write fishery data to
#'@param showPlot - flag to show plots
#'
#'@import ggplot2
#'@import reshape2
#'
#'@export
#'
writeSim.TCSAM.Fisheries<-function(mc,mp,mr,fnFshs,showPlot=TRUE){
    #model dimensions
    d <- mc$dims;
        
    #--total catch abundance/biomass (millions/1000s mt) [NOT mortality]
    cpN_fyxms <-mr$F_list$cpN_fyxms; #captured abundance by f, y, x, m, s
    cpB_fyxms <-mr$F_list$cpB_fyxms; #captured biomass by f, y, x, m, s
    cpN_fyxmsz<-mr$F_list$cpN_fyxmsz;#size-specific numbers captured
    NC_fyx<-dimArray(mc,'f.y.x');  #captured numbers caught by f, x, y 
    BC_fyx<-dimArray(mc,'f.y.x');  #captured biomass caught by f, x, y 
    nr_f<-dimArray(mc,'f');#number of data rows that will be output for captured catch data
    for (f in d$f$nms){
        for (y in d$y$nms) {
            for (x in d$x$nms) {
                NC_fyx[f,y,x]<-sum(cpN_fyxms[f,y,x,,]);
                BC_fyx[f,y,x]<-sum(cpB_fyxms[f,y,x,,]);
            }
        }
    }
        
    #--Retained catch abundance/biomass (millions/1000s mt)
    rmN_fyxms <-mr$F_list$rmN_fyxms; #retained abundance by f, y, x. m, s
    rmB_fyxms <-mr$F_list$rmB_fyxms; #retained biomass by f, y, x, m, s
    rmN_fyxmsz<-mr$F_list$rmN_fyxmsz;#size-specific numbers retained
    NR_fyx<-dimArray(mc,'f.y.x');  #retained numbers caught by f, x, y 
    BR_fyx<-dimArray(mc,'f.y.x');  #retained biomass caught by f, x, y 
    nr_f<-dimArray(mc,'f');#number of data rows that will be output for retained catch data
    for (f in d$f$nms){
        for (y in d$y$nms) {
            for (x in d$x$nms) {
                NR_fyx[f,y,x]<-sum(rmN_fyxms[f,y,x,,]);
                BR_fyx[f,y,x]<-sum(rmB_fyxms[f,y,x,,]);
            }
            if (!is.na(NR_fyx[f,y,1])){nr_f[f]<-nr_f[f]+1;}
        }
    }
    
    #--Discarded catch abundance/biomass (millions/1000s mt) [NOT mortality]
    dsN_fyxms <-mr$F_list$dsN_fyxms; #discarded abundance by f, y, x. m, s
    dsB_fyxms <-mr$F_list$dsB_fyxms; #discarded biomass by f, y, x, m, s
    dsN_fyxmsz<-mr$F_list$dsN_fyxmsz;#size-specific numbers discarded

    ND_fyx<-dimArray(mc,'f.y.x');  #discard numbers caught by f, x, y (NOT mortality)
    BD_fyx<-dimArray(mc,'f.y.x');  #discard biomass caught by f, x, y (NOT mortality)
    nd_fx<-dimArray(mc,'f.x');     #number of rows to be output for discard catch data
    for (f in d$f$nms){
        for (y in d$y$nms) {
            for (x in d$x$nms){
                ND_fyx[f,y,x]<-sum(rmN_fyxms[f,y,x,,]);
                BD_fyx[f,y,x]<-sum(rmB_fyxms[f,y,x,,]);
                if (sum(ND_fyx[f,y,x],na.rm=TRUE)>0){nd_fx[f,x]<-nd_fx[f,x]+1;}
            }#x
        }#y
    }#f
    
    #write results to data file
    for (f in d$f$nms){
        cat("writing fishery data to '",fnFshs[[f]],"'\n",sep='');
        conn<-file(fnFshs[[f]],open="w");
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
        cat(anyTot,"    #has observed total catch\n",file=conn);
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
                cat(nr_f[f],"    	#number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat("1		#number of factor combinations\n",file=conn);
                cat("MALE ALL_MATURITY ALL_SHELL\n",file=conn);
                cat("#year    value	cv_m\n",file=conn);
                for (y in d$y$nms){
                    if (!is.na(NR_fyx[f,y,1])) cat(y,NR_fyx[f,y,1],fsh$output$ret$abundance$err,'\n',sep='  ',file=conn);
                }#y
            }
            if (fsh$output$ret$biomass$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat(toupper(fsh$output$ret$biomass$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$ret$biomass$errType),"\t\t#likelihood type\n",file=conn);
                cat(nr_f[f],"    	#number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat("1		#number of factor combinations\n",file=conn);
                cat("MALE ALL_MATURITY ALL_SHELL\n",file=conn);
                cat("#year    value	cv_m\n",file=conn);
                for (y in d$y$nms){
                    if (!is.na(BR_fyx[f,y,1]))  cat(y,BR_fyx[f,y,1],fsh$output$ret$biomass$err,'\n',sep='  ',file=conn);
                }#y
            }
            if (fsh$output$ret$sizecomps$flag){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat(toupper(fsh$output$ret$sizecomps$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$ret$sizecomps$errType),"\t\t#likelihood type\n",file=conn);
                cat(nr_f[f],"     #number of years of data\n",file=conn);
                cat("???         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(d$s$n,"   #number of shell factor combinations\n",file=conn);
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
            }
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
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"	#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value	cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(ND_fyx[f,y,x])) cat(y,ND_fyx[f,y,x],fsh$output$dsc$abundance$err,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$dsc$biomass$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat(toupper(fsh$output$dsc$biomass$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$dsc$biomass$errType),"\t\t#likelihood type\n",file=conn);
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value    cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(BD_fyx[f,y,x])) cat(y,BD_fyx[f,y,x],fsh$output$dsc$biomass$err,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$dsc$sizecomps$flag){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat(toupper(fsh$output$dsc$sizecomps$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$dsc$sizecomps$errType),"\t\t#likelihood type\n",file=conn);
                cat(nd_fx[f,1],"     #number of years of data\n",file=conn);
                cat("???         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(d$x$n*d$s$n,"   #number of factor combinations\n",file=conn);
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
            }
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
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value	cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(NC_fyx[f,y,x])) cat(y,NC_fyx[f,y,x],fsh$output$tot$abundance$err,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$tot$biomass$flag){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat(toupper(fsh$output$tot$biomass$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$tot$biomass$errType),"\t\t#likelihood type\n",file=conn);
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value    cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(BC_fyx[f,y,x])) cat(y,BC_fyx[f,y,x],fsh$output$tot$biomass$err,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$tot$sizecomps$flag){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat(toupper(fsh$output$tot$sizecomps$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(fsh$output$tot$sizecomps$errType),"\t\t#likelihood type\n",file=conn);
                cat(nd_fx[f,1],"     #number of years of data\n",file=conn);
                cat("???         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(d$x$n*d$s$n,"   #number of factor combinations\n",file=conn);
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
#     return(invisible(list(NR_fyx=NR_fyx,BR_fyx=BR_fyx,
#                           ND_fyx=ND_fyx,BD_fyx=BD_fyx,
#                           NT_fyx=NT_fyx,BT_fyx=BT_fyx)));
}
    
    