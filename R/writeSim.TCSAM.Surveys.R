#'
#'@title Write model survey results to output file.
#'
#'@description Function to write model survey results to output file.
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param mr - model results list object
#'@param fnSrvs - files to write survey data to
#'@param out.dir - folder to write survey data to
#'@param showPlot - flag to show plots
#'
#'@return NULL 
#'
#'@import ggplot2
#'@import reshape2
#'
#'@export
#'
writeSim.TCSAM.Surveys<-function(mc,mp,mr,fnSrvs,out.dir='.',showPlot=TRUE){
    #model dimensions
    d <- mc$dims;
    
    #--Survey abundance/biomass numbers
    N_vyxms<-mr$S_list$N_vyxms; #survey abundance by y,x,m,s
    B_vyxms<-mr$S_list$B_vyxms; #survey biomass by y,x,m,s
    N_vyxmsz<-mr$S_list$N_vyxmsz;#survey abundance by y,x,m,s,z
    
    ny_v<-dimArray(mc,'v');#number of years of 'data'
    nc_v<-dimArray(mc,'v');#number of factor combinations for 'data'
    for (v in 1:d$v$n){
        for (y in d$y$nms) {
            if (sum(N_vyxms[v,y,,,],na.rm=TRUE)>0){ny_v[v]<-ny_v[v]+1;}
        }#y
        for (x in d$x$nms) {
            for (m in d$m$nms) {
                for (s in d$s$nms) {
                    if (any(N_vyxms[v,,x,m,s]>0,na.rm=TRUE)){nc_v[v]<-nc_v[v]+1;}
                }#s
            }#m
        }#x
    }#v
    for (v in d$v$nms){
        cat("writing survey data to '",file.path(out.dir,fnSrvs[[v]]),"'\n",sep='');
        conn<-file(file.path(out.dir,fnSrvs[[v]]),open="w");
        srv<-mc$params$surveys[[v]];
        cat("\n\n",file=conn)
        cat("#####################################################################\n",file=conn);
        cat("#TCSAM2015 Model File for",v,"\n",file=conn);
        cat("#####################################################################\n",file=conn);
        cat("SURVEY_DATA     #required keyword\n",file=conn);
        cat(gsub('[[:blank:]]',"_",v),"    #survey name\n",file=conn);
        cat("#------------SURVEY CATCH DATA------------#\n",file=conn);
        if (srv$output$abundance$flag|srv$output$biomass$flag|srv$output$sizecomps$flag){
            #total catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(srv$output$abundance$flag,"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(srv$output$biomass$flag,  "   #has aggregate catch biomass (weight)\n",file=conn);
            cat(srv$output$sizecomps$flag,"   #has size frequency data\n",file=conn);
            if (srv$output$abundance$flag){
                cv     <-srv$output$abundance$err;
                errType<-srv$output$abundance$errType;
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat(toupper(srv$output$abundance$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(srv$output$abundance$errType),"\t\t#likelihood type\n",file=conn);
                cat(ny_v[v],"  #number of years\n",file=conn);
                cat("MILLIONS   #catch (numbers) units\n",file=conn);
                cat(nc_v[v],"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (m in d$m$nms){
                        for (s in d$s$nms){
                            if (sum(N_vyxms[v,,x,m,s],na.rm=TRUE)>0){
                                cat(toupper(x),toupper(gsub('[[:blank:]]','_',m)),toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                                cat("#year    value    cv\n",file=conn);
                                for (y in d$y$nms){
                                    val <- N_vyxms[v,y,x,m,s];
                                    if (!is.na(val)) {
                                        if (srv$output$abundance$addErr){
                                            val <- addError(val,cv=cv,type=errType);
                                        }
                                        cat(y,val,cv,'\n',file=conn);
                                    }
                                }#y
                            }
                        }#s
                    }#m
                }#x
            }#abundance$flag
            if (srv$output$biomass$flag){
                cv     <-srv$output$biomass$err;
                errType<-srv$output$biomass$errType;
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat(toupper(srv$output$biomass$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(srv$output$biomass$errType),"\t\t#likelihood type\n",file=conn);
                cat(ny_v[v],"  #number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(nc_v[v],"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    for (m in d$m$nms){
                        for (s in d$s$nms){
                            if (sum(N_vyxms[v,,x,m,s],na.rm=TRUE)>0){
                                cat(toupper(x),toupper(gsub('[[:blank:]]','_',m)),toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                                cat("#year    value    cv\n",file=conn);
                                for (y in d$y$nms){
                                    val <- B_vyxms[v,y,x,m,s];
                                    if (!is.na(val)) {
                                        if (srv$output$biomass$addErr){
                                            val <- addError(val,cv=cv,type=errType);
                                        }
                                        cat(y,val,cv,'\n',file=conn);
                                    }
                                }#y
                            }
                        }#s
                    }#m
                }#x
            }#biomass$flag
            if (srv$output$sizecomps$flag){
                cv     <-srv$output$sizecomps$err;
                errType<-srv$output$sizecomps$errType;
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat(toupper(srv$output$sizecomps$aggType),"\t\t#objective function fitting option\n",file=conn);
                cat(toupper(srv$output$sizecomps$errType),"\t\t#likelihood type\n",file=conn);
                cat(ny_v[v],"     #number of years of data\n",file=conn);
                cat("MILLIONS         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(nc_v[v],"   #number of factor combinations\n",file=conn);
                mdfr<-melt(N_vyxmsz[v,,,,,],value.name='var');
                ddfr<-dcast(mdfr,x+m+s+y~z,fun.aggregate=sum,na.rm=TRUE,value.var='var');
                for (x in d$x$nms){
                    for (m in d$m$nms){
                        for (s in d$s$nms){
                            idx<-(ddfr$x==x)&(ddfr$m==m)&(ddfr$s==s);
                            dp<-ddfr[idx,];
                            if (sum(dp[,4+(1:d$z$n)],na.rm=TRUE)>0){
                                cat(toupper(x),toupper(gsub('[[:blank:]]','_',m)),toupper(gsub('[[:blank:]]','_',s)),'\n',file=conn);
                                cat("#year  ss  ",d$z$nms,'\n',file=conn);
                                ss<-srv$output$sizecomps$err
                                for (r in 1:nrow(dp)){
                                    vals<-dp[r,4+(1:d$z$n)];
                                    if (sum(vals,na.rm=TRUE)>0){
                                        vals<-addError(vals,ss=ss,type=errType)
                                        cat(dp[r,4],ss,' ',file=conn);
                                        #for (j in 1:d$z$n) cat(' ',vals[j],file=conn);
                                        cat(vals,file=conn);
                                        cat('\n',file=conn)
                                    }
                                }#r
                            }
                        }#s
                    }#m
                }#x
            }#sizecomps$flag
        } else {
            cat("#-----no survey catch data\n",file=conn);
        }
        close(conn);
    }#v
    return(NULL)
}