#'
#'@title Write model fisheries results to output file.
#'
#'@description Function to write model fisheries results to output file.
#'
#'@param mc - model configuration list object
#'@param mp - model processes list object
#'@param mr - model results list object
#'@param conn - open file connection for output
#'@param showPlot - flag to show plots
#'
#'@import ggplot2
#'@import reshape2
#'
#'@export
#'
writeSim.TCSAM.Fisheries<-function(mc,mp,mr,conn,showPlot=TRUE){
    #model dimensions
    d <- mc$dims;
    #--Retained catch abundance/biomass (1000s mt)
    NR_fy<-dimArray(mc,'f.y');     #retained abundance by f, y
    BR_fy<-dimArray(mc,'f.y');     #retained biomass by f, y
    W_yxmsz<-mp$W_yxmsz;           #weight-at-size retained
    NR_fyxmsz<-mr$F_list$NR_fyxmsz;#numbers retained
    for (f in d$f$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        NR_fy[f,y]<-NR_fy[f,y]+sum(NR_fyxmsz[f,y,x,m,s,]);
                        BR_fy[f,y]<-BR_fy[f,y]+sum(W_yxmsz[y,x,m,s,]*NR_fyxmsz[f,y,x,m,s,]);
                    }#s
                }#m
            }#x
        }#y
    }#f
    if (showPlot){
        mdfr<-melt(NR_fy,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=f,shape=f),data=mdfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Retained Catch Abundance (millions)")
        p <- p + guides(color=guide_legend('fishery',order=1),shape=guide_legend(''))
        print(p);
        mdfr<-melt(BR_fy,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=f,shape=f),data=mdfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Retained Catch Biomass (1000s t)")
        p <- p + guides(color=guide_legend('fishery',order=1),shape=guide_legend(''))
        print(p);
    }    
    #calc number of data rows that will be output for retained catch data
    nr_f<-dimArray(mc,'f');
    for (f in d$f$nms){
        for (y in d$y$nms) {if (!is.na(NR_fy[f,y])){nr_f[f]<-nr_f[f]+1;}}
    }
    
    #--Discard catch numbers
    ND_fyx<-dimArray(mc,'f.y.x');  #discard numbers by f, x, y (NOT mortality)
    BD_fyx<-dimArray(mc,'f.y.x');  #discard biomass by f, x, y (NOT mortality)
    NC_fyxmsz<-mr$F_list$NC_fyxmsz;#numbers captured
    for (f in d$f$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        ND_fyx[f,y,x]<-ND_fyx[f,y,x]+sum(NC_fyxmsz[f,y,x,m,s,]-NR_fyxmsz[f,y,x,m,s,]);
                        BD_fyx[f,y,x]<-BD_fyx[f,y,x]+sum(W_yxmsz[y,x,m,s,]*(NC_fyxmsz[f,y,x,m,s,]-NR_fyxmsz[f,y,x,m,s,]));
                    }#s
                }#m
            }#x
        }#y
    }#f
    if (showPlot){
        mdfr<-melt(ND_fyx,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=f,shape=x),data=mdfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Discard Catch Abundance (millions)")
        p <- p + guides(color=guide_legend('fishery',order=1),shape=guide_legend('sex'))
        print(p);
        mdfr<-melt(BD_fyx,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=f,shape=x),data=mdfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Discard Catch Biomass (1000s t)")
        p <- p + guides(color=guide_legend('fishery',order=1),shape=guide_legend('sex'))
        print(p);
    }
    #calc number of rows to be output for discard catch data
    nd_fx<-dimArray(mc,'f.x');
    for (f in d$f$nms){
        for (y in d$y$nms) {
            for (x in d$x$nms){
                if (sum(ND_fyx[f,y,x],na.rm=TRUE)>0){
                    nd_fx[f,x]<-nd_fx[f,x]+1;
                }
            }#x
        }#y
    }#f
    for (f in d$f$nms){
        fsh<-mc$params$fisheries[[f]];
        cat("\n\n",file=conn)
        cat("#####################################################################\n",file=conn);
        cat("#TCSAM2015 Model File for",f,"\n",file=conn);
        cat("#####################################################################\n",file=conn);
        cat("FISHERY_DATA     #required keyword\n",file=conn);
        cat(gsub('[[:blank:]]',"_",f),"    #fishery name\n",file=conn);
        cat("FALSE","   #has effort data?\n",file=conn);
        cat(any(fsh$output$ret),"   #has retained catch?\n",file=conn);
        cat(any(fsh$output$dsc),"   #has observed discard catch\n",file=conn);
        cat(any(fsh$output$tot),"    #has observed total catch\n",file=conn);
        cat("#------------EFFORT DATA-----------#\n",file=conn);
        cat("#-----no effort data\n",file=conn);
        cat("#------------RETAINED CATCH DATA------------#\n",file=conn);
        if (any(fsh$output$ret)){
            #retained catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(fsh$output$ret[1],"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(fsh$output$ret[2],"   #has aggregate catch biomass (weight)\n",file=conn);
            cat(fsh$output$ret[3],"   #has size frequency data\n",file=conn);
            if (fsh$output$ret[1]){
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat("BY_TOTAL     #objective function fitting option\n",file=conn);
                cat("NORM2        #likelihood type\n",file=conn);
                cat(nr_f[f],"    	#number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat("1		#number of factor combinations\n",file=conn);
                cat("MALE ALL_MATURITY ALL_SHELL\n",file=conn);
                cat("#year    value	cv_m\n",file=conn);
                for (y in d$y$nms){
                    if (!is.na(NR_fy[f,y])) cat(y,NR_fy[f,y],0.01,'\n',sep='  ',file=conn);
                }#y
            }
            if (fsh$output$ret[2]){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat("BY_TOTAL     #objective function fitting option\n",file=conn);
                cat("NORM2        #likelihood type\n",file=conn);
                cat(nr_f[f],"    	#number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat("1		#number of factor combinations\n",file=conn);
                cat("MALE ALL_MATURITY ALL_SHELL\n",file=conn);
                cat("#year    value	cv_m\n",file=conn);
                for (y in d$y$nms){
                    if (!is.na(BR_fy[f,y]))  cat(y,BR_fy[f,y],0.02,'\n',sep='  ',file=conn);
                }#y
            }
            if (fsh$output$ret[3]){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat("BY_TOTAL    #objective function fitting option\n",file=conn);
                cat("MULTINOMIAL #likelihood type\n",file=conn);
                cat(nr_f[f],"     #number of years of data\n",file=conn);
                cat("???         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(d$s$n,"   #number of shell factor combinations\n",file=conn);
                mdfr<-melt(NR_fyxmsz[f,,,,,],value.name='var');
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
                                    cat(dp[r,3],50,file=conn);
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
        if (any(fsh$output$dsc)){
            #discarded catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(fsh$output$dsc[1],"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(fsh$output$dsc[2],"   #has aggregate catch biomass (weight)\n",file=conn);
            cat(fsh$output$dsc[3],"   #has size frequency data\n",file=conn);
            if (fsh$output$dsc[1]){
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat("BY_SEX     #objective function fitting option\n",file=conn);
                cat("NORM2        #likelihood type\n",file=conn);
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"	#number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value	cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(ND_fyx[f,y,x])) cat(y,ND_fyx[f,y,x],0.05,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$dsc[2]){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat("BY_TOTAL     #objective function fitting option\n",file=conn);
                cat("NORM2        #likelihood type\n",file=conn);
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value    cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(BD_fyx[f,y,x])) cat(y,BD_fyx[f,y,x],0.05,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$dsc[3]){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat("BY_SEX    #objective function fitting option\n",file=conn);
                cat("MULTINOMIAL #likelihood type\n",file=conn);
                cat(nd_fx[f,1],"     #number of years of data\n",file=conn);
                cat("???         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(d$x$n*d$s$n,"   #number of factor combinations\n",file=conn);
                mdfr<-melt(NC_fyxmsz[f,,,,,]-NR_fyxmsz[f,,,,,],value.name='var');
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
                                    cat(dp[r,3],50,file=conn);
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
        if (any(fsh$output$tot)){
            #total catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(fsh$output$tot[1],"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(fsh$output$tot[2],"   #has aggregate catch biomass (weight)\n",file=conn);
            cat(fsh$output$tot[3],"   #has size frequency data\n",file=conn);
            if (fsh$output$tot[1]){
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat("BY_SEX     #objective function fitting option\n",file=conn);
                cat("NORM2        #likelihood type\n",file=conn);
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value	cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(ND_fyx[f,y,x])) cat(y,ND_fyx[f,y,x],0.05,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$tot[2]){
                cat("#------------AGGREGATE CATCH ABUNDANCE (BIOMASS)------------#\n",file=conn);
                cat("AGGREGATE_BIOMASS #required keyword\n",file=conn);
                cat("BY_TOTAL     #objective function fitting option\n",file=conn);
                cat("NORM2        #likelihood type\n",file=conn);
                cat(nd_fx[f,1],"  #number of years\n",file=conn);
                cat("THOUSANDS_MT         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value    cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(BD_fyx[f,y,x])) cat(y,BD_fyx[f,y,x],0.05,'\n',file=conn);
                    }#y
                }#x
            }
            if (fsh$output$tot[3]){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat("BY_SEX    #objective function fitting option\n",file=conn);
                cat("MULTINOMIAL #likelihood type\n",file=conn);
                cat(nd_fx[f,1],"     #number of years of data\n",file=conn);
                cat("???         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(d$x$n*d$s$n,"   #number of factor combinations\n",file=conn);
                mdfr<-melt(NC_fyxmsz[f,,,,,],value.name='var');
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
                                    cat(dp[r,3],50,file=conn);
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
    }#f
}
    
    