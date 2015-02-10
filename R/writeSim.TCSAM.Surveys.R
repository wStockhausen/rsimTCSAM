#'
#'@title Write model survey results to output file.
#'
#'@description Function to write model survey results to output file.
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
writeSim.TCSAM.Surveys<-function(mc,mp,mr,conn,showPlot=TRUE){
    #model dimensions
    d <- mc$dims;
    #--Survey abundance numbers
    N_vyx<-dimArray(mc,'v.y.x'); #survey numbers by y, x
    B_vyx<-dimArray(mc,'v.y.x'); #survey biomass by y, x
    W_yxmsz<-mp$W_yxmsz;         #weight-at-size
    N_vyxmsz<-mr$N_vyxmsz;       #survey abundance by y,x,m,s,z
    for (v in d$v$nms){
        for (y in d$y$nms){
            for (x in d$x$nms){
                for (m in d$m$nms){
                    for (s in d$s$nms){
                        N_vyx[v,y,x]<-N_vyx[v,y,x]+sum(N_vyxmsz[v,y,x,m,s,]);
                        B_vyx[v,y,x]<-B_vyx[v,y,x]+sum(W_yxmsz[y,x,m,s,]*N_vyxmsz[v,y,x,m,s,]);
                    }#s
                }#m
            }#x
        }#y
    }#v
    if (showPlot){
        mdfr<-melt(N_vyx,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=v,shape=x),data=mdfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Survey Abundance (millions)")
        p <- p + guides(color=guide_legend('survey',order=1),shape=guide_legend('sex'))
        print(p);
        mdfr<-melt(B_vyx,value.name='val');
        p <- ggplot(aes(x=y,y=val,color=v,shape=x),data=mdfr);
        p <- p + geom_point();
        p <- p + geom_line();
        p <- p + labs(x='year',y="Survey Biomass (1000s t)")
        p <- p + guides(color=guide_legend('survey',order=1),shape=guide_legend('sex'))
        print(p);
    }
    n_v<-dimArray(mc,'v');
    for (v in 1:d$v$n){
        for (y in d$y$nms) {
            if (sum(N_vyx[v,y,],na.rm=TRUE)>0){n_v[v]<-n_v[v]+1;}
        }
    }
    for (v in d$v$nms){
        srv<-mc$params$surveys[[v]];
        cat("\n\n",file=conn)
        cat("#####################################################################\n",file=conn);
        cat("#TCSAM2015 Model File for",v,"\n",file=conn);
        cat("#####################################################################\n",file=conn);
        cat("SURVEY_DATA     #required keyword\n",file=conn);
        cat(gsub('[[:blank:]]',"_",v),"    #survey name\n",file=conn);
        cat("#------------SURVEY CATCH DATA------------#\n",file=conn);
        if (any(srv$output)){
            #total catch
            cat("CATCH_DATA  #required keyword\n",file=conn);
            cat(srv$output[1],"   #has aggregate catch abundance (numbers)\n",file=conn);
            cat(srv$output[2],"   #has aggregate catch biomass (weight)\n",file=conn);
            cat(srv$output[3],"   #has size frequency data\n",file=conn);
            if (srv$output[1]){
                cat("#------------AGGREGATE CATCH ABUNDANCE (NUMBERS)------------#\n",file=conn);
                cat("AGGREGATE_ABUNDANCE #required keyword\n",file=conn);
                cat("BY_SEX     #objective function fitting option\n",file=conn);
                cat("NORM2        #likelihood type\n",file=conn);
                cat(n_v[v],"  #number of years\n",file=conn);
                cat("MILLIONS         #catch (numbers) units\n",file=conn);
                cat(d$x$n,"    #number of factor combinations\n",file=conn);
                for (x in d$x$nms){
                    cat(toupper(x),"ALL_MATURITY ALL_SHELL\n",file=conn);
                    cat("#year    value    cv_m\n",file=conn);
                    for (y in d$y$nms){
                        if (!is.na(N_vyx[v,y,x])) cat(y,N_vyx[v,y,x],srv$error[1],'\n',file=conn);
                    }#y
                }#x
            }
            if (srv$output[2]){
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
                        if (!is.na(B_vyx[v,y,x])) cat(y,B_vyx[v,y,x],srv$error[2],'\n',file=conn);
                    }#y
                }#x
            }
            if (srv$output[3]){
                cat("#------------NUMBERS-AT-SIZE DATA-----------#\n",file=conn);
                cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=conn);
                cat("BY_SEX    #objective function fitting option\n",file=conn);
                cat("MULTINOMIAL #likelihood type\n",file=conn);
                cat(n_v[v],"     #number of years of data\n",file=conn);
                cat("???         #units\n",file=conn);
                cat(d$zc$n,"  #NUMBER OF SIZE BIN CUTPTS\n",file=conn);
                cat("#SIZE BIN CUTPTS (mm CW)\n",file=conn);																																	
                cat(d$zc$nms,"\n",file=conn);
                cat("#--------------\n",file=conn);
                cat(d$x$n*d$m$n*d$s$n,"   #number of factor combinations\n",file=conn);
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
                                for (r in 1:nrow(dp)){
                                    if (sum(dp[r,4+(1:d$z$n)],na.rm=TRUE)>0){
                                        cat(dp[r,4],50,file=conn);
                                        for (j in 4+(1:d$z$n)) cat(' ',dp[r,j],file=conn);
                                        cat('\n',file=conn)
                                    }
                                }#r
                            }
                        }#s
                    }#m
                }#x
            }
        } else {
            cat("#-----no survey catch data\n",file=conn);
        }
    }#v
}