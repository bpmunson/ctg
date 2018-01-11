library(qvalue)

Var<-function(x) mean(x^2)-mean(x)^2 #scalar version
vVar<-function(x) apply(x^2,1,mean)-apply(x,1,mean)^2 #vector version

Cov<-function(x,y) mean(x*y)-mean(x)*mean(y)
vCov<-function(x,y) apply(t(x)*y,2,mean)-apply(x,1,mean)*mean(y) #x is a matrix and y is a vector

sqrtsum<-function(y) sqrt(sum(y^2))


# x1 is log2 frequencies for the 1st replicate of all timepts
# x2 is log2 frequencies for the 2nd replicate of all timepts
# ab1 is abundance thresholds for all 1st replicates
# ab2 is abundance thresholds for all 2nd replicates
fit_ac_fc<-function(x1,ab1,x2,ab2) { #badx is TRUE when x-value is bad

   er_ac<-1
   l<-0
   nx<-nrow(x1)


   good1<-t(t(x1)>ab1) # constructs as rows, timept for 1st replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
   good2<-t(t(x2)>ab2) # constructs as rows, timept for 2nd replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
   useless1<-apply(good1,1,sum)<2 # 1 = sum over rows--i.e., constructs.  Any construct that isn't above abundance threshold in 1st replicate in at least 2 timepoints has 1/TRUE in "useless" matrix
   useless2<-apply(good2,1,sum)<2

    # For each constructs that isn't above abundance thresholds in at least two timepoints for this replicate, this sets that construct's "good" values to FALSE for *all* timepoints in this replicate
   good1[useless1,]<-FALSE #remove singletons
   good2[useless2,]<-FALSE #remove singletons

    # Note: here apply(goodX,1,sum) is NOT uselessX because goodX was changed above
    #allbad is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
   allbad<-apply(good1,1,sum)<2 & apply(good2,1,sum)<2 #in this case I have nothing to use in either experiment

    # nt = number of timepoints
   lambda1<-rep(0,nt)
   lambda2<-rep(0,nt)
   ac1<-x1[,1] #just a guess # log2 frequencies for all constructs for first timepoint in this replicate
   ac2<-x2[,1] #just a guess
   fc<-rep(0,nx) # nx = number of constructs

   #BPM add, cache replicate fitnesses
   f1s<-rep(0,nx)
   f2s<-rep(0,nx)

    # for 1 to number of constructs
   for (i in 1:nx) {
       # if this construct doesn't have at least two timepoints above abundance threshold in at least one replicate, ignore it and move on
      if (allbad[i]) next #from now on there is at least one good experiment

    # apparently "f" stands for fitness--which is covariance (see below)--and # stands for replicate
    # v stands for variance and # " "
      f1<-0
      f2<-0
      v1<-0
      v2<-0

      g1<-good1[i,] # g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the first replicate
      if (sum(g1)>1) { #it's a good experiment # if there are at least two good timepoints for this construct in replicate 1
         mx1<-mean(x1[i,g1]) # get the mean of the log2 frequencies for all the good timepoints for this construct in replicate 1
          # time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
         mt1<-mean(time[g1]) # get mean of timepoints (e.g. mean number of days) for all the good timepoints for this construct in replicate 1
         v1<-Var(time[g1]) # get variance of timepoints " "
         f1<-Cov(x1[i,g1],time[g1]) # f1 = covariance of log2 frequencies for all the good timepoints for this construct in replicate 1 with the timepoints for those good timepoints
      }

      # do the exact same thing as above, but for replicate 2
      g2<-good2[i,]
      if (sum(g2)>1) { #it's a good experiment
         mx2<-mean(x2[i,g2])
         mt2<-mean(time[g2])
         v2<-Var(time[g2])
         f2<-Cov(x2[i,g2],time[g2])
      }

    # fc[i] is the combined fitness (across replicates) for construct i
      fc[i]<-(f1+f2)/(v1+v2) #the combined fitness from replicate 1+2
      #fc remains defined up to an additive constant

      #BPM add: cache replicate fitnesses
      f1s[i]<-f1/v1
      f2s[i]<-f2/v2

      if (sum(g1)>1) { # if there are at least two good timepoints for this construct in replicate 1
          # ac1 = mean of log2 freqs for this construct for good timepoints for rep 1 - (mean of timepts for good timepoints for this construct for rep 1)*combined fitness across replicates for this construct
         ac1[i]<-mx1-fc[i]*mt1
      }

    # same as above but for replicate 2
      if (sum(g2)>1) {
         ac2[i]<-mx2-fc[i]*mt2
      }
   }

    # ac is the initial condition (in log2 frequency) for construct c
    # this is normalizing ac1:
    # Roman's methods say "By definition, log2 relative frequencies satisfy the constraint
    # sum over c of (2^xc) = 1 at all times"
    # ac1 is the set of log2 frequencies for all constructs in replicate 1 at "initial conditions", so
    # ac1 must satisfy the above constraint.  IFF the above constraint is satisfied, -log2(1) = 0.
    # If the above constraint isn't satisfied, the value of -log2(sum over c of (2^ac)) is subtracted from
    # ac to ensure the constraint is satisfied.
   alpha<- -log2(sum(2^ac1))
   ac1<-ac1+alpha #enforce normalization at time=0, sum(2^ac)=1

    # same as above but for replicate 2
   alpha<- -log2(sum(2^ac2))
   ac2<-ac2+alpha #enforce normalization at time=0, sum(2^ac)=1

    # ok, the *expected* log2 frequency for construct x at time t is:
    # xc(t) = ac + fc*t - log2(sum over c of 2^(ac + fc*t))
    # The below code calls the term starting with "-log2" "lambda".
    # Note that lambda is being calculated separately for each combination of timepoint+replicate

    # for 1 to number of timepoints
   for (i in 1:nt) {
      lambda1[i]<- -log2(sum(2^(ac1+fc*time[i])))
      lambda2[i]<- -log2(sum(2^(ac2+fc*time[i])))
   } #these are initial estimates of lambda(t)

    # x1 is log2 frequencies for the 1st replicate of all timepts
   xfit1<-x1 #for size # I think this means that xfitX is being set to xX not because we're using *any* of the
    # xX values but just to initialize xfitX to the desired size (which is the same as the size of xX)

    # for 1 to number of timepoints
   for (j in 1:nt) {
       # the expected value of xc(t) at this timepoint t, calculated as a column for all constructs c
      xfit1[,j]<-ac1+fc*time[j]+lambda1[j]
   }

    # same as above but for replicate 2
   xfit2<-x2 #for size
   for (j in 1:nt) {
      xfit2[,j]<-ac2+fc*time[j]+lambda2[j]
   }


   sdfc<-rep(0.1,nx) #standard error of fc
   tstat<-rep(0,nx) # presumably t statistic
   df<-rep(0,nx) # df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant construct are good, -1 if just one is, -2 if neither are.  Inits to 0 for all
   p_t<-rep(1,nx) # p value from t test ... initialize to 1 (not significant) for everything

    # for 1 to number of constructs
   for (i in 1:nx) {
       # if this construct doesn't have enough "good" measurements, skip it
      if (allbad[i]) next


      g1<-good1[i,] # g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the first replicate
      g2<-good2[i,]

        # df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant construct are good, -1 if just one is, -2 if neither are
        # I suspect that "df" is "degrees of freedom" for each construct
      df[i]<-sum(g1)+sum(g2)-2

        # sqrtsum<-function(y) sqrt(sum(y^2))
        # outside of this function, this sdfc value is used only in the output of the construct file
        # I think that sdfc is "standard deviation of fc" for each construct.

        # So, Roman says:
        # std err of fc = sqrt of sum over t of (lower-case-epsilon for construct c, as a function of t)^2
        # divided by sqrt of (nc - 2)*sum over t of (t^2 - (mean of t)^2)
        # Note that lower-case-epsilon is Xc(t) - xc(t) = xX - xfitX
        # so xfitX - xX = - lower-case-epsilon ... but since it is being squared and then square-rooted, I
        # suppose the negation doesn't matter.
        # nc - 2 is the number of degrees of freedom
        # where nc = number of data points = 2*nt minus any number of points below the threshold
        # (note the description above seems to assume 2 replicates) ... seems to me this must mean
        # number of data points *for this construct c*, not total.
        # nt in above is number of timepoints, as here ...
        # Roman also says that tc [i.e., the t statistic for construct c] = fc /SE(fc)
        # and the internet tells me that SE(x) = SD(x)/sqrt(n) ... but maybe the sqrt(n) is just the most usual
        # sqrt of degrees of freedom, and could be something else in a more complex system.
        # So, the value being calculated directly below is SD(x), which is why it doesn't have the
        # nc - 2 term that is in the denominator of the SE(x) calculation (in the manuscript, Roman says that
        # fc's sd = sqrt(nc-2)*SE(fc), where nc-2 = degrees of freedom.  Farther below, after the calculation
        # of sdfc, we get the calculation of tc, which is fc/(sdfc/sqrt(df)), where the sdfc/sqrt(df) term
        # is the calculation of SE(fc).

    # result is vector with one std dev of fc for each construct
      sdfc[i]<-sqrtsum( c(xfit1[i,g1],xfit2[i,g2]) - c(x1[i,g1],x2[i,g2]) ) /
          sqrtsum( c(time[g1],time[g2]) - mean(c(time[g1],time[g2])) )

   }
   #find median sd
   has_sd<-df>0
   median_sd<-median(sdfc[has_sd])
   sdfc[!has_sd]<-median_sd #just so it isn't 0

   # for 1 to number of constructs
   for (i in 1:nx) {
       # don't try to calculate t statistic and p value for any fc that doesn't have a stderr
      if (!has_sd[i]) next
        # calc t statistic of fc
      tstat[i]<-fc[i]/(sdfc[i]/sqrt(df[i]))
        # pt is R function, presumably to do t-test :)  Result must be one-tailed, hence *2
      p_t[i]<-2*pt(-abs(tstat[i]),df=df[i]) #raw p-values from t-test
   }

    # ok, lfdr stands for "local fdr", and the lfdr function comes from the qvalue bioconductor package
    # that is installed way above.  The first input is the vector of p-values, and the second input is the
    # estimated proportion of true null p-values, where pi0.method is "the method for automatically choosing tuning
    # parameter in the estimation of π0, the proportion of true null hypotheses."
   lfdr_fc<-rep(1,nx) # nx = number of constructs; default value of lfdr is set to one for all of them
   l<-lfdr(p_t[has_sd],pi0.method="bootstrap")
   lfdr_fc[has_sd]<-l # for constructs that have an sd, the lfdr of the fc could be calculated and is now set
       # to its calculated vale instead of the default

    # ac1 is the initial condition (in log2 frequency) for each construct c for replicate 1
    # ac2 is the initial condition (in log2 frequency) for each construct c for replicate 2
    # fc is the fitness of each construct c (calculated across both replicates)
    # sdfc is the std deviation of the fitness of each construct c (calculated across both replicates)
    # p_t is the raw p value of the fc of each construct c (calculated across both replicates)
    # lfdr_fc is the local FDR of each construct (calculated across both replicates)
    # df is the degrees of freedom of each construct c (calculated across both replicates)
    # allbad is a boolean value for each construct c that is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
   #vl<-list(ac1,ac2,fc,sdfc,p_t,lfdr_fc,df,allbad)
   vl<-list(ac1,ac2,fc,sdfc,p_t,lfdr_fc,df,allbad,f1s,f2s)# BPM added return of replicate fitnesses
   return(vl)
}
