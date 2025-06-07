
### Function for single stage design
singleStageII <- function(p0,p1,alpha_star,pow_star,gamma,plot=T){

  theta <- (p1*(1-p0))/(p0*(1-p1));theta
  CPOW<-0
  m1 <- 0

  if(plot) plot(0,0,col=0,xlim=c(0,150),ylim=c(0,1))

  while(CPOW < pow_star){

    m1 <- m1 + 1
    n1 <- max(round(m1*gamma),1)
    zn <- (n1+m1)
    cond_pow <- rep(NA,zn)
    pow_zz <- rep(NA,zn)
    cond_alp <- rep(NA,zn)
    #AI <- rep(NA,zn)

    for(z in 1:zn){

      mp <- min(z,max(n1,m1))
      a <- seq(1,min(z,m1))
      alpha <- a
      cont <- T

      #st.a <- max(1,round(z*n1/(2*(n1+m1))-0.49))
      st.a <- 1
      st.a <- min(st.a,length(a))

      for(i in st.a:length(a)){
        am <- round((a[i]) + 1e-07);am
        if((am>=mp)|z==zn) alpha[i]<-1
        if(am<mp & cont==T) {
          x<-am:min(z,m1);x
          y<-z-x;y
          cd <- (x*(n1-y))/(y*(m1-x));cd;x;y;n1;m1
          ind <- which(cd>=1)
          if(length(ind)>0) alpha[i] <- sum(dhyper(x[ind],m1,n1,z))
          if(alpha[i]<alpha_star) cont <- F
        }
      }


      ## smallest a to satisfy alpha_star
      a.int <- which(alpha<alpha_star)
      #AI[z] <- a.int

      pow_z <- 0
      if(length(a.int)>0){
        a_z <- a[min(a.int)]
        ## calculating power
        pow_z <- sum(pmf(a_z:mp,z,m1,n1,theta))
      }

      if (plot) points(m1,pow_z)
      ## Marginal power
      gmi <- max(0,z-n1)
      gma <- min(z,m1)
      g_z <- sum(dbinom(gmi:gma,m1,p0)*(dbinom(z-(gmi:gma),n1,p1)))
      g_z0 <- sum(dbinom(gmi:gma,m1,p0)*(dbinom(z-(gmi:gma),n1,p0)))

      pow_zz[z] <- pow_z
      cond_pow[z] <- pow_z*g_z;
      cond_alp[z] <- 0
      if(length(a.int)>0) cond_alp[z] <- alpha[a_z]*g_z0

    }
    m1;n1;pow_zz;cond_pow
    ### saving paramers
    CPOW <- sum(cond_pow);CPOW;m1
    CALP <- sum(cond_alp)

    if(plot) points(m1,CPOW,pch=20,col=2,cex=2);abline(h=pow_star)
  }

  res  <- data.frame("p0"=p0,"p1"=p1,"pow_star"=pow_star,"alpha_star"=alpha_star,"gamma"=gamma,"m1"=m1,"n1"=n1,"ss"=(n1+m1),"Power"=CPOW,"Alpha"=CALP,"Theta"=theta)

}


## Function for two-stage design
randPhaseII <- function(p0,p1,alpha_star,pow_star,gamma){

  ### Setting up parameters
  theta <- (p1*(1-p0))/(p0*(1-p1))
  N1 <- 300
  N2 <- 300
  a_1 <- 0
  Ncond <- 0
  end <- 0

  Power <- matrix(NA,N1,N2)
  Alpha <- matrix(NA,N1,N2)
  PET <- matrix(NA,N1)
  stage1 <- matrix(NA,N1,N2)
  stage2 <- matrix(NA,N1,N2)

  start <-round((20/theta))

  start2 <- start
  ss1 <-max(1,start)
  pow.prev <- 0
  ss2.prev <- 0
  cond <-1
  plot(0,0,xlim=c(0,200),ylim=c(0,500),col=0)
  COUNT <- 0
  min.ss <- 999

  while(end!=1){
    COUNT<- COUNT+1

    ss2 <- max(1,start2)
    m1 <- ss1
    n1 <- max(round(m1*gamma),1)
    zn1 <- m1+n1
    force <- 0

    ncp<-0
    count <- 0


    while(Ncond < pow_star & force==0){
      count <- count+1

      m2 <- ss2
      n2 <- max(round(m2*gamma),1)
      zn2 <- m2+n2

      Cpow.z <- matrix(0,(zn1+1),(zn2+1))
      Calp.z <- matrix(0,(zn1+1),(zn2+1))
      pet.z  <- matrix(0,(zn1+1))

      all.des <- mapply(desParmT,z1=rep(0:zn1,zn2+1),z2=as.numeric(gl(zn2+1,zn1+1))-1,m1=m1,n1=n1,m2=m2,n2=n2,p0=p0,p1=p1,alpha_star=alpha_star)


      # Setting up next round
      Ncond <- sum(unlist(all.des[2,]))
      ncp <- c(ncp,Ncond)
      ## Failsafe for massively large studies
      if(ss2>250) Ncond<-1

      # Saving parameters
      PET[ss1] <- sum(unlist(unique(all.des[3,])))
      Power[ss1,ss2] <- Ncond
      Alpha[ss1,ss2] <- sum(unlist(all.des[1,]))
      stage1[ss1,ss2] <- n1+m1
      stage2[ss1,ss2] <- n1+m1+n2+m2

      ## setting exit conditions
      if((ss2 > (ss1*(3/gamma)))|(ss1+ss2)>400){
        force <- 1
      }

      # Setting up next round
      if(Ncond<pow_star&force==0) {

        incr <- int(ss2.prev[max(round(cond),1)],ss2,pow.prev[max(round(cond),1)],Ncond,pow_star);incr

        ### parameters for next increase
        ss2.prev <- c(ss2.prev,ss2)
        pow.prev <- c(pow.prev,Ncond)
        cond <- cond+0.33
        ss2 <- ss2 + incr
        #ss2 <- ss2 +1
      }
      Ncond;pow_star

    }

    ss.cond <- stage2[ss1,ss2];ss.cond;min.ss

    points(ss1,Ncond*400,pch=20,col=2);abline(h=pow_star*400,lty=2)
    points(ss1,stage1[ss1,ss2],pch=15,col=3)
    points(ss1,stage2[ss1,ss2]-stage1[ss1,ss2],pch=15,col=4)
    points(ss1,ss.cond,pch=15,col=6)

    if(Ncond>pow_star) min.ss <- min(min.ss,ss.cond)

    stop.con <- min(min.ss+50,c(min.ss*1.25))
    Ncond;ss1;ss2;ss.cond;start2;min.ss;stop.con

    if(ss.cond>stop.con) end <- 1

    # resetting
    Ncond <- 0
    pow.prev <- 0
    ss2.prev <- 0
    cond <- 1
    ss1 <- ss1+1
    start2 <- round(max(ss1*0.66,round(ss2-(5+50/ss1)*((1/gamma)+2*force))));start2;end;force
    force <- 0
  }

  ### The following code actually selects the design
  ss.cond <- (Power>pow_star)+0

  ### Optimum Sample Size
  mPET <- matrix(PET,N1,N2)
  optSS <- (mPET*stage1+(1-mPET)*stage2)*ss.cond

  optSS[-which(is.na(optSS))]
  SOss <- min(optSS[-which(is.na(optSS)|optSS==0)])
  SO.id <- which(optSS==SOss)[1]
  #####

  ### MinMax Sample Size
  tot.SS <-  c(stage2*ss.cond+0)
  tot.SS <- min(tot.SS[-which(is.na(tot.SS)|tot.SS==0)])

  SSpow <- min(Power[which(stage2*ss.cond==tot.SS)])
  SS.id <- which(Power==SSpow)[1]
  #####
  n1 <- round(stage1[c(SS.id,SO.id)]*(gamma/(gamma+1)))
  m1 <- round(stage1[c(SS.id,SO.id)]*(1/(gamma+1)))
  n2 <- round(stage2[c(SS.id,SO.id)]*(gamma/(gamma+1))) - n1
  m2 <- round(stage2[c(SS.id,SO.id)]*(1/(gamma+1))) - m1

  ss <- n1+n2+m1+m2
  pow <- Power[c(SS.id,SO.id)]
  alp <- Alpha[c(SS.id,SO.id)]
  pet <- mPET[c(SS.id,SO.id)]
  des <- c("MM","OPT")

  res <- data.frame(p0=p0,p1=p1,pow_star=pow_star,alpha_star=alpha_star,gamma=gamma,n1=n1,m1=m1,n2=n2,m2=m2,ss=ss,Power=round(pow,4),Alpha=round(alp,4),PET=round(pet,3),DES=des)
}

## Anxilliary Functions
### Probability mass function
pmf <- function(x,z,m,n,theta){

  y <- z - x
  mn <- max(0,(z-n))
  mp <- min(z,m)

  num <- choose(m,x)*choose(n,y)*(theta^x)
  denom <- sum(choose(m, mn:mp)*choose(n,z-(mn:mp))*(theta^(mn:mp)))
  res <-num/denom

  inf.id <- which(res=="Inf"|res=="NaN"|is.na(res))
  if(length(inf.id)>0) res[inf.id] <- 1
  res

}


# Design parameters
desParmT <- function(m1,n1,m2,n2,z1,z2,p0,p1,alpha_star){


  theta <- (p1*(1-p0))/(p0*(1-p1))
  mn1 <- max(0,(z1-n1))
  mn2 <- max(0,(z2-n2))
  mp1 <- min(z1,m1)
  mp2 <- min(z2,m2)


  #### This is what we use if it's unequal allocation
  a <- c(1,0)
  x_1 <- mn1:mp1
  y_1 <- z1 - (x_1)
  x_2 <- mn2:mp2
  y_2 <- z2 - (x_2)


  X <- t(matrix(x_1,length(x_1),length(x_2)))+x_2
  N <- n1+n2
  M <- m1+m2
  Y <- (z1+z2) - X


  theta1_hat <- abs( (x_1*(n1-y_1)) / (y_1*(m1-x_1)) )
  theta_hat <-  abs( (X*(N-Y)) / (Y*(M-X)) )

  if(is.na(theta1_hat[length(theta1_hat)])|theta1_hat[length(theta1_hat)]==Inf){
    theta1_hat[length(theta1_hat)]<-1000
  }

  if(is.na(theta_hat[nrow(theta_hat),ncol(theta_hat)])|theta_hat[nrow(theta_hat),ncol(theta_hat)]==Inf){
    theta_hat[nrow(theta_hat),ncol(theta_hat)]<-1000
  }


  ### Stage 1 condition for continuation (set at THETA>=1 for now)
  ind.s1 <- theta1_hat>=a[1]

  ff <- t(dhyper((mn1:mp1)[ind.s1],m1,n1,z1)%*%t(dhyper(mn2:mp2,m2,n2,z2)))

  TH <- sort(unique(c(theta_hat)))
  alpha <- rep(1,length(unique(c(theta_hat))))

  ### Cutting down the number of evaluations that need to be carried out
  a.st <- 1
  THcon <- which(TH>=1)
  THcon

  if(length(THcon)>0){
    a.st <- min(THcon)
    cont <- T
    ### Choosing 'a'
    for(a.i in a.st:length(TH)){
      if(cont){
        ind <- ( (theta_hat>=TH[a.i])+0 )[,ind.s1]
        if(sum(ind)==0|TH[a.i]<=1) alpha[a.i]<-1
        if(sum(ind)>0&TH[a.i]>1) alpha[a.i] <- sum(ff*ind)
        if(alpha[a.i]<alpha_star) cont <-F
      }
    }
  }


  a[2]<-Inf # An inordinately large number
  a.id <- which(alpha<alpha_star)
  if(length(a.id)>0) a[2] <- TH[min(a.id)]
  ####

  ### Power for chosen a
  if(a[2]!=Inf){
    ind <- ( (theta_hat>=a[2])+0 )[,ind.s1] ;ind
    ffp <- t(pmf((mn1:mp1)[ind.s1],z1,m1,n1,theta)%*%t(pmf(mn2:mp2,z2,m2,n2,theta)))
  }
  ##########


  ##### Calculation conditional probabilities
  ## marginal probabilities of observing z (H1)
  g1 <- sum(dbinom(mn1:mp1,m1,p0)*dbinom(z1-(mn1:mp1),n1,p1))
  g2 <- sum(dbinom(mn2:mp2,m2,p0)*dbinom(z2-(mn2:mp2),n2,p1))

  ## marginal probabilities of observing z (H0)
  g01 <- sum(dbinom(mn1:mp1,m1,p0)*dbinom(z1-(mn1:mp1),n1,p0))
  g02 <- sum(dbinom(mn2:mp2,m2,p0)*dbinom(z2-(mn2:mp2),n2,p0))

  ## 1-b  for every combination of 'z'
  Cpow.z <- 0
  Calp.z <- 0
  if(a[2]!=Inf){
    Cpow.z <- sum(ffp*ind)*g1*g2
    Calp.z <- sum(ff*ind)*g01*g02}

  # Probability of Early Termination	(under null hypothesis)
  pet.z <- sum(dhyper((mn1:mp1)[!ind.s1],m1,n1,z1))*g01

  res<-data.frame("alp"=Calp.z,"pow"=Cpow.z,"z"=pet.z)
  res

}


# Calculating the increment between sample sizes
int <- function(ssn1,ssn2,pow1,pow2,pow_star){
  grad <- (pow2-pow1)/(ssn2-ssn1)
  targ <- (pow_star - pow2)/grad
  inc <- max(1,round(targ))
  inc <- min(inc,15)
  inc
}

# Estimating the allocation ratio
binAlloc <- function(p0,p1){
  A <- sqrt( (p0*(1-p0))/(p1*(1-p1)) )
  gam <- (1+A)^(-1)
  gam/(1-gam)
}
