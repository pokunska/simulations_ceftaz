data {
  int<lower=1> nt;
  int<lower=1> nObsC;
  int<lower=1> nObsT;
  int<lower=1> nSubjects;
  int nIIV;
  array[nObsC] int<lower=1> iObsC;
  array[nObsT] int<lower=1> iObsT;
  array[nSubjects] int<lower=1> start;
  array[nSubjects] int<lower=1> end;
  array[nt] int<lower=1> cmt;
  array[nt] int<lower=1> cmtC;
  array[nt] int<lower=1> cmtT;
  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;
  vector<lower=0>[nObsC] cObsC;
  vector<lower=0>[nObsT] cObsT;
  int<lower=0, upper=1> runestimation; //   a switch to evaluate the likelihood
}

transformed data {
  vector[nObsC] logCObsC = log(cObsC);
  vector[nObsT] logCObsT = log(cObsT);
  int nTheta = 5;   // explanation ntheta = 8???
  int nCmt = 3;   // explanation ncmt = 2???
  array[nSubjects] int nti;
  array[nCmt] real biovar;
  array[nCmt] real tlag;

  for (i in 1 : nSubjects) nti[i] = end[i] - start[i] + 1;

  for (i in 1 : nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters {

  real<lower=0, upper=500> CLCHat;  
  real<lower=0, upper=500> QCHat;
  real<lower=0, upper=3500> V1CHat;
  real<lower=0, upper=3500> V2CHat;
  real<lower=0, upper=500> CLTHat;
  real<lower=0, upper=500> QTHat;
  real<lower=0, upper=3500> V1THat;
  real<lower=0, upper=3500> V2THat;
  vector<lower=0>[2] sigma;
  vector<lower=3>[2] nu; // normality constant

  // Inter-Individual variability
  vector<lower=0.01, upper=2>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
  cholesky_factor_corr[nIIV] L;   
}

transformed parameters {

  vector<lower=0>[nIIV] thetaHat;
  matrix<lower=0>[nSubjects, nIIV] etaM; // variable required for Matt's trick
  array[nTheta] real<lower=0> thetaC;
  array[nTheta] real<lower=0> thetaT;
  matrix<lower=0>[nCmt, nt] xT;
  matrix<lower=0>[nCmt, nt] xC;
  row_vector<lower=0>[nt] cHatC;
  row_vector<lower=0>[nt] cHatT;
  row_vector<lower=0>[nObsC] cHatObsC;
  row_vector<lower=0>[nObsT] cHatObsT;
  
  thetaHat[1] = CLCHat;
  thetaHat[2] = QCHat;
  thetaHat[3] = V1CHat;
  thetaHat[4] = V2CHat;
  thetaHat[5] = CLTHat;
  thetaHat[6] = QTHat;
  thetaHat[7] = V1THat;
  thetaHat[8] = V2THat;
  
 // Matt's trick to use unit scale 
  etaM =  exp(diag_pre_multiply(omega, L * etaStd))'; 
  
  for(j in 1:nSubjects)
  {
    thetaC[1] = thetaHat[1] * etaM[j, 1] ; // CL
    thetaC[2] = thetaHat[2] * etaM[j, 2] ; // Q
    thetaC[3] = thetaHat[3] * etaM[j, 3] ; // V1
    thetaC[4] = thetaHat[4] * etaM[j, 4] ; // V2
    thetaC[5] = 0; // ka
    thetaT[1] = thetaHat[5] * etaM[j, 5] ; // CL
    thetaT[2] = thetaHat[6] * etaM[j, 6] ; // Q
    thetaT[3] = thetaHat[7] * etaM[j, 7] ; // V1
    thetaT[4] = thetaHat[8] * etaM[j, 8] ; // V2
    thetaT[5] = 0; // ka
    
    //amount 
    
    xC[,start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]], 
                                       amt[start[j]:end[j]],
                                       rate[start[j]:end[j]],
                                       ii[start[j]:end[j]],
                                       evid[start[j]:end[j]],
                                       cmtC[start[j]:end[j]],
                                       addl[start[j]:end[j]],
                                       ss[start[j]:end[j]],
                                       thetaC, biovar, tlag);
    xT[,start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]], 
                                       amt[start[j]:end[j]],
                                       rate[start[j]:end[j]],
                                       ii[start[j]:end[j]],
                                       evid[start[j]:end[j]],
                                       cmtT[start[j]:end[j]],
                                       addl[start[j]:end[j]],
                                       ss[start[j]:end[j]],
                                       thetaT, biovar, tlag); 
                                       
    //concentrations 
    
    cHatC[start[j]:end[j]] = xC[2,start[j]:end[j]] ./ thetaC[3]; // divide by V1
    cHatT[start[j]:end[j]] = xT[2,start[j]:end[j]] ./ thetaT[3]; // divide by V1
  }

  cHatObsC  = cHatC[iObsC];
  cHatObsT  = cHatT[iObsT];
}

model{
  //Informative Priors
      
  CLCHat ~ lognormal(log(5.882),0.25);
  QCHat  ~ lognormal(log(2.545),0.25);
  V1CHat ~ lognormal(log(10.64),0.25);
  V2CHat ~ lognormal(log(4.227),0.25);
  CLTHat ~ lognormal(log(20.8),0.25);
  QTHat  ~ lognormal(log(4.06),0.25);
  V1THat ~ lognormal(log(12.9),0.25);
  V2THat ~ lognormal(log(5.06),0.25);
  L~lkj_corr_cholesky(10);

  nu ~ gamma(2,0.1);

 // Inter-individual variability (see transformed parameters block
 // for translation to PK parameters)
  to_vector(etaStd) ~ normal(0, 1);
  omega ~ lognormal(log(0.4),0.25);
  sigma ~ lognormal(log(0.10), 0.25);

  if(runestimation==1){
    logCObsC ~ student_t(nu[1],log(cHatObsC), sigma[1]);
    logCObsT ~ student_t(nu[2],log(cHatObsT), sigma[2]);
  }
}

generated quantities{

  matrix[nIIV, 900] etaStdPred;
  matrix<lower=0>[900, nIIV] etaPredM;
  array[nTheta, 900] real<lower=0> thetaPredC;
  array[nTheta, 900] real<lower=0> thetaPredT;
    

    for(i in 1:900){
      for(j in 1:nIIV){ 
        etaStdPred[j, i] = normal_rng(0, 1);
      }
    }

    etaPredM = exp(diag_pre_multiply(omega, L * etaStdPred))';

    for(j in 1:900){
     
    thetaPredC[1,j] = thetaHat[1] * etaPredM[j, 1] ; // CL
    thetaPredC[2,j] = thetaHat[2] * etaPredM[j, 2] ; // Q
    thetaPredC[3,j] = thetaHat[3] * etaPredM[j, 3] ; // V1
    thetaPredC[4,j] = thetaHat[4] * etaPredM[j, 4] ; // V2
    thetaPredC[5,j] = 0; // ka
    thetaPredT[1,j] = thetaHat[5] * etaPredM[j, 5] ; // CL
    thetaPredT[2,j] = thetaHat[6] * etaPredM[j, 6] ; // Q
    thetaPredT[3,j] = thetaHat[7] * etaPredM[j, 7] ; // V1
    thetaPredT[4,j] = thetaHat[8] * etaPredM[j, 8] ; // V2
    thetaPredT[5,j] = 0; // ka
    
  }
 
} 

