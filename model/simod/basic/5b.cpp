

$PROB CEFTOZOLAN/TAZOBAKTAM

$CMT DEP1 PER1 DEP2 PER2

$PARAM
 X1 = 1;
 X2 = 1;
 X3 = 1;
 X4 = 1;
 X6 = 1;
 X7 = 1;
 X8 = 1;
 X9 = 1;
 WT = 70;
 WTCLC = 0.79;
 WTQC = 1;
 WTV1C = 0.87;
 WTV2C = 0.48;
 WTCLT = 0.65;
 WTQT = 0.75;
 WTV1T = 0.74;
 WTV2T = 0.83;

$MAIN
double CLC=X1*pow(WT/70,WTCLC);
double QC=X2*pow(WT/70,WTQC);
double V1C=X3*pow(WT/70,WTV1C);
double V2C=X4*pow(WT/70,WTV2C);
double CLT=X6*pow(WT/70,WTCLT);
double QT=X7*pow(WT/70,WTQT);
double V1T=X8*pow(WT/70,WTV1T);
double V2T=X9*pow(WT/70,WTV2T);

double K10 = CLC/V1C;
double K12 = QC/V1C;
double K21 = QC/V2C;
double K30 = CLT/V1T;
double K34 = QT/V1T;
double K43 = QT/V2T;

$DES
dxdt_DEP1= -K12*DEP1+K21*PER1-K10*DEP1;
dxdt_PER1=  K12*DEP1-K21*PER1 ;
dxdt_DEP2= -K34*DEP2+K43*PER2-K30*DEP2;
dxdt_PER2=  K34*DEP2-K43*PER2;

$TABLE
double  CEF=DEP1/V1C;
double TAZ=DEP2/V1T;

$CAPTURE @etas 1:LAST
CEF TAZ
 
