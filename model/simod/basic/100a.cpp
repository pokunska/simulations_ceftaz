

$PROB  
# Model jednokompartmentowy, 
# Podanie dożylne CEFT i TAZ, 
# Populacja pediatryczna ze skalowaniem allometrycznym
# WT na CL i V dla ceftolozanu i tazobaktamu

$CMT DEPOT1 PER1 DEPOT2

$PARAM @annotated
TVCLC   : 6.37    : Typical clearence of cef for 70 kg subjects (L/h)  
TVVCC   : 30.2    : Typical central volume of cef for 70 kg subjects (L) 
TVCLT   : 15.5    : Typical clearence of taz for 70 kg subjects (L/h) 
TVVCT   : 55.4    : Typical central volume of taz for 70 kg subjects (L)
TVSD    : 0.812   : SDCLT:SDv1T
WT      : 70      : Body mass of a typical subject, kg #znajdź równanie w Covariate analysis w pracy z której korzystasz i zobacz na jaką masę ciała jest ustawione to równanie
WTCLC    : 0.79    : Allometric exponent for CL cef
WTVCC    : 0.87    : Allometric exponent for VC cef
WTCLT    : 0.65    : Allometric exponent for CL taz
WTVCT    : 0.74    : Allometric exponent for VC taz



$MAIN
double V1C=TVVCC*pow(WT/70,WTVCC)*exp(EV1C) ;      
double CLC=TVCLC*pow(WT/70,WTCLC)*exp(ECL1C) ;
double V1T=TVVCT*pow(WT/70,WTVCT)*exp(EV1T)     ;  
double CLT=TVCLT*pow(WT/70,WTCLT)*exp(TVSD*EV1T) ;


$DES

double AC1=DEPOT1/V1C;
double AT1=DEPOT2/V1T;

dxdt_DEPOT1=  - CLC * AC1;
dxdt_PER1=  0 ;
dxdt_DEPOT2=  - CLT * AT1;

$OMEGA @annotated @block @correlation
EV1C: 0.5       : ETA on V1C
ECL1C: 0.1  0.257        : ETA on CLC
EV1T: 0.1  0.1  0.507   : ETA on V1T

$TABLE

double  CEF=DEPOT1/V1C;
double TAZ=DEPOT2/V1T;

$CAPTURE @etas 1:LAST
CEF TAZ
 
