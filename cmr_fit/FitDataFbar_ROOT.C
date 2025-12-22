//
// FitData.cc (v.3.0) release: May, 10, 2003 
// Laboratory for Nuclear Science, Massachusetts Institute of Technology
//   Sean Stave     <stave@lns.mit.edu>
//   Kenneth Jensen <sanctity@disenchanted.mit.edu>
//   Itaru Nakagawa <itaru@lns.mit.edu>
//

//S. Stave 12-08-04: Modifying FitData.cc to run in ROOT
//  Must run: gSystem->Load("libMinuit.so");
//S. Stave 1-3-05: Making changes so that the Fbars are fit
//S. Stave 3-28-05: Dumping original unfit mpoles and fbar
//S. Stave 10-17-05: Adding the rest of the I=3/2 and I=1/2 multipoles to the database

#include "Riostream.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TComplex.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TChain.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "FitDataFbar_ROOT.h"


//Needed for calculation of the range of the correlation coeff
Double_t Z(Double_t r){
  return 0.5*log((1.0+r)/(1.0-r));
}
Double_t rho(Double_t mu){
  return (exp(2.0*mu)-1.0)/(exp(2.0*mu)+1.0);
}

//
// Function name : Usage(char *argv[])
//
// Description   : show usage
// Input         : char *argv[]
// Return        : 
//
int
Usage(char *argv[]){

  cout << "\n Usage:\n  " << argv[0] << " [-hxcLrsdOQL] [-m <int>] [-e <int>] " 
          "[-C <int>] [-f <file>] \n\t  [-W <W>] [-q <Q2>]" << endl;
  cout << "\n Description: " << endl;
  cout << "\t Fit data with E1+, M1+, L1+, (M1-, L1-, L0+, E0+)" << endl;
  cout << "\t as free parameters. " << endl;
  cout << "\n Options:" << endl;
  cout << "\t -C <int>  confidence level of parameter errors. int=0,1,2,3"<< endl; 
  cout << "\t            0: total chi2=chi2_min+1, 1[def]: 1-sigma/parameter"<<endl; 
  cout << "\t -e <int>  1:ExpErr^2 = dSta^2" << endl;  
  cout << "\t    \t   2:ExpErr^2 = dSys^2" << endl;
  cout << "\t    \t   3:ExpErr^2 = dSta^2+dSys^2" << endl;
  cout << "\t    \t   7:ExpErr^2 = dSta^2+dSys^2+dMod^2" << endl;
  cout << "\t    \t   [def]:dSta^2+dSys^2"      <<endl;
  cout << "\t -q <Q2>   output the best fit result at Q2=<Q2> (GeV/c)^2" << endl;
  cout << "\t -W <W>    output the best fit result at W=<W> (MeV)" << endl;
  //  cout << "\t -O \t   input data file in old format" << endl;
  cout << "\t -f <file> input file    [def]:FitData.inp" << endl;
  //  cout << "\t -Q \t   run MINUIT in quiet mode [def]:off" << endl;
  cout << "\t -1 \t   run MINUIT only once. Calculate and print out XS for each data" << endl;
  cout << "\t -v \t   run MINUIT in verbose mode [def]:off" << endl;
  cout << "\t -m <int>  fitting mode (Use -L to list all modes)" << endl;
  cout << "\t    \t    111: -1 : Dump out original model prediction. No fit." << endl;
  cout << "\t    \t    112: 31 : par(3)=E1+,M1+,L1+  [def] " << endl;
  cout << "\t    \t    124: 56 : par(5)=E1+,M1+,L1+,M1-,L1-" << endl;
  cout << "\t    \t    127: 71 : par(7)=E1+,M1+,L1+,M1-,L1-,E0+,L0+" << endl;
  //  Following option was substituted with -m 111 option 
  //  cout << "\t -M \t   output unfitted MAID results "<< endl;
  cout << "\t -S \t   par error weighted by sqrt(chi2)" << endl;
  cout << "\t -P \t   output parameter code and exit " << endl;
  cout << "\t -r \t   output relative parameters (input for MAID homepage)" 
       << endl;
  cout << "\t -c \t   output chi2 of each observable" << endl;
  cout << "\t -2 \t   output chi2 of each observable(more detail)" << endl;
  cout << "\t -d \t   output in simple data mode (1)" << endl;
  cout << "\t    \t    Mode,Q2,W,E1+,dE1+,M1+,dM1+,S1+,dS1+,EMR,dEMR,CMR,dCMR," << endl;
  cout << "\t    \t    chi2,nDoF,ch2DoF" << endl;
  cout << "\t -X \t   output in extended format" << endl;
  cout << "\t -D \t   output in \"RfunK -w <wfile>\" format" << endl;
  cout << "\t -L \t   show list of all available modes" << endl;
  cout << "\t -p \t   use phi dependence " << endl;
  cout << "\t -s \t   show observable flags " << endl;
  cout << "\t -h \t   show this help    " << endl;
  cout << "\t -x \t   show examples  " << endl;
  cout << endl;
  exit(0);


}


//
// Function name : Example(char *argv[])
//
// Description   : show examples
// Input         : char *argv[]
// Return        : 
//
int 
Example(char *argv[]){

  cout << "\n Exapmle: " << endl;
  cout << "  o Perform 7 parameter fit with FitData.Y1998.inp data file. " << endl;
  cout << "\t" << argv[0] << " -m 127 -f FitData.Y1998.inp" << endl << endl;
  cout << "  o Reproduce Kamalov Fitting Results." << endl;
  cout << "\t" << argv[0] << " -f FitData.Kamalov.inp" << endl << endl;
  cout << "  o Following Executions Cause Same Results." << endl;
  cout << "\t" << argv[0] << " -m 127" << endl;
  cout << "\t" << argv[0] << " -m 71" << endl << endl;
  cout << "  o Dump out fitting result in RfunK weight format and recalc observables." << endl;
  cout << "\t" << argv[0] << " -m 63 -D > weight.dat " << endl ;
  cout << "\t" << "RfunK -w weight.dat " << endl;
  cout << endl;

  exit(0);


}

//
// Function name : ListMode()
//
// Description   : show the list of available modes
// Input         : 
// Return        : 
//
void
ListMode(){

  cout << endl;
  cout << "\t  Mode :ParCode" << endl;
  cout << "\t   111 :  -1 : Dump out original model prediction. No fit." << endl;
  cout << "\t   112 :  31 : par(3)=E1+,M1+,L1+  [def] " << endl;
  cout << "\t   113 :  41 : par(4)=E1+,M1+,L1+,L0+ " << endl;
  cout << "\t   114 :  42 : par(4)=E1+,M1+,L1+,E0+ " << endl;
  cout << "\t   115 :  51 : par(5)=E1+,M1+,L1+,E0+,L0+" << endl;
  cout << "\t   116 :  43 : par(4)=E1+,M1+,L1+,L1-" << endl;
  cout << "\t   117 :  52 : par(5)=E1+,M1+,L1+,L1-,L0+" << endl;
  cout << "\t   118 :  53 : par(5)=E1+,M1+,L1+,L1-,E0+" << endl;
  cout << "\t   119 :  61 : par(6)=E1+,M1+,L1+,L1-,E0+,L0+" << endl;
  cout << "\t   120 :  44 : par(4)=E1+,M1+,L1+,M1-" << endl;
  cout << "\t   121 :  54 : par(5)=E1+,M1+,L1+,M1-,L0+" << endl;
  cout << "\t   122 :  55 : par(5)=E1+,M1+,L1+,M1-,E0+" << endl;
  cout << "\t   123 :  62 : par(6)=E1+,M1+,L1+,M1-,E0+,L0+" << endl;
  cout << "\t   124 :  56 : par(5)=E1+,M1+,L1+,M1-,L1-" << endl;
  cout << "\t   125 :  63 : par(6)=E1+,M1+,L1+,M1-,L1-,L0+" << endl;
  cout << "\t   126 :  64 : par(6)=E1+,M1+,L1+,M1-,L1-,E0+" << endl;
  cout << "\t   127 :  71 : par(7)=E1+,M1+,L1+,M1-,L1-,E0+,L0+" << endl;
  cout << endl;
  exit(0);


  return ;

}


//
// Function name : 
//
// Description   : 
// Input         : 
// Return        : 
//
int 
getOutputKinema(int NDATA){

  int Q2_FLAG = 0;
  int W_FLAG = 0;

  for ( int i=0 ; i<NDATA ; ++i) {
    if ( q2[0]-q2[i] ) ++Q2_FLAG ;
    if ( w[0]-w[i] ) ++W_FLAG ;
  }

  Q2 = Q2_FLAG ? DEF_Q2 : q2[0] ;
  W  = W_FLAG  ? DEF_W  :  w[0] ;

  return 0;
}

//
// Function name : main(int argc, char *argv[]))
//
// Description   : handle options and call main functions
// Input         : argc, argv
// Return        : 
//
int 
//go(int argc, char *argv[]) {
go(int Mode=112, int ParOut=1, int QuietMode=1,int CnfLevel=1,int mod=0,int PhiDep=0,int p12_tmp=-1) {
  p12=p12_tmp;
  //  int Mode       = 112;
  int ShowChiSQR = 0;
  //  int ParOut     = 0;
  //  int QuietMode  = 1;
  //Default out
  int DataOut    = 0;
  //Custom output
  //  int DataOut    = 4;
  //Custom output for other set of mpoles
  //  int DataOut    = 5;

  //  int CnfLevel   = 1;
  //  int ExpErr     = 3;  // stat + sys
  int ExpErr     = 1;//stat only
  int RtnPcode   = 0;
  //int Ch2Weight  = 1;
  int Ch2Weight  = 0;
  int SingleExe  = 0;
  DumpChi2 = 1;

  //Current data format (extra 3/2 and 1/2 vars)
  OldData=0;
  //Unused format
  //OldData=1;
  //Old format (with CGLNF and S&P mpoles)
  //  OldData=2;

  if (mod==0)
    DATA_FILE="FitData_Fbar.inp";
  if (mod==1)
    DATA_FILE="FitData_Fbar.aznu00.Y2000.Bates.inp";
  if (mod==2)
    DATA_FILE="FitData_Fbar.maid00.Y2000.Bates.inp";
  if (mod==3)
    DATA_FILE="FitData_Fbar.sl2000.Y2000.Bates.inp";
  if (mod==4)
    DATA_FILE="FitData_Fbar.M2003.Y2000.Bates.inp";
  if (mod==41)
    DATA_FILE="FitData_Fbar.M2003.Y2000.Bates.noALT.inp";
  if (mod==5)
    DATA_FILE="FitData_Fbar.DMT.Y2000.Bates.inp";
  if (mod==104)
    DATA_FILE="CLAS.M2003.inp";
  if (mod==204)
    //    DATA_FILE="Mainz.M2003.inp";
    DATA_FILE="FitData_Fbar.M2003.MainzQ06.inp";
  if (mod==304)
    DATA_FILE="Mainz_phi.M2003.inp";
  if (mod==404)
    //    DATA_FILE="Mainz_phi.M2003online.inp";
    DATA_FILE="Mainz_phi.082605.M2003online.inp";
  if (mod==504)
    DATA_FILE="M2003.CLAS126_pseudo.inp";

  if (mod==999)
    DATA_FILE="FitData.tmp";



  /*
  int opt;
   while (EOF != (opt = getopt(argc, argv, "C:c12e:f:hm:dDXrsSxvOLPpQM?q:W:"))) {
    switch (opt) {
    case 'W':
      DEF_W = atof(optarg);
      break;
    case 'q':
      DEF_Q2 = atof(optarg);
      break;
    case 'C':
      CnfLevel = atoi(optarg);
      break;
    case 'O':
      OldData = 1;
      break;
    case '1':
      SingleExe = 1;
      break;
    case 'c':
      ShowChiSQR = 1;
      break;
    case '2':
      DumpChi2 = 1;
      break;
    case 'e':
      ExpErr = atoi(optarg);
      break;
    case 'f':
      DATA_FILE = optarg;
      break;
    case 'm':
      Mode = atoi(optarg);
      break;
    case '?':
      RtnPcode = 1;
      break;
    case 'x':
      Example(argv);
      break;
    case 'S':
      Ch2Weight = 1;
      break;
    case 'r':
      ParOut = 1;
      break;
    case 'd':
      DataOut = 1;
      break;
    case 'D':
      DataOut = 2;
      break;
    case 'X':
      DataOut = 3;
      break;
    case 's':
      ShowObsFlag();
      break;
    case 'L':
      ListMode();
      break;
    case 'P':
      RtnPcode = 1;
      break;
    case 'p':
      PhiDep = 1;
      break;
    case 'Q':
      QuietMode = 1;
      break;
    case 'v':
      QuietMode = 0;
      break;
    case 'M':
      MAID_OUT = 1;
      break;
    case 'h':
    case '*':
      Usage(argv);
    }

   }
  */   
   fit_data(Mode, ShowChiSQR, ParOut, QuietMode, DataOut, 
	    CnfLevel, ExpErr, RtnPcode, Ch2Weight, SingleExe, PhiDep);
  
  return 0;
} // end-main()

//added by H.Atac
void FitDataFbar_ROOT(){
go(112,1,1,1,999,1);
}

//
// Function name : fcn
//
// Description   : handle options and call main functions
// Input         : argc, argv
// Return        : 
//
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
         Int_t iflag ) {

  npar*=1;
  gin[0]*=1.0;
  iflag*=1;
  Double_t md, mf, mpion=0, mnucl=0;
  Double_t Fb1=0, Fb2=0, Fb3=0, Fb4=0, Fb5=0, Fb6=0;
  Double_t L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, 
    dL1=0, dL2=0, dL3=0, dL4=0, dL5=0, dL6=0, dL7=0,
    dL8=0, dL9=0, dL10=0,
    rlt[NUM_DATA],rt[NUM_DATA],rl[NUM_DATA],rtt[NUM_DATA],r0[NUM_DATA],
    alt[NUM_DATA],plt[NUM_DATA],Pxe[NUM_DATA],Py[NUM_DATA],Pze[NUM_DATA],
    rltp[NUM_DATA],rlt0y[NUM_DATA],rtlp_xp0[NUM_DATA],rttp_zp0[NUM_DATA], E_pi_cm[NUM_DATA],
    p_pi_cm[NUM_DATA], k_gamma_cm[NUM_DATA], Ppikgamma[NUM_DATA], 
    QOmegacm[NUM_DATA], qcm[NUM_DATA], Qqcm[NUM_DATA];
  Double_t Sig0[NUM_DATA], SigLT[NUM_DATA], SigLTp[NUM_DATA], SigTT[NUM_DATA];
  TComplex dF1[NUM_DATA], dF2[NUM_DATA], dF3[NUM_DATA], dF4[NUM_DATA],
    dF5[NUM_DATA], dF6[NUM_DATA]; 

  /*
  par[0] = 0.855;
  par[1] = 1.00377;
  par[3] = 0.97496;
  */



//    printf("Mpoles: %f %f %f\n",e1p[21].Im(),m1p[21].Im(),l1p[21].Im());
//    printf("Mpoles: %f %f %f\n",e1p3[6].Im(),m1p3[6].Im(),l1p3[6].Im());

  L1 = par[0];
  L2 = par[1];
  L3 = par[2];
  //  L1 = par[0]/e1p[21].Im();
  //  L2 = par[1]/m1p[21].Im();
  //  L3 = par[2]/l1p[21].Im();
  //  L1 = par[0]/3.0;
  //  L2 = par[1]/3.0;
  //  L3 = par[2]/3.0;
  L4 = par[3];
  L5 = par[4];
  L6 = par[5];
  L7 = par[6];

  Fb1=par[7];
  Fb2=par[8];
  Fb3=par[9];
  Fb4=par[10];
  Fb5=par[11];
  Fb6=par[12];

  L8 = par[13];
  L9 = par[14];
  L10 = par[15];

  /* Since the CGLN F's are the full values and the lambdas are the
     resonant piece only, a conversion must be done to see how the CGLN's are
     actually modified.  Part of the term needs to be subtracted off or
     added.*/
  Int_t i;
  for (i=0; i<NUM_DATA; i++) {
  //pi0,p
  if (obs[i]<100){
    dL1 = (L1 - 1.0)*2.0/3.0;
    dL2 = (L2 - 1.0)*2.0/3.0;
    dL3 = (L3 - 1.0)*2.0/3.0;
    dL4 = L4 - 1;
    dL5 = L5 - 1;
    dL6 = L6 - 1;
    dL7 = L7 - 1;
    dL8 = L8 - 1;
    dL9 = L9 - 1;
    dL10 = L10 - 1;
  }
  //pi+,n
  if (obs[i]>100){
    dL1 = (L1 - 1)*(-sqrt(2.0)/3.0);
    dL2 = (L2 - 1)*(-sqrt(2.0)/3.0);
    dL3 = (L3 - 1)*(-sqrt(2.0)/3.0);
    dL4 = L4 - 1;
    dL5 = L5 - 1;
    dL6 = L6 - 1;
    dL7 = L7 - 1;
    dL8 = sqrt(2.0)*(L8 - 1);
    dL9 = sqrt(2.0)*(L9 - 1);
    dL10 = sqrt(2.0)*(L10 - 1);
  }
    if (OldData==1) {
      dF1[i] = f1[i] + dL6*e0p[i] + (dL2*m1p[i] + dL1*e1p[i])*3.*cos(t_pi_rad[i]);
      dF2[i] = f2[i] + 2.*dL2*m1p[i] + dL4*(m1m0[i]+m1m1[i]/3.);
      dF3[i] = f3[i] + 3.*(dL1*e1p[i] - dL2*m1p[i]);
      dF4[i] = f4[i];
      dF5[i] = f5[i] + dL7*l0p[i] + 6.*cos(t_pi_rad[i])*dL3*l1p[i];
      dF6[i] = f6[i] - 2.*dL3*l1p[i] + dL5*(l1m0[i]+l1m1[i]/3.);
    } else {
      dF1[i] = Fb1*(f1[i] - (e0p[i] + (m1p[i] + e1p[i])*3.*cos(t_pi_rad[i]))) + 
	(e0p[i] + (m1p[i] + e1p[i])*3.*cos(t_pi_rad[i])) +
	dL6*e0p[i] + (dL2*m1p3[i] + dL1*e1p3[i])*3.*cos(t_pi_rad[i])
	+ (dL9*m1p1[i] + dL8*e1p1[i])*3.0*cos(t_pi_rad[i]);
      dF2[i] = Fb2*(f2[i] - (m1m[i]+2.0*m1p[i])) +
	(m1m[i]+2.0*m1p[i])+
	2.*dL2*m1p3[i] + dL4*m1m[i] 
	+ 2.0*dL9*m1p1[i];
      dF3[i] = Fb3*(f3[i] - (3.0*(e1p[i]-m1p[i]))) +
	(3.0*(e1p[i]-m1p[i])) +
	3.*(dL1*e1p3[i] - dL2*m1p3[i])
	+ 3.0*(dL8*e1p1[i] - dL9*m1p1[i]);
      dF4[i] = Fb4*f4[i];
      dF5[i] = Fb5*(f5[i] - (l0p[i]+6.0*l1p[i]*cos(t_pi_rad[i]))) +
	(l0p[i]+6.0*l1p[i]*cos(t_pi_rad[i])) +
	//Trying Nikos' suggestion to fit the real part of l0p
	//	dL7*l0p[i].Re() + 6.*cos(t_pi_rad[i])*dL3*l1p3[i]
	dL7*l0p[i] + 6.*cos(t_pi_rad[i])*dL3*l1p3[i]
	+ 6.*cos(t_pi_rad[i])*dL10*l1p1[i];
      dF6[i] = Fb6*(f6[i]-(l1m[i]-2.0*l1p[i])) +
	(l1m[i]-2.0*l1p[i]) +
	(dL5*l1m[i] - 2.*dL3*l1p3[i]) 
	- 2.0*dL10*l1p1[i];
    }
    //  printf("Re: %f %f %f %f %f %f\n",(dF1[i]-f1[i]).Re(),(dF2[i]-f2[i]).Re(),(dF3[i]-f3[i]).Re(),
    //	 (dF4[i]-f4[i]).Re(),(dF5[i]-f5[i]).Re(),(dF6[i]-f6[i]).Re());
    //  printf("Im: %f %f %f %f %f %f\n",(dF1[i]-f1[i]).Im(),(dF2[i]-f2[i]).Im(),(dF3[i]-f3[i]).Im(),
    //	 (dF4[i]-f4[i]).Im(),(dF5[i]-f5[i]).Im(),(dF6[i]-f6[i]).Im());
  //  printf("Re: %f %f %f %f %f %f\n",dF1[i].Re(),dF2[i].Re(),dF3[i].Re(),dF4[i].Re(),dF5[i].Re(),dF6[i].Re());
  //  printf("Im: %f %f %f %f %f %f\n",dF1[i].Im(),dF2[i].Im(),dF3[i].Im(),dF4[i].Im(),dF5[i].Im(),dF6[i].Im());
  }


  //dummy complex variable 
  TComplex a;
  for (int i=0; i<NUM_DATA; i++) {

    rt[i] = pow(a.Abs(dF1[i]),2) + pow(a.Abs(dF2[i]),2) 
      + pow(sin(t_pi_rad[i]),2) / 2.0
      * (pow(a.Abs(dF3[i]),2) + pow(a.Abs(dF4[i]),2))
      + (pow(sin(t_pi_rad[i]),2) 
	 * (a.Conjugate(dF2[i]) * dF3[i] + a.Conjugate(dF1[i]) * dF4[i] 
	    + cos(t_pi_rad[i]) * a.Conjugate(dF3[i]) * dF4[i])
	 - 2.0 * cos(t_pi_rad[i]) * a.Conjugate(dF1[i]) * dF2[i]).Re();
    
    rl[i]=(pow(a.Abs(dF5[i]),2)+pow(a.Abs(dF6[i]),2)+2.0*cos(t_pi_rad[i])*
	a.Conjugate(dF5[i])*dF6[i]).Re();
    
    rtt[i]=0.5*pow(sin(t_pi_rad[i]),2)*(pow(a.Abs(dF3[i]),2)+
				  pow(a.Abs(dF4[i]),2))+
      pow(sin(t_pi_rad[i]),2)*(a.Conjugate(dF2[i])*dF3[i] +
			    a.Conjugate(dF1[i])*dF4[i]+cos(t_pi_rad[i])*
			    a.Conjugate(dF3[i])*dF4[i]).Re();
    
    rlt[i]=-sin(t_pi_rad[i])*(a.Conjugate(dF2[i])*dF5[i]
			+a.Conjugate(dF3[i])*dF5[i]+a.Conjugate(dF1[i])*dF6[i] 
			+a.Conjugate(dF4[i])*dF6[i]+cos(t_pi_rad[i])*
			(a.Conjugate(dF4[i])*dF5[i]+a.Conjugate(dF3[i])*dF6[i])).Re();
    
    rltp[i]=-sin(t_pi_rad[i])*(a.Conjugate(dF2[i])*dF5[i]
			 +a.Conjugate(dF3[i])*dF5[i]+a.Conjugate(dF1[i])*dF6[i] 
			 +a.Conjugate(dF4[i])*dF6[i]+cos(t_pi_rad[i])*
			 (a.Conjugate(dF4[i])*dF5[i]+a.Conjugate(dF3[i])*dF6[i])).Im();

    rlt0y[i]=(-a.Conjugate(dF1[i])*dF5[i]+a.Conjugate(dF2[i])*dF6[i]+
                  cos(t_pi_rad[i])*(a.Conjugate(dF2[i])*dF5[i] -
				    a.Conjugate(dF1[i])*dF6[i]) +
                  pow(sin(t_pi_rad[i]),2)*(a.Conjugate(dF3[i])*dF6[i] -
					   a.Conjugate(dF4[i])*dF5[i])).Im();

    rtlp_xp0[i]=(-a.Conjugate(dF2[i])*dF5[i]+a.Conjugate(dF1[i])*dF6[i]+
		     cos(t_pi_rad[i])*(a.Conjugate(dF1[i])*dF5[i] -
				       a.Conjugate(dF2[i])*dF6[i])).Re();

    rttp_zp0[i]=(-2.*a.Conjugate(dF1[i])*dF2[i] + 
                     cos(t_pi_rad[i])*(pow(a.Abs(dF1[i]),2) +
				       pow(a.Abs(dF2[i]),2))
		     -pow(sin(t_pi_rad[i]),2)*
                     (a.Conjugate(dF1[i])*dF3[i] +a.Conjugate(dF2[i])*dF4[i])).Re();
  
    /* Convert mpoles to proper units */
    /* These response functions are defined in MAID convention */
    rlt[i]     *= MPU;
    rltp[i]     *= MPU;
    rt[i]      *= MPU; 
    rl[i]      *= MPU;
    rtt[i]     *= MPU; 
    //    rltp[i]    *= MPU;
    rlt0y[i]   *= MPU; 
    rtlp_xp0[i]*= MPU;
    rttp_zp0[i]*= MPU;

    //Set proper masses for the channels
    //pi0
    if (obs[i]<100){
      md=MASS_PROTON;
      mf=MASS_PI0;
      mpion=MASS_PI0;
      mnucl=MASS_PROTON;
    }
    //pi+
    if (obs[i]>100){
      md=MASS_PIPLUS;
      mf=MASS_NEUTRON;
      mpion=MASS_PIPLUS;
      mnucl=MASS_NEUTRON;
    }


    /* Convert to the epsilon_l fixed notation
       Multiply by |qcm*|/|omega*| */
    omega_cm[i]=(w[i]*w[i]-q2[i]*1e6-pow(MASS_PROTON,2))/(2.0*w[i]);
    E_pi_cm[i]=(pow(w[i],2)+pow(mpion,2)-pow(mnucl,2))/(2.0*w[i]);
    p_pi_cm[i]=sqrt(pow(E_pi_cm[i],2)-pow(mpion,2));
    k_gamma_cm[i]=0.5*(w[i]-pow(MASS_PROTON,2)/w[i]);
    Ppikgamma[i] = p_pi_cm[i]/k_gamma_cm[i];
    QOmegacm[i] = sqrt(q2[i]*1e6)/omega_cm[i] ;
    qcm[i] = sqrt(q2[i]*1e6 + omega_cm[i]*omega_cm[i]);
    Qqcm[i] = sqrt(q2[i]*1e6)/qcm[i];

    //MAID convention-original
    plt[i]=sqrt(2.0*eps[i]*q2[i]*1000000.0/pow(omega_cm[i],2)*(1.0+eps[i]));
    alt[i]=-plt[i]*rlt[i]/
      (rt[i]+eps[i]*q2[i]*1e6/pow(omega_cm[i],2)*rl[i]+
       eps[i]*rtt[i]);

    //Mertz convention - untested
    /*
    plt[i]=sqrt(2.0*eps[i]*q2[i]*1000000.0/pow(q_cm[i],2)*(1.0+eps[i]));
    alt[i]=-plt[i]*rlt[i]/
      (rt[i]+eps[i]*q2[i]*1e6/pow(omega_cm[i],2)*rl[i]+
       eps[i]*rtt[i]);
    */
    // Note: r0,rlt are not redefined in Mertz convention, which is not consistent
    //       with MAID convention. Sig0 and SigLT are MAID sigma.
    r0[i]  = rt[i]+rl[i]*q2[i]*1e6*eps[i]/pow(omega_cm[i],2);
    rlt[i] = rlt[i]*qcm[i]/omega_cm[i];
    rltp[i] = rltp[i]*qcm[i]/omega_cm[i];

    Sig0[i]  = Ppikgamma[i] * r0[i];
    SigLT[i] = Ppikgamma[i] * Qqcm[i] * rlt[i];
    SigLTp[i]= Ppikgamma[i] * Qqcm[i] * rltp[i];
    SigTT[i] = Ppikgamma[i] * rtt[i];

//----- Added by A Blomberg
//    printf("i: %d, s0: %g, sLT: %g, sLT': %g, sTT: %g\n",
//             i,Sig0[i],SigLT[i],SigLTp[i],SigTT[i]);

    /* q_pi/omega_cm cancelled for both numerator and denominator          
       Negative sign came from identity -Rtl_y'0=Rtl_0y */

    Py[i]=(sqrt(q2[i]*1e6)/omega_cm[i])*
      sqrt(2.0*eps[i]*(1.0+eps[i]))*
      (-rlt0y[i])/
      (rt[i] +eps[i]*q2[i]*1e6/(pow(omega_cm[i],2))*rl[i]);
	  
    /* Inserting negative sign to go from x'0 to transverse */
    Pxe[i]=(sqrt(q2[i]*1e6)/omega_cm[i])*
      sqrt(2*eps[i]*(1.0-eps[i]))*
      (-rtlp_xp0[i])/
      (rt[i] +eps[i]*q2[i]*1e6/(pow(omega_cm[i],2))*rl[i]);

    /* Inserting negative sign to go from z'0 to L */
    Pze[i]=sqrt(1.0-pow(eps[i],2))*(-rttp_zp0[i])/
      (rt[i] +eps[i]*q2[i]*1e6/(pow(omega_cm[i],2))*rl[i]);

  }


  chisq = 0 ;
  for (int j=0; j<NCTG; ++j) {
    chi[j]=0;
    ctr[j]=0;
  }

  //  double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;

  for (int i=0; i<NUM_DATA; ++i) {
    //Mod is to check for pi+ or pi0
    // Mertz Data //
    if ((obs[i] % 100) == 1) {
      chi[1] += pow((alt[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = alt[i];
      ++ctr[1];
      if (dump_chi2==1) printf("  ALT th: %f Chi2: %f\n",t_pi_rad[i]*r2d,pow((alt[i]-data_exp[i]),2)/TotErr2[i]);
    }

    if ((obs[i] % 100) == 2) {
      chi[2] += pow((rlt[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = rlt[i];
      ++ctr[2];
      if (dump_chi2==1) printf("  RLT th: %f Chi2: %f\n",t_pi_rad[i]*r2d,pow((rlt[i]-data_exp[i]),2)/TotErr2[i]);
    }

    if ((obs[i] % 100) == 3) {
      chi[3] += pow((r0[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = r0[i];
      ++ctr[3];
      if (dump_chi2==1) printf("   R0 th: %f Chi2: %f\n",t_pi_rad[i]*r2d,pow((r0[i]-data_exp[i]),2)/TotErr2[i]);
    }

    // Kunz & Nikos Data//
    if ((obs[i] % 100) == 4) {
      chi[4] += pow((Sig0[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = Sig0[i];
      ++ctr[4];
      if (dump_chi2==1) printf(" Sig0 th: %f Chi2: %f\n",t_pi_rad[i]*r2d,pow((Sig0[i]-data_exp[i]),2)/TotErr2[i]);
    }

    if ((obs[i] % 100) == 5) {
      chi[5] += pow((SigLT[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = SigLT[i];
      ++ctr[5];
      if (dump_chi2==1) printf(" SigLT th: %f Chi2: %f\n",t_pi_rad[i]*r2d, pow((SigLT[i]-data_exp[i]),2)/TotErr2[i]);
    }

    if ((obs[i] % 100) == 6) {
      chi[6] += pow((SigTT[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = SigTT[i];
      ++ctr[6];
      if (dump_chi2==1) printf(" SigTT th: %f Chi2: %f\n",t_pi_rad[i]*r2d,pow((SigTT[i]-data_exp[i]),2)/TotErr2[i]);
    }

    if ((obs[i] % 100) == 7) {
      chi[7] += pow((SigLTp[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = SigLTp[i];
      ++ctr[7];
      if (dump_chi2==1) printf("SigLTp th: %f Chi2: %f\n",t_pi_rad[i]*r2d,pow((SigLTp[i]-data_exp[i]),2)/TotErr2[i]);
    }

    // Recoil Polarizations //
    if ((obs[i] % 100) == 10) {
      chi[10] += pow((Py[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = Py[i];
      ++ctr[10];
    }

    if ((obs[i] % 100) == 11) {
      chi[11] += pow((Pxe[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = Pxe[i];
      ++ctr[11];
    }

    if ((obs[i] % 100) == 12) {
      chi[12] += pow((Pze[i]-data_exp[i]),2)/TotErr2[i];
      Theory[i] = Pze[i];
      ++ctr[12];
    }

    //Phi dependent spectrometer XS
    if ((obs[i] % 100) == 99) {
      chi[99] += pow((Sig0[i]+
		      sqrt(2*eps[i]*(1+eps[i]))*cos(ph_pi_rad[i])*SigLT[i]+
		      eps[i]*cos(2.0*ph_pi_rad[i])*SigTT[i]
		      -data_exp[i]),2)/TotErr2[i];
      if (dump_chi2==1) printf("SpecXS th: %f ph: %f Chi2: %f\n",
			       t_pi_rad[i]*r2d,
			       ph_pi_rad[i]*r2d,
			       pow((Sig0[i]+
				    sqrt(2*eps[i]*(1+eps[i]))*cos(ph_pi_rad[i])*SigLT[i]+
				    eps[i]*cos(2.0*ph_pi_rad[i])*SigTT[i]
				    -data_exp[i]),2)/TotErr2[i]
			       );
      ++ctr[99];
    }
    //Phi dependent spectrometer XS - helicity part
    if ((obs[i] % 100) == 98) {
      chi[98] += pow((-2.0*Pe*sqrt(2*eps[i]*(1-eps[i]))*
		      sin(ph_pi_rad[i])*SigLTp[i]
		      -data_exp[i]),2)/TotErr2[i];
      if (dump_chi2==1) printf("SpecXSh th: %f ph: %f Chi2: %f\n",
			       t_pi_rad[i]*r2d,
			       ph_pi_rad[i]*r2d,
			       pow((-2.0*Pe*sqrt(2*eps[i]*(1-eps[i]))*
				    sin(ph_pi_rad[i])*SigLTp[i]
				    -data_exp[i]),2)/TotErr2[i]
			       );
      ++ctr[98];
    }

    //Multipole checks
    /*
    if (dump_chi2==1){
      tmp1=e0p1[i].Re(); tmp2=e0p3[i].Re(); tmp3=e0p[i].Re();
      tmp4=e0p1[i].Im(); tmp5=e0p3[i].Im(); tmp6=e0p[i].Im();
      printf("E0+: %7.3f %7.3f %7.3f %7.3f %6.4f %7.3f %7.3f %7.3f %7.3f %6.4f\n",
	     tmp1,tmp2,tmp1+(2.0/3.0)*tmp2,tmp3,(tmp1+(2.0/3.0)*tmp2)/tmp3,
	     tmp4,tmp5,tmp4+(2.0/3.0)*tmp5,tmp6,(tmp4+(2.0/3.0)*tmp5)/tmp6);

      tmp1=l0p1[i].Re(); tmp2=l0p3[i].Re(); tmp3=l0p[i].Re();
      tmp4=l0p1[i].Im(); tmp5=l0p3[i].Im(); tmp6=l0p[i].Im();
      printf("L0+: %7.3f %7.3f %7.3f %7.3f %6.4f %7.3f %7.3f %7.3f %7.3f %6.4f\n",
	     tmp1,tmp2,tmp1+(2.0/3.0)*tmp2,tmp3,(tmp1+(2.0/3.0)*tmp2)/tmp3,
	     tmp4,tmp5,tmp4+(2.0/3.0)*tmp5,tmp6,(tmp4+(2.0/3.0)*tmp5)/tmp6);

      tmp1=e1p1[i].Re(); tmp2=e1p3[i].Re(); tmp3=e1p[i].Re();
      tmp4=e1p1[i].Im(); tmp5=e1p3[i].Im(); tmp6=e1p[i].Im();
      printf("E1+: %7.3f %7.3f %7.3f %7.3f %6.4f %7.3f %7.3f %7.3f %7.3f %6.4f\n",
	     tmp1,tmp2,tmp1+(2.0/3.0)*tmp2,tmp3,(tmp1+(2.0/3.0)*tmp2)/tmp3,
	     tmp4,tmp5,tmp4+(2.0/3.0)*tmp5,tmp6,(tmp4+(2.0/3.0)*tmp5)/tmp6);

      tmp1=m1p1[i].Re(); tmp2=m1p3[i].Re(); tmp3=m1p[i].Re();
      tmp4=m1p1[i].Im(); tmp5=m1p3[i].Im(); tmp6=m1p[i].Im();
      printf("M1+: %7.3f %7.3f %7.3f %7.3f %6.4f %7.3f %7.3f %7.3f %7.3f %6.4f\n",
	     tmp1,tmp2,tmp1+(2.0/3.0)*tmp2,tmp3,(tmp1+(2.0/3.0)*tmp2)/tmp3,
	     tmp4,tmp5,tmp4+(2.0/3.0)*tmp5,tmp6,(tmp4+(2.0/3.0)*tmp5)/tmp6);

      tmp1=l1p1[i].Re(); tmp2=l1p3[i].Re(); tmp3=l1p[i].Re();
      tmp4=l1p1[i].Im(); tmp5=l1p3[i].Im(); tmp6=l1p[i].Im();
      printf("L1+: %7.3f %7.3f %7.3f %7.3f %6.4f %7.3f %7.3f %7.3f %7.3f %6.4f\n",
	     tmp1,tmp2,tmp1+(2.0/3.0)*tmp2,tmp3,(tmp1+(2.0/3.0)*tmp2)/tmp3,
	     tmp4,tmp5,tmp4+(2.0/3.0)*tmp5,tmp6,(tmp4+(2.0/3.0)*tmp5)/tmp6);

      tmp1=m1m1[i].Re(); tmp2=m1m3[i].Re(); tmp3=m1m[i].Re();
      tmp4=m1m1[i].Im(); tmp5=m1m3[i].Im(); tmp6=m1m[i].Im();
      printf("M1-: %7.3f %7.3f %7.3f %7.3f %6.4f %7.3f %7.3f %7.3f %7.3f %6.4f\n",
	     tmp1,tmp2,tmp1+(2.0/3.0)*tmp2,tmp3,(tmp1+(2.0/3.0)*tmp2)/tmp3,
	     tmp4,tmp5,tmp4+(2.0/3.0)*tmp5,tmp6,(tmp4+(2.0/3.0)*tmp5)/tmp6);

      tmp1=l1m1[i].Re(); tmp2=l1m3[i].Re(); tmp3=l1m[i].Re();
      tmp4=l1m1[i].Im(); tmp5=l1m3[i].Im(); tmp6=l1m[i].Im();
      printf("L1-: %7.3f %7.3f %7.3f %7.3f %6.4f %7.3f %7.3f %7.3f %7.3f %6.4f\n",
	     tmp1,tmp2,tmp1+(2.0/3.0)*tmp2,tmp3,(tmp1+(2.0/3.0)*tmp2)/tmp3,
	     tmp4,tmp5,tmp4+(2.0/3.0)*tmp5,tmp6,(tmp4+(2.0/3.0)*tmp5)/tmp6);
    }
    */
    //    printf("%f %f %f %f %f %f %f\n",w[i],q2[i],t_pi_rad[i]*r2d,Sig0[i],SigTT[i],SigLT[i],SigLTp[i]);

  } // end-for-i-loop


  for (int j=0; j<NCTG; ++j) {
    chisq += chi[j];
  }

  f = chisq;

}


//
// Function name : ShowChi2Old(Double_t chisq, double chi[], int ctr[], 
//               : int NCTG, int NPAR){
//
// Description   : show chi2 (Old Version)
// Input         : chisq, chi[NCTG], ctr[NCTG], NCTG, NPAR
// Return        : 
//
/* This Routine is no longer being used.
void
ShowChi2Old(Double_t chisq, double chi[], int ctr[], int NCTG, int NPAR){

  int ctrTot = 0;
  for (int j=0; j<NCTG; ++j) ctrTot += ctr[j];

  printf("\n");
  printf("------------------------------------------------------\n");
  printf("                    Mertz      Kunz/Nikos            \n");
  printf("               ALT    RLT     R0  Sigma  Recoil  Total\n");
  printf("chi2_min    ");
  for (int i=0; i<NCTG; ++i)
    printf(" %5.1f ", chi[i]);
  printf(" %5.1f \n", chisq);
  printf("#data       ");

  for (int i=0; i<NCTG; ++i)
    printf(" %5d ", ctr[i]);
  printf(" %5d ", ctrTot);
  printf("\n");

  printf("chi2_min/DoF");
  for (int i=0; i<NCTG; ++i) {
    if ((ctr[i]-NPAR)>0) {
      printf(" %5.1f ", chi[i]/(ctr[i]-NPAR));
    }else {
      printf("   --- ");
    }
  }
  printf(" %5.1f ", chisq/(ctrTot-NPAR));
  printf("\n");

  printf("Chi2/DoF %44.1f", (chisq+DeltaChi2)/(ctrTot-NPAR));
  printf("\n");

  printf("------------------------------------------------------\n");
  printf("\n");

  return;

} */

//
// Function name : ShowChi2(Double_t chisq, double chi[], int ctr[], 
//               : int NCTG, int NPAR)
//
// Description   : show chi2
// Input         : chisq, chi[NCTG], ctr[NCTG], NCTG, NPAR
// Return        : 
//
void
ShowChi2(Double_t chisq, double chi[], int ctr[], int NCTG, int NPAR){

  int ctrTot = 0;
  for (int j=1; j<=NCTG; ++j) ctrTot += ctr[j];


  double chi2_dof[NCTG];
  for (int j=1; j<=NCTG; ++j) chi2_dof[j] = ctr[j]-NPAR>0 ? chi[j]/(ctr[j]-NPAR) : -1 ;
 
  printf("\n-------------------------------------------- \n");

  printf(" Mertz:    #data  chi2_min  chi2_min/D.o.F  \n");
  if (chi2_dof[1]>0) {
    printf("  ALT   %7d %8.1f %12.1f\n", ctr[1], chi[1], chi2_dof[1]);
  } else {
    printf("  ALT   %7d %8.1f          ---\n", ctr[1], chi[1]);
  }
  if (chi2_dof[2]>0) {
    printf("  RLT   %7d %8.1f %12.1f\n", ctr[2], chi[2], chi2_dof[2]);
  } else {
    printf("  RLT   %7d %8.1f          ---\n", ctr[2], chi[2]);
  }
  if (chi2_dof[3]>0) {
    printf("  R0    %7d %8.1f %12.1f\n", ctr[3], chi[3], chi2_dof[3]);
  } else {
    printf("  R0    %7d %8.1f          ---\n", ctr[3], chi[3]);
  }

  printf(" Kunz & Nikos:   \n");
  if (chi2_dof[4]>0) {
    printf("  Sig0  %7d %8.1f %12.1f\n", ctr[4], chi[4], chi2_dof[4]);
  } else {
    printf("  Sig0  %7d %8.1f          ---\n", ctr[4], chi[4]);
  }
  if (chi2_dof[5]>0) {
    printf("  SigLT %7d %8.1f %12.1f\n", ctr[5], chi[5], chi2_dof[5]);
  } else {
    printf("  SigLT %7d %8.1f          ---\n", ctr[5], chi[5]);
  }
  if (chi2_dof[6]>0) {
    printf("  SigTT %7d %8.1f %12.1f\n", ctr[6], chi[6], chi2_dof[6]);
  }else{
    printf("  SigTT %7d %8.1f          ---\n", ctr[6], chi[6]);
  }

  printf(" Recoil :   \n");
  if (chi2_dof[10]>0) { 
    printf("  Py    %7d %8.1f %12.1f\n", ctr[10], chi[10], chi2_dof[10]);
  } else {
    printf("  Py    %7d %8.1f          ---\n", ctr[10], chi[10]);
  }
  if (chi2_dof[11]>0) {
    printf("  Px    %7d %8.1f %12.1f\n", ctr[11], chi[11], chi2_dof[11]);
  } else {
    printf("  Px    %7d %8.1f          ---\n", ctr[11], chi[11]);
  }
  if (chi2_dof[12]>0) {
    printf("  Pz    %7d %8.1f %12.1f\n", ctr[12], chi[12], chi2_dof[12]);
  } else {
    printf("  Pz    %7d %8.1f          ---\n", ctr[12], chi[12]);
  }
  printf(" Spec XS :   \n");
  if (chi2_dof[99]>0) { 
    printf("  SpecXS%7d %8.1f %12.1f\n", ctr[99], chi[99], chi2_dof[99]);
  } else {
    printf("  SpecXS%7d %8.1f          ---\n", ctr[99], chi[99]);
  }
  printf(" Total :  %5d %8.1f %12.1f\n", ctrTot, chisq, chisq/(ctrTot-NPAR) );
  printf("-------------------------------------------- \n ");

  return;


}



 //
// Function name : ShowObsFlag()
//
// Description   : show observable flags
// Input         : 
// Return        : 
//
void
ShowObsFlag(){

  cout << "Available Observable Flags. Positive for (e,e'p) and Negative for (e,e'pi+)" << endl;
  cout << endl;
  cout << "\t ALT    :  1" << endl;
  cout << "\t RLT    :  2" << " (Note RLT and R0 are in Mertz's convention" << endl;
  cout << "\t R0     :  3" << "  and not consistent with MAID.)" << endl;
  cout << "\t Sig0   :  4" << endl;
  cout << "\t SigLT  :  5" << endl;
  cout << "\t SigTT  :  6" << endl;
  cout << "\t SigLTp :  7" << endl;
  cout << "\t Py     : 10" << endl;
  cout << "\t Px     : 11" << endl;
  cout << "\t Pz     : 12" << endl;
  cout << "\t SpecXS : 99" << endl;
  cout << endl;
  exit(0);


}


//
// Function name : fit_data(int Mode, int ShowChiSQR, int ParOut
//                          int DataOut, int CnfLevel, int ExpErr, int RtnPcode)
// Description   : control fitting mode
// Input         : int Mode, int ShowChiSQR, int ParOut
//               : int QuietMode, int DataOut, int CnfLevel,
//               : int ExpErr, int RtnPcode, int Ch2Weight, int SingleExe
// Return        : 
//
void 
fit_data(int Mode, int ShowChiSQR, int ParOut, int QuietMode, int DataOut, 
	 int CnfLevel, int ExpErr, int RtnPcode, int Ch2Weight, int SingleExe,
	 int PhiDep) {

  if (Mode%1000 == 111) {
    MAID_OUT = 1 ;
  } else if (Mode%1000<111) {
     Mode = getMODE(Mode);
     if (Mode == -111) {
       cerr << "Warning: Invalid Fitting Mode. Try -L option. " << endl;
       MAID_OUT = 1;
       Mode     = 111;
     }
     //  } else if (Mode>127) {
  } else if (Mode%1000>127) {
     cerr << "Error: Invalid Fitting Mode =" << Mode << endl;
     cerr << "       Use -L option to list available modes" << endl;
    exit(-1);
  }

  if (RtnPcode) ReturnPcodeAndExit(Mode);

  // get input data from file 
  read_data(ExpErr,PhiDep);

  // perform fit 
  int NPAR = minimize(Mode, QuietMode, CnfLevel, SingleExe, PhiDep);


  // Show Fitting Results for Individual Data 
  if (SingleExe) PrintFit();


  // get resulting parameters
  getParameters(Mode, NPAR, ParOut, DataOut, Ch2Weight);


  if (ShowChiSQR) ShowChi2(chisq, chi, ctr, NCTG, NPAR);

} // end-of-fit_data()


//
// Function name : read_data(int ExpErr)
// Description   : read input data file
// Input         : DATA_FILE
// Return        : 
//

void
read_data(int ExpErr, int PhiDep)
{
  ifstream infile(DATA_FILE);
  
  if ( infile.fail() ) {
    cout << "unable to find file:" << DATA_FILE << "\n" << flush;
    exit(-1);
  }

  int i=0;
  int ch;
  while ( ( ch = infile.peek()) != EOF ) {
    
    if (OldData==1) {
    infile >> obs[i] >> q2[i] >> w[i] >> t_pq[i] >> eps[i] >> data_exp[i] 
	   >> stat_err[i] >> inst_err[i] >> mod_err[i] >> e1pr[i] >> e1pi[i] 
	   >> m1pr[i] >> m1pi[i] >> l1pr[i] >> l1pi[i] >> m1m0r[i] >> m1m0i[i] 
	   >> l1m0r[i] >> l1m0i[i] >> m1m1r[i] >> m1m1i[i] >> l1m1r[i] 
	   >> l1m1i[i] >> f1r[i] >> f1i[i] >> f2r[i] >> f2i[i] >> f3r[i] 
	   >> f3i[i] >> f4r[i] >> f4i[i] >> f5r[i] >> f5i[i] >> f6r[i] 
	   >> f6i[i];  // 35 inputs arguments
    }else if (PhiDep){
      /*
      infile >> obs[i] >> q2[i] >> w[i] >> t_pq[i] >> ph_pq[i] >> eps[i] >> data_exp[i] 
	   >> stat_err[i] >> inst_err[i] >> mod_err[i] 
	   >> e0pr[i] >> e0pi[i] >> l0pr[i] >> l0pi[i] 
           >> e1pr[i] >> e1pi[i] >> m1pr[i] >> m1pi[i] >> l1pr[i] >> l1pi[i] 
	   >> m1mr[i] >> m1mi[i] >> l1mr[i] >> l1mi[i] 
	   >>  f1r[i] >> f1i[i] >> f2r[i] >> f2i[i] >> f3r[i] 
	   >> f3i[i] >> f4r[i] >> f4i[i] >> f5r[i] >> f5i[i] >> f6r[i] >> f6i[i]; // 36 input columns
      */
      infile >> obs[i] >> q2[i] >> w[i] >> t_pq[i] >> ph_pq[i] >> eps[i] >> data_exp[i] 
	     >> stat_err[i] >> inst_err[i] >> mod_err[i] 
	     >> e0pr[i] >> e0pi[i] >> l0pr[i] >> l0pi[i] 
	     >> e1p3r[i] >> e1p3i[i] >> m1p3r[i] >> m1p3i[i] >> l1p3r[i] >> l1p3i[i] 
	     >> m1mr[i] >> m1mi[i] >> l1mr[i] >> l1mi[i] 
	     >>  f1r[i] >> f1i[i] >> f2r[i] >> f2i[i] >> f3r[i] 
	     >> f3i[i] >> f4r[i] >> f4i[i] >> f5r[i] >> f5i[i] >> f6r[i] >> f6i[i]
	     >> e1pr[i] >> e1pi[i] >> m1pr[i] >> m1pi[i] >> l1pr[i] >> l1pi[i] ; // 36+6 input columns
      if (OldData==0){
	infile >> e0p3r[i] >> e0p3i[i] >> l0p3r[i] >> l0p3i[i] >> m1m3r[i] >> m1m3i[i] >> l1m3r[i] >> l1m3i[i] >>
	  e0p1r[i] >> e0p1i[i] >>	l0p1r[i] >> l0p1i[i] >>	e1p1r[i] >> e1p1i[i] >>	m1p1r[i] >> m1p1i[i] >>
	  l1p1r[i] >> l1p1i[i] >>	m1m1r[i] >> m1m1i[i] >>	l1m1r[i] >> l1m1i[i];
      }
    }else{
      infile >> obs[i] >> q2[i] >> w[i] >> t_pq[i] >> eps[i] >> data_exp[i] 
	     >> stat_err[i] >> inst_err[i] >> mod_err[i] 
	     >> e0pr[i] >> e0pi[i] >> l0pr[i] >> l0pi[i] 
	     >> e1p3r[i] >> e1p3i[i] >> m1p3r[i] >> m1p3i[i] >> l1p3r[i] >> l1p3i[i] 
	     >> m1mr[i] >> m1mi[i] >> l1mr[i] >> l1mi[i] 
	     >>  f1r[i] >> f1i[i] >> f2r[i] >> f2i[i] >> f3r[i] 
	     >> f3i[i] >> f4r[i] >> f4i[i] >> f5r[i] >> f5i[i] >> f6r[i] >> f6i[i]
	     >> e1pr[i] >> e1pi[i] >> m1pr[i] >> m1pi[i] >> l1pr[i] >> l1pi[i] ; // 35+6 input columns
      if (OldData==0){
	infile >> e0p3r[i] >> e0p3i[i] >> l0p3r[i] >> l0p3i[i] >> m1m3r[i] >> m1m3i[i] >> l1m3r[i] >> l1m3i[i] >>
	  e0p1r[i] >> e0p1i[i] >>	l0p1r[i] >> l0p1i[i] >>	e1p1r[i] >> e1p1i[i] >>	m1p1r[i] >> m1p1i[i] >>
	  l1p1r[i] >> l1p1i[i] >>	m1m1r[i] >> m1m1i[i] >>	l1m1r[i] >> l1m1i[i];
      }
    }

    e1p3[i]  = TComplex(e1p3r[i],  e1p3i[i]);
    m1p3[i]  = TComplex(m1p3r[i],  m1p3i[i]);
    l1p3[i]  = TComplex(l1p3r[i],  l1p3i[i]);
    m1m[i]  = TComplex(m1mr[i],  m1mi[i]);
    l1m[i]  = TComplex(l1mr[i],  l1mi[i]);
    e0p[i]  = TComplex(e0pr[i],  e0pi[i]);
    l0p[i]  = TComplex(l0pr[i],  l0pi[i]);

    e1p[i]  = TComplex(e1pr[i],  e1pi[i]);
    m1p[i]  = TComplex(m1pr[i],  m1pi[i]);
    l1p[i]  = TComplex(l1pr[i],  l1pi[i]);
      
    e0p3[i]  = TComplex(e0p3r[i],  e0p3i[i]);
    l0p3[i]  = TComplex(l0p3r[i],  l0p3i[i]);
    m1m3[i]  = TComplex(m1m3r[i],  m1m3i[i]);
    l1m3[i]  = TComplex(l1m3r[i],  l1m3i[i]);

    e0p1[i]  = TComplex(e0p1r[i],  e0p1i[i]);
    l0p1[i]  = TComplex(l0p1r[i],  l0p1i[i]);
    e1p1[i]  = TComplex(e1p1r[i],  e1p1i[i]);
    m1p1[i]  = TComplex(m1p1r[i],  m1p1i[i]);
    l1p1[i]  = TComplex(l1p1r[i],  l1p1i[i]);
    m1m1[i]  = TComplex(m1m1r[i],  m1m1i[i]);
    l1m1[i]  = TComplex(l1m1r[i],  l1m1i[i]);

    //    m1m0[i] = TComplex(m1m0r[i], m1m0i[i]);
    //    l1m0[i] = TComplex(l1m0r[i], l1m0i[i]);
    //    m1m1[i] = TComplex(m1m1r[i], m1m1i[i]);
    //    l1m1[i] = TComplex(l1m1r[i], l1m1i[i]);

    f1[i] = TComplex(f1r[i],f1i[i]);
    f2[i] = TComplex(f2r[i],f2i[i]);
    f3[i] = TComplex(f3r[i],f3i[i]);
    f4[i] = TComplex(f4r[i],f4i[i]);
    f5[i] = TComplex(f5r[i],f5i[i]);
    f6[i] = TComplex(f6r[i],f6i[i]);

    //If obs<100, then is pi0 and angle in the file is for the proton
    //If obs>100, then is pi+ and angle in the file is for the pion
    if (obs[i]<100){
      t_pi_rad[i] = (180.0 - t_pq[i])*M_PI/180.0;
      ph_pi_rad[i] = (180.0 + ph_pq[i])*M_PI/180.0;
    }
    else if (obs[i]>100){
      t_pi_rad[i] = t_pq[i]*M_PI/180.0;
      ph_pi_rad[i] = ph_pq[i]*M_PI/180.0;
    }


    // Calculate total experimental errors
    if (ExpErr==3){
      TotErr2[i] = stat_err[i]*stat_err[i] + inst_err[i]*inst_err[i] ;}
    if (ExpErr==1){
      TotErr2[i] = stat_err[i]*stat_err[i] ;}
    if (ExpErr==2){
      TotErr2[i] = inst_err[i]*inst_err[i] ;}
    if (ExpErr==7){
      TotErr2[i] = stat_err[i]*stat_err[i] + inst_err[i]*inst_err[i] 
	+ mod_err[i]*mod_err[i];}

    //Filtering of data - Needed for CLAS to eliminate zeroes.  Makes a small difference for the other fits
    //S. Stave 8-29-05: Apparently changes the NUM_DATA to one too little
    //Probably because the input loop takes the last empty line, adds one and the filter will drop it but it never used to.  Solution should be to leave the filter and not subtract 1 from it
    /*
    */
    //if (fabs(data_exp[i])>1.0e-6
    if (fabs(data_exp[i])>0.0
	//For filtering on W
	//	&&w[i]>1220-1
	//	&&w[i]<1240+1
	)
      i++ ;

  } // end-of-while loop

  getOutputKinema(i-1+1);

  for (int j=0; j<i ; ++j){
    // pick up original multipoles of a model
    if ((q2[j]==Q2)&&(w[j]==W)) {

      // resonance multipoles
      ORG_rE1p = e1pr[j];
      ORG_iE1p = e1pi[j];
      ORG_rM1p = m1pr[j];
      ORG_iM1p = m1pi[j];
      double os=(pow(w[j],2)-q2[j]*1.0e6-pow(938.27,2))/(2*w[j]);
      ORG_rS1p = l1pr[j]*sqrt(pow(os,2)+q2[j]*1.0e6)/os;
      ORG_iS1p = l1pi[j]*sqrt(pow(os,2)+q2[j]*1.0e6)/os;

      ORG_E1p = e1p[j];
      ORG_M1p = m1p[j];
      ORG_L1p = l1p[j];

      // non-resonance multipoles
      ORG_E0p = e0p[j];
      ORG_S0p = l0p[j]*sqrt(pow(os,2)+q2[j]*1.0e6)/os;
      ORG_L0p = l0p[j];
      ORG_M1m = m1m[j];
      ORG_S1m = l1m[j]*sqrt(pow(os,2)+q2[j]*1.0e6)/os;
      ORG_L1m = l1m[j];

      //F's and Fbars
      if (fabs(t_pq[j]-90.0)<1e-6){
	ORG_F1=f1[j];
	ORG_F2=f2[j];
	ORG_F3=f3[j];
	ORG_F4=f4[j];
	ORG_F5=f5[j];
	ORG_F6=f6[j];

	ORG_F1b= f1[j] - (e0p[j] + (m1p[j] + e1p[j])*3.*cos(t_pi_rad[j]));
	ORG_F2b= f2[j] - (m1m[j]+2.0*m1p[j]);
	ORG_F3b= f3[j] - (3.0*(e1p[j]-m1p[j]));
	ORG_F4b= f4[j];
	ORG_F5b= f5[j] - (l0p[j]+6.0*l1p[j]*cos(t_pi_rad[j])); 
	ORG_F6b= f6[j]-(l1m[j]-2.0*l1p[j]);
      }
	
    }

  }


  if (ORG_iE1p*ORG_iM1p*ORG_iS1p==0) {
     cerr << "Error : Data are not available at Q2=" << Q2 
     << " , W=" << W << "\n\tDesired Q2 and W have to be included in "
     << " the input data file." << endl;
    exit(-1);
  }

  infile.close();
  NUM_DATA = i;

  return ;
}


//
// Function name : getFixedParameter(int Mode, Double_t arglist[NUM_PARMS], 
//               :                    int &nFix)
//
// Description   : perform bitwise operation on Mode to interpret fix parameters 
// Input         : int Mode, Double_t arglist[NUM_PARMS], int &nFix
// Return        : 
//
void 
getFixedParameter(int Mode, Double_t arglist[NUM_PARMS], int &nFix){

  nFix=0; // number of fixed parameters
  int free;
  int tester=1;

  for (int i=NUM_PARMS-9; i>0; --i){
    free = Mode & tester ? 1 : 0 ;
    if (!free) {
      arglist[nFix] = i ;
      ++nFix;
    }
    tester *= 2;

  }
  return ;

}


//
// Function name : getDeltaChi2(int CnfLevel, int NFree)
//
// Description   : search delta_chi2 for given number of parameters 
//               : from a confidence level table
// Input         : int CnfLevel, int NFree
// Return        : delta_chi2
//
double 
getDeltaChi2(int CnfLevel, int NFree){

  if (NFree>7) {
     cerr << "Error: Number of Free Parameters " << NFree 
       << "exceeds the range of delta_chi2 database." << endl;
    exit(-1);
  } 
  if ((CnfLevel<0)||(CnfLevel>3)) {
        cerr << "Warning: Confidence Level " << CnfLevel << " is not available"
       << " in the database. \n\t delta_chi2 = 1 forced" << endl << endl;
    CnfLevel = 1;
  }

  return CnfLevel==0 ? 0 : delta_chi2[NFree][CnfLevel] ;

}

//Function to return a certain place from a number (6th digit of 9 digit number, let's say)
int strip(int x,int y,int l){

  return int(x/pow(10.0,l-y))-10*int(x/pow(10.0,l-y+1));

}


//
// Function name : minimize(int Mode, int QuietMode, int CnfLevel, int SingleExe)
//
// Description   : handle free and fix parameters and perform fitting
// Input         : int Mode, int QuietMode, int CnfLevel)
// Return        : NPAR
//

int 
minimize(int Mode, int QuietMode, int CnfLevel, int SingleExe, int PhiDep)
{

  Double_t arglist[NUM_PARMS];
  Int_t ierflg = 0;
  int NFree = 0 ; 

  PhiDep*=1;

  //  gMinuit= new TMinuit(NUM_DATA);
  gMinuit= new TMinuit(NUM_PARMS);
  gMinuit->SetFCN( fcn );

  if (QuietMode) {
    arglist[0]=-1;
    gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);
    arglist[0]=1;
    gMinuit->mnexcm("SET NOWarnings", arglist, 1, ierflg);
  }

  initArgList(arglist);


  /* set starting values and steps */
  for (int i=0; i<NUM_PARMS; i++) {
    gMinuit->mnparm( i, PARM_NAMES[i], PARM_START[i], PARM_STEP[i], 
		     0, 0, ierflg );
  }

  /* Fix and Free Parameters Manipulation */
  int nFix; // number of parameters to be fixed
  getFixedParameter(Mode%1000, arglist, nFix);
  
  if (Mode%1000!=127) gMinuit->mnexcm("FIX", arglist, nFix, ierflg);
  NFree = NUM_PARMS-9-nFix;

  /* Confidence level of parameter error */
  DeltaChi2 = getDeltaChi2(CnfLevel, NFree); 
  arglist[0] = DeltaChi2;
  if (arglist[0]) gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  /* MAID Output */
  Double_t tmp[7]={1,2,3,4,5,6,7};
  if (MAID_OUT==1){gMinuit->mnexcm("FIX", tmp, 7, ierflg);}

  //Fix the right Fbar parameters
  int i;
  for (i=0;i<6;i++){
    if (strip(Mode,i+1,9)==0) {
      tmp[0]=i+8;
      gMinuit->mnexcm("FIX", tmp, 1, ierflg);
    }
    else NFree++;
  }

  //Allow individual p1/2 mpoles to vary
  if (p12>=0){
    gMinuit->mnparm( p12+13, PARM_NAMES[p12+13], 
		     PARM_START[p12+13], 0.1, 
		     0, 0, ierflg );
    NFree++;
  }

    

  /* Minimize It! */
  initArgList(arglist);
  arglist[0] = 100000;            // do at least 1000 function calls
  arglist[1] = 0.1;             // tolerance = 0.1

  if (!SingleExe) {
    //  gMinuit->mnexcm("HESSE", arglist, 2, ierflg );
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg );
  } else {
    //Just call FCN once for debugging
    gMinuit->mnexcm("call", arglist, 2, ierflg );
  }

  if (NFree==4){
    gMinuit->mnemat(&err_mat[0][0],4);
  }
  if (NFree==3){
    gMinuit->mnemat(&err_mat3[0][0],3);
  }

  //Do not want the chi2 output all the time.  Only once at the end.
  dump_chi2=DumpChi2;
  gMinuit->mnexcm("call", arglist, 2, ierflg );
  dump_chi2=0;

  return NFree ;

}


//
// Function name : minimize_org(int NPAR)
//
// Description   : handle free and fix parameters and perform fitting
//               : the minimze function is upgraded to minimze() though,
//               : original function is kept here to reproduce previous results
// Input         : int NPAR
// Return        : 
//

void minimize_org(int NPAR) 
{
  cout << "performing minimization of all parameters...\n" << flush;
  gMinuit = new TMinuit(NUM_DATA);
  gMinuit->SetFCN( fcn );

  Double_t arglist[10];
  Int_t ierflg = 0;

  /* Set starting values and steps */
  cout << "setting starting values and steps...\n" << flush;
  for (int i=0; i<NUM_PARMS; i++) {
    gMinuit->mnparm( i, PARM_NAMES[i], PARM_START[i], PARM_STEP[i], 
		     0, 0, ierflg );
  }
  
  if (NPAR==5) {
    arglist[0] = 6;
    arglist[1] = 7;
    gMinuit->mnexcm("FIX", arglist, 2, ierflg);
  }else if (NPAR==3){
    arglist[0] = 4;
    arglist[1] = 5;
    arglist[2] = 6;
    arglist[3] = 7;
    gMinuit->mnexcm("FIX", arglist, 4, ierflg);
  }

  /* minimize it */
  cout << "running migrad...\n" << flush;
  arglist[0] = 10000;            // do at least 1000 function calls
  arglist[1] = 0.01;             // tolerance = 0.1
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg );

}


//
// Function name : PrintFit
//
// Description   : Print fitting resuls for all individual data points
// Input         : 
// Return        : 
//
void 
PrintFit(){

  printf("----------------------------------------------------------\n");
  printf("   obs   Q2      W     th*piq    eps     exp    Theory\n");
  printf("----------------------------------------------------------\n");


  for (int i=0; i<NUM_DATA; ++i) {
    printf("%5d",         obs[i]);
    printf("%8.3f",        q2[i]);
    printf("%8.1f",         w[i]);
    printf("%8.2f",  180-t_pq[i]);
    printf("%8.3f",       eps[i]);
    printf("%8.3f",  data_exp[i]); 
    printf("%9.4f",    Theory[i]);
    printf("\n");
  }

  return;

}

//Routine to find the eigenvectors and put them in order
void print_eig(int NPAR,Double_t var[]){
  Double_t eigen[4][4];
  Double_t eigen_err[4];
  Int_t mapit[4]={0,1,2,3};
  Double_t tmp_vec[4];
  Double_t eig_err[4];
  Int_t i,j;
  
  if (NPAR==3){
    Double_t err2_vec[9]={err_mat3[0][0],err_mat3[0][1],err_mat3[0][2],
			  err_mat3[1][0],err_mat3[1][1],err_mat3[1][2],
			  err_mat3[2][0],err_mat3[2][1],err_mat3[2][2]};
    
    TMatrixD err2(3,3,err2_vec);
    
    TVectorD eig_vals(3);
    
    TMatrixD curv = err2.Invert();
    err2.Invert();
    TMatrixD eig = curv.EigenVectors(eig_vals);
    
    Double_t eig_val_vec[9]={eig_vals[0],0,0,
			     0,eig_vals[1],0,
			     0,0,eig_vals[2]};
    
    TMatrixD eig_val_mat(3,3,eig_val_vec);
    
    //    printf("Eigenvectors:\n");
    //    eig.Print();
    //    printf("Eigenvalues:\n");
    //    eig_vals.Print();
    
    printf("Errors on eigenvectors:\n");
    printf("d1': %f\n",1.0/sqrt(eig_vals[0]));
    printf("d2': %f\n",1.0/sqrt(eig_vals[1]));
    printf("d3': %f\n",1.0/sqrt(eig_vals[2]));
  
    for (i=0;i<3;i++) eig_err[i]=1.0/sqrt(eig_vals[i]);
    
    //Map the eigenvectors back to the parameter order
    //Simply taking the largest value as the primary direction in the parameter space
    for (i=0;i<3;i++){
      tmp_vec[0]=fabs(eig[i][0]);
      tmp_vec[1]=fabs(eig[i][1]);
      tmp_vec[2]=fabs(eig[i][2]);
      mapit[i]=TMath::LocMax(3,tmp_vec);
    }
    
    eigen_err[0]=eig_err[mapit[0]];
    eigen_err[1]=eig_err[mapit[1]];
    eigen_err[2]=eig_err[mapit[2]];
    
    for (i=0;i<3;i++){
      eigen[0][i]=eigen_err[0]*eig[i][mapit[0]];
      eigen[1][i]=eigen_err[1]*eig[i][mapit[1]];
      eigen[2][i]=eigen_err[2]*eig[i][mapit[2]];
    }
    
    printf("Re-ordered eigenvectors:\n");
    for (i=0;i<3;i++)
      printf("%f, %f, %f\n",eigen[0][i]/eigen_err[0],
	     eigen[1][i]/eigen_err[1],
	     eigen[2][i]/eigen_err[2]);
    
    printf("Re-ordered eigenvector errors:\n");
    printf("%f, %f, %f\n",eigen_err[0],eigen_err[1],eigen_err[2]);
    
    printf("plot_obs_err_fits.C format:\n");
    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	printf("v[%i][%i]=%f;\n",i,j,eigen[j][i]/eigen_err[j]);
      }
    }
    for (i=0;i<3;i++) printf("v_err[%i]=%f;\n",i,eigen_err[i]);
    for (i=0;i<3;i++) printf("var[%i]=%f;\n",i,var[i]);
  }

  if (NPAR==4){
    Double_t err2_vec[16]={err_mat[0][0],err_mat[0][1],err_mat[0][2],err_mat[0][3],
			   err_mat[1][0],err_mat[1][1],err_mat[1][2],err_mat[1][3],
			   err_mat[2][0],err_mat[2][1],err_mat[2][2],err_mat[2][3],
			   err_mat[3][0],err_mat[3][1],err_mat[3][2],err_mat[3][3]};
    
    TMatrixD err2(4,4,err2_vec);
    
    TVectorD eig_vals(4);
    
    TMatrixD curv = err2.Invert();
    err2.Invert();
    TMatrixD eig = curv.EigenVectors(eig_vals);
    
    Double_t eig_val_vec[16]={eig_vals[0],0,0,0,
			      0,eig_vals[1],0,0,
			      0,0,eig_vals[2],0,
			      0,0,0,eig_vals[3]};
    
    TMatrixD eig_val_mat(4,4,eig_val_vec);
    
    printf("Eigenvectors:\n");
    eig.Print();
    printf("Eigenvalues:\n");
    eig_vals.Print();
    
    printf("Errors on eigenvectors:\n");
    printf("d1': %f\n",1.0/sqrt(eig_vals[0]));
    printf("d2': %f\n",1.0/sqrt(eig_vals[1]));
    printf("d3': %f\n",1.0/sqrt(eig_vals[2]));
    printf("d4': %f\n",1.0/sqrt(eig_vals[3]));
  
    for (i=0;i<4;i++) eig_err[i]=1.0/sqrt(eig_vals[i]);
    
    
    //Map the eigenvectors back to the parameter order
    //Simply taking the largest value as the primary direction in the parameter space
    for (i=0;i<4;i++){
      tmp_vec[0]=fabs(eig[i][0]);
      tmp_vec[1]=fabs(eig[i][1]);
      tmp_vec[2]=fabs(eig[i][2]);
      tmp_vec[3]=fabs(eig[i][3]);
      mapit[i]=TMath::LocMax(4,tmp_vec);
    }
    
    eigen_err[0]=eig_err[mapit[0]];
    eigen_err[1]=eig_err[mapit[1]];
    eigen_err[2]=eig_err[mapit[2]];
    eigen_err[3]=eig_err[mapit[3]];
    
    for (i=0;i<4;i++){
      eigen[0][i]=eigen_err[0]*eig[i][mapit[0]];
      eigen[1][i]=eigen_err[1]*eig[i][mapit[1]];
      eigen[2][i]=eigen_err[2]*eig[i][mapit[2]];
      eigen[3][i]=eigen_err[3]*eig[i][mapit[3]];
    }
    
    printf("Re-ordered eigenvectors:\n");
    for (i=0;i<4;i++)
      printf("%f, %f, %f, %f\n",eigen[0][i]/eigen_err[0],
	     eigen[1][i]/eigen_err[1],
	     eigen[2][i]/eigen_err[2],
	     eigen[3][i]/eigen_err[3]);
    
    printf("Re-ordered eigenvector errors:\n");
    printf("%f, %f, %f, %f\n",eigen_err[0],eigen_err[1],eigen_err[2],eigen_err[3]);

    printf("plot_obs_err_fits.C format:\n");
    for (i=0;i<4;i++){
      for (j=0;j<4;j++){
	printf("v[%i][%i]=%f;\n",i,j,eigen[j][i]/eigen_err[j]);
      }
    }
    for (i=0;i<4;i++) printf("v_err[%i]=%f;\n",i,eigen_err[i]);
    for (i=0;i<4;i++) printf("var[%i]=%f;\n",i,var[i]);
      
  }



}


//
// Function name : getParameters()
//
// Description   : get parameters and errors
// Input         : int Mode, int NPAR, int ParOut, int DataOut, Ch2Weight
// Return        : 
//
void
getParameters(int Mode, int NPAR, int ParOut, int DataOut, int Ch2Weight){
  Double_t var[NUM_PARMS], verr[NUM_PARMS];
  TString Tag[NUM_PARMS];
  Int_t ivarbl;
  Double_t bnd1, bnd2;
  Double_t Ch2DoF = NUM_DATA-NPAR ? chisq/(NUM_DATA-NPAR) : -1 ;

  for (int i=0; i<NUM_PARMS; i++) {
    gMinuit->mnpout( i, Tag[i], var[i], verr[i], bnd1, bnd2, ivarbl);
    if (Ch2Weight) verr[i] *= sqrt(Ch2DoF) ;
  }

  // fitting result in relative unit to original amplitudes 
  double fitE1p = var[0];
  double fitM1p = var[1];
  double fitS1p = var[2];
  double fitM1m = var[3];
  double fitS1m = var[4];
  double fitE0p = var[5];
  double fitS0p = var[6];

  double E1p  = ORG_iE1p * var[0];
  double dE1p = ORG_iE1p * verr[0];
  double M1p  = ORG_iM1p * var[1];
  double dM1p = ORG_iM1p * verr[1];
  double S1p  = ORG_iS1p * var[2];
  double dS1p = ORG_iS1p * verr[2];
  double EMR  = E1p/M1p * 100;
  double dEMR = 100*Error(E1p, dE1p, M1p, dM1p);
  double CMR  = S1p/M1p * 100;
  double dCMR = 100*Error(S1p, dS1p, M1p, dM1p);

  //dummy complex
  TComplex a;

  if (DataOut==1) {

    printf("%4d",Mode);
    printf("%7.3f%7.1f ", Q2, W);
    printf("%7.2f%7.2f ", E1p, dE1p);
    printf("%7.2f%7.2f ", M1p, dM1p);
    printf("%7.2f%7.2f ", S1p, dS1p);
    printf("%7.2f%7.2f", EMR, Abs(dEMR));
    printf("%7.2f%7.2f", CMR, Abs(dCMR));
    printf("%7.2f%4d", chisq, NUM_DATA-NPAR);
    NUM_DATA-NPAR ? printf("%7.2f", chisq/(NUM_DATA-NPAR)) : printf("  -1");
    printf("\n");

  } else if (DataOut==2) {

    printf("%4d%4d",Mode,getPCODE(Mode));
    printf("%7.3f%7.1f ", Q2, W);
    printf("%7.3f%7.3f ", ORG_rE1p, ORG_iE1p);
    printf("%7.3f%7.3f ", ORG_rM1p, ORG_iM1p);
    printf("%7.3f%7.3f ", ORG_rS1p, ORG_iS1p);
    printf("%7.3f%7.3f%7.3f ", fitE1p, fitM1p, fitS1p);
    printf("%7.3f%7.3f%7.3f%7.3f ", fitE0p, fitS0p, fitM1m, fitS1m);
    printf("\n");


  } else if (DataOut==3) {


    printf("%4d",getPCODE(Mode));
    printf("%7.3f%7.1f ", Q2, W);
    printf("%6.2f%6.2f ", E1p, dE1p);
    printf("%6.2f%6.2f ", M1p, dM1p);
    printf("%6.2f%6.2f ", S1p, dS1p);
    printf("%6.2f%6.2f", EMR, Abs(dEMR));
    printf("%6.2f%6.2f", CMR, Abs(dCMR));
    printf("%6.2f%6.2f", ORG_E0p.Re(), ORG_E0p.Im());
    printf("%7.2f%6.2f", var[5]==1 ? -99 : var[5], verr[5]);
    printf("%6.2f%6.2f", ORG_S0p.Re(), ORG_S0p.Im());
    printf("%7.2f%6.2f", var[6]==1 ? -99 : var[6], verr[6]);
    printf("%6.2f%6.2f", ORG_M1m.Re(), ORG_M1m.Im());
    printf("%7.2f%6.2f", var[3]==1 ? -99 : var[3], verr[3]);
    printf("%7.2f%6.2f", ORG_S1m.Re(), ORG_S1m.Im());
    printf("%7.2f%6.2f", var[4]==1 ? -99 : var[4], verr[4]);
    printf("\n");

  } else if (DataOut==4){
    //Basic fit output for the 3 res pars and the 3 p1/2 pars
    printf("%f %f %f %f %f %f ",var[0],verr[0],var[1],verr[1],var[2],verr[2]);
    printf("%f %f %f %f %f %f ",var[13],verr[13],var[14],verr[14],var[15],verr[15]);
    printf("%f ",chisq);
    //Adding L0+
    printf("%f %f \n",var[6],verr[6]);
  } else if (DataOut==5){
    //Basic fit output for the 3 res pars and the 3 p1/2 pars
    printf("%f %f %f %f %f %f ",var[0],verr[0],var[1],verr[1],var[2],verr[2]);
    printf("%f %f %f %f %f %f ",var[3],verr[3],var[4],verr[4],var[5],verr[5]);
    printf("%f %f %f %f %f %f ",var[6],verr[6],var[7],verr[7],var[8],verr[8]);
    printf("%f %f %f %f %f %f ",var[9],verr[9],var[10],verr[10],var[11],verr[11]);
    printf("%f %f %f %f %f %f ",var[12],verr[12],var[13],verr[13],var[14],verr[14]);
    printf("%f %f ",var[15],verr[15]);
    printf("%f \n",chisq);

  } else {

    printf("\n");
    printf("-------------------------------------------\n");
    printf("  Q^2 = %5.3f (GeV/c)^2, W = %5.1f (MeV)\n", Q2, W);
    printf("-------------------------------------------\n");
    if (ParOut) {
      printf("           par        error       error[%%]\n");
      for (int i=0; i<NUM_PARMS; i++) {
	cout << "  " << Tag[i] ;
	if (fabs(var[i])>1e-8)
	  printf(" %10.6f  %10.6f  %10.2f\n", var[i], verr[i], verr[i]/var[i]*100);
	else
	  printf(" %10.6f  %10.6f  %10s\n", var[i], verr[i], "---");

      }
    } else {
      printf("           Mult        error       error[%%]\n");
      printf("  E1+ %10.2f  %10.2f  %11.1f\n", E1p, dE1p, dE1p/E1p*100);
      printf("  M1+ %10.2f  %10.2f  %11.1f\n", M1p, dM1p, dM1p/M1p*100);
      printf("  S1+ %10.2f  %10.2f  %11.1f\n", S1p, dS1p, dS1p/S1p*100);
    }
    printf("  EMR  = %7.2f +/- %4.2f\n", EMR, Abs(dEMR));
    printf("  CMR  = %7.2f +/- %4.2f\n", CMR, Abs(dCMR));
    printf("-------------------------------------------\n");
    printf("NUM_DATA: %i NPAR: %i\n",NUM_DATA,NPAR);
    printf("chi2 = %5.2f, d.o.f=%3d,",chisq,NUM_DATA-NPAR);
    NUM_DATA-NPAR ? printf(" chi2/d.o.f = %5.2f\n", chisq/(NUM_DATA-NPAR)) : printf("chi2/d.o.f = n/a\n");
    printf("-------------------------------------------\n");




// note above calculations for CMR are done with the 
// wrong multipoles  (i.e. charge channel)
  double os=(pow(w[0],2)-q2[0]*1.0e6-pow(938.27,2))/(2*w[0]);
  double E1a  = e1p3i[0] * var[0];
  double dE1a = e1p3i[0] * verr[0];
  double M1a  = m1p3i[0] * var[1];
  double dM1a = m1p3i[0] * verr[1];
  double S1a  = l1p3i[0] * sqrt(pow(os,2)+q2[0]*1.0e6)/os * var[2];
  double dS1a = l1p3i[0] * sqrt(pow(os,2)+q2[0]*1.0e6)/os * verr[2];
  double EMRa  = E1a/M1a * 100;
  double dEMRa = 100*Abs(Error(E1a, dE1a, M1a, dM1a));
  double CMRa  = S1a/M1a * 100;
  double dCMRa = 100*Abs(Error(S1a, dS1a, M1a, dM1a));

// note above assumes imaginary parts go to zero!
// calc below without that assumption... same result checked
/*  double fEMR = 100*( (var[0]*var[1]*(e1p3r[0]*m1p3r[0]+e1p3i[0]*m1p3i[0]))
                  / (var[1]*var[1]*(m1p3r[0]*m1p3r[0]+m1p3i[0]*m1p3i[0]) ) );
  double fCMR = 100*( (sqrt(pow(os,2)+q2[0]*1.0e6)/os*
                        var[2]*var[1]*(l1p3r[0]*m1p3r[0]+l1p3i[0]*m1p3i[0]))
                  / (var[1]*var[1]*(m1p3r[0]*m1p3r[0]+m1p3i[0]*m1p3i[0]) ) );
*/
    printf("EMR adam %7.2f +/- %4.2f\n",EMRa,dEMRa);
//    printf("EMR full %7.2f +/- %4.2f\n",fEMR,0.0);

    printf("CMR adam %7.2f +/- %4.2f\n",CMRa,dCMRa);
//    printf("CMR full %7.2f +/- %4.2f\n",fCMR,0.0);

    printf("M1  adam %7.2f +/- %4.2f\n",M1a,dM1a);
    printf("\n");
  }


  //All extras: eigenvectors, error mats etc
  //xoxo
  if (0){

    for (int j=0;j<16;j++){
      printf("par[%i]=%f;\n",j,var[j]);
    }
    
    Int_t extra;
    if (Mode==113) extra=6;
    else if (Mode==114) extra=5;
    else if (Mode==116) extra=4;
    else if (Mode==120) extra=3;
    else if (Mode==100000112) extra=7;
    else if (Mode==10000112) extra=8;
    else if (Mode==1000112) extra=9;
    else if (Mode==100112) extra=10;
    else if (Mode==10112) extra=11;
    else if (Mode==1112) extra=12;
    else if (p12==0) extra=13;
    else if (p12==1) extra=14;
    else if (p12==2) extra=15;
    else extra=3;
    
    printf("plot_4par_fits.C format:\n");
    printf("%f %f %f %f %f %f %f %f %f %i\n",
	   var[0],verr[0],
	   var[1],verr[1],
	   var[2],verr[2],
	   var[extra],verr[extra],
	   chisq,NUM_DATA-NPAR);
    
    Double_t r12,r23,r13;
    
    if (NPAR==3){
      printf("Error matrix (3x3):\n");
      printf("%e %e %e \n",
	     err_mat3[0][0],err_mat3[0][1],err_mat3[0][2]);
      printf("%e %e %e \n",
	     err_mat3[1][0],err_mat3[1][1],err_mat3[1][2]);
      printf("%e %e %e \n",
	     err_mat3[2][0],err_mat3[2][1],err_mat3[2][2]);
      printf("Corr coeff:\n");
      r12=err_mat3[0][1]/sqrt(err_mat3[0][0]*err_mat3[1][1]);
      r13=err_mat3[0][2]/sqrt(err_mat3[0][0]*err_mat3[2][2]);
      r23=err_mat3[1][2]/sqrt(err_mat3[1][1]*err_mat3[2][2]);
      printf("r12: %10.7f\n",r12);
      printf("r13: %10.7f\n",r13);
      printf("r23: %10.7f\n",r23);
      printf("Range: %10.7f < r12 < %10.7f\n",
	     rho(Z(r12)-1.0/sqrt(1.0*(NUM_DATA-NPAR-3))),
	     rho(Z(r12)+1.0/sqrt(1.0*(NUM_DATA-NPAR-3))));
      printf("Range: %10.7f < r13 < %10.7f\n",
	     rho(Z(r13)-1.0/sqrt(1.0*(NUM_DATA-NPAR-3))),
	     rho(Z(r13)+1.0/sqrt(1.0*(NUM_DATA-NPAR-3))));
      printf("Range: %10.7f < r23 < %10.7f\n",
	     rho(Z(r23)-1.0/sqrt(1.0*(NUM_DATA-NPAR-3))),
	     rho(Z(r23)+1.0/sqrt(1.0*(NUM_DATA-NPAR-3))));
      
      /*
	printf("%10s %10.7f %10.7f\n","---",
	err_mat3[0][1]/sqrt(err_mat3[0][0]*err_mat3[1][1]),
	err_mat3[0][2]/sqrt(err_mat3[0][0]*err_mat3[2][2]));
	printf("%10.7f %10s %10.7f\n",
	err_mat3[0][1]/sqrt(err_mat3[0][0]*err_mat3[1][1]),
	"---",
	err_mat3[1][2]/sqrt(err_mat3[1][1]*err_mat3[2][2]));
	printf("%10.7f %10.7f %10s\n",
	err_mat3[0][2]/sqrt(err_mat3[0][0]*err_mat3[2][2]),
	err_mat3[1][2]/sqrt(err_mat3[1][1]*err_mat3[2][2]),
	"---");
      */
    }
    if (NPAR==4){
      printf("Error matrix (4x4):\n");
      printf("%e %e %e %e\n",
	     err_mat[0][0],err_mat[0][1],err_mat[0][2],err_mat[0][3]);
      printf("%e %e %e %e\n",
	     err_mat[1][0],err_mat[1][1],err_mat[1][2],err_mat[1][3]);
      printf("%e %e %e %e\n",
	     err_mat[2][0],err_mat[2][1],err_mat[2][2],err_mat[2][3]);
      printf("%e %e %e %e\n",
	     err_mat[3][0],err_mat[3][1],err_mat[3][2],err_mat[3][3]);
    }
    
    Double_t var_tmp[4];
    var_tmp[0]=var[0];
    var_tmp[1]=var[1];
    var_tmp[2]=var[2];
    var_tmp[3]=var[extra];
    
    print_eig(NPAR,var_tmp);
    
    printf("Abs orig mpoles: E0+ L0+ E1+ M1+ L1+ M1- L1-\n");
    printf("%f %f %f %f %f %f %f\n",
	   a.Abs(ORG_E0p),
	   a.Abs(ORG_L0p),
	   a.Abs(ORG_E1p),
	   a.Abs(ORG_M1p),
	   a.Abs(ORG_L1p),
	   a.Abs(ORG_M1m),
	   a.Abs(ORG_L1m));
    
    printf("Abs orig Fbars: F1b F2b f3b f4b f5b f6b\n");
    printf("%f %f %f %f %f %f \n",
	   a.Abs(ORG_F1b),
	   a.Abs(ORG_F2b),
	   a.Abs(ORG_F3b),
	   a.Abs(ORG_F4b),
	   a.Abs(ORG_F5b),
	   a.Abs(ORG_F6b));
  }//end print extras

  return;
}



//
// Function name : getPCODE(int Mode)
//
// Description   : return Parameter Code
// Input         : int Mode
// Return        : Parameter Code
//
int 
getPCODE(int Mode){

  int i=0;
  while (Mode != PCODE[0][i]) ++i;

  return PCODE[1][i];

}

//
// Function name : getMODE(int Pcode)
//
// Description   : return Mode matches with Parameter Code
// Input         : int Pcode
// Return        : Mode. No fit mode if Pcode doesn't match 
//
int 
getMODE(int Pcode){

  int i=0;
  while ( (Pcode != PCODE[1][i])&&(i<=NMODE) ) ++i;

  return i!=NMODE+1 ? PCODE[0][i] : -111 ;

}

//
// Function name : ReturnPcodeAndExit(int Mode)
//
// Description   : return Pcode and Exit
// Input         : int Mode
// Return        : 
//
int 
ReturnPcodeAndExit(int Mode){

  cout << getPCODE(Mode) << endl;
  exit(0);

  return 1;
}

void 
initArgList(Double_t arglist[NUM_PARMS]){

  for (int i=0; i<NUM_PARMS; ++i) arglist[i] = 0 ;

  return;
}


double 
Error(double  x, double  dx, double y, double dy){

  return x/y*sqrt(dx*dx/x/x + dy*dy/y/y);

}


//
// Function name : Abs(double x)
//
// Description   : return absolute of x
// Input         : x
// Return        : abs(x)
//
double 
Abs(double X) { 
  return X>=0 ? X : -(X) ; 
}









