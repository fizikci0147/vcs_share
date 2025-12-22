#ifndef FIT_DATA_H
#define FIT_DATA_H

/* constants */
const int MAX_DATA = 4000;
const int NUM_PARMS = 7+6+3;
const TString PARM_NAMES[NUM_PARMS] = {"E1+", "M1+", "L1+", 
				       "M1-", "L1-", "E0+", "L0+",
				       "Fb1","Fb2","Fb3","Fb4","Fb5","Fb6",
				       "E1+(p1/2)","M1+(p1/2)","L1+(p1/2)"};
Double_t PARM_START[NUM_PARMS] = { 1 ,1,1,
				  1, 1, 1, 1, 
				  1,1,1,1,1,1,
				  1,1,1};
//9-15-05: W=1221 M2003 fit from data in Q06_W1221_fit_data_Sep1105.sxc
//const Double_t PARM_START[NUM_PARMS] = {1.317, 0.974, 0.753, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1};
//const Double_t PARM_START[NUM_PARMS] = {1.3907, 0.9756, 1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1,1};
Double_t PARM_STEP[NUM_PARMS] = {0.1, 0.1, 0.1, 
				 0.1, 0.1, 0.1, 0.1, 
				 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
				 0.0,0.0,0.0};
//const Double_t PARM_STEP[NUM_PARMS] = {0.0, 0.0, 0.0, 
//				       0.0, 0.0, 0.0, 0.0, 
//				       0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
//				       0.0,0.0,0.0};
int p12;
const int NMODE = 17 ;
const int PCODE[2][NMODE] ={{111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127},
			    {-1,  31, 41, 42, 51, 43, 52, 53, 61, 44, 54, 55, 62, 56, 63, 64, 71}};

/* Particle Mass */
const double MASS_PROTON  = 938.2723128;
const double MASS_NEUTRON = 939.5656;
const double MASS_PIPLUS  = 139.56995;
const double MASS_PI0     = 134.9764;

/* Plank Constant [MeV*cm] */
const double hc           = 197.327053;
/* Mass Pion Unit */
const double MPU          = 0.01*(hc/MASS_PIPLUS)*(hc/MASS_PIPLUS);

const double r2d=180.0/3.14159265359;

//const double M_PI=3.14159265359;
Double_t Pe=-0.73;
/* for comparison of original model prediction */
//double DEF_Q2 = 0.127;
//Bates
//double DEF_Q2 = 0.126;
//double DEF_Q2 = 0.127;
//double DEF_W  = 1232;
//Mainz
//double DEF_Q2 = 0.06;
//double DEF_W  = 1221;

double DEF_Q2 = 0.09;
double DEF_W  = 1232;


//double DEF_W  = 1225;
//double DEF_W  = 1155;
//CLAS
//double DEF_Q2 = 0.16;
//double DEF_W  = 1220;
//double DEF_Q2 = 0.0;
//double DEF_W  = 1084.83;
double W = 0;
double Q2 = 0;
double ORG_iM1p = 0;
double ORG_iE1p = 0;
double ORG_iS1p = 0;
double ORG_rM1p = 0;
double ORG_rE1p = 0;
double ORG_rS1p = 0;

double ORG_iM1m = 0;
double ORG_iL1m = 0;
double ORG_iL0p = 0;
double ORG_iE0p = 0;

TComplex ORG_E0p(0,0);
TComplex ORG_S0p(0,0);
TComplex ORG_L0p(0,0);
TComplex ORG_M1m(0,0);
TComplex ORG_S1m(0,0);
TComplex ORG_L1m(0,0);

TComplex ORG_E1p(0,0);
TComplex ORG_M1p(0,0);
TComplex ORG_L1p(0,0);

TComplex ORG_F1b(0,0);
TComplex ORG_F2b(0,0);
TComplex ORG_F3b(0,0);
TComplex ORG_F4b(0,0);
TComplex ORG_F5b(0,0);
TComplex ORG_F6b(0,0);

TComplex ORG_F1(0,0);
TComplex ORG_F2(0,0);
TComplex ORG_F3(0,0);
TComplex ORG_F4(0,0);
TComplex ORG_F5(0,0);
TComplex ORG_F6(0,0);


TMinuit *gMinuit;


/* default data file name */ 
TString DATA_FILE = "FitData_Fbar.inp";

/* DataFormat Indicator */
int OldData = 0;

/* MAID Output flag */
int MAID_OUT = 0;

/* Fitting Results Parameters */
int NUM_DATA; // total number of experimental data points
int NCTG=100;   // number of categories to classify observables 
int ctr[100];   // number of experimental data counter for each category
double chi[100];
Double_t chisq = 0;
Double_t DeltaChi2 = 0;  // DeltaChi2 

/* Multi-parameter confidence region, DeltaChi2 
   Reference: MINUIT Reference Manual v94.1 */
double delta_chi2[8][4] = {{0.00, 0.00,  0.00,  0.00},
			   {0.00, 1.07,  3.84,  6.63},
			   {0.00, 2.41,  5.99,  9.21},
			   {0.00, 3.67,  7.82, 11.36},
			   {0.00, 4.88,  9.49, 13.28},
			   {0.00, 6.06, 11.07, 15.09},
			   {0.00, 7.23, 12.59, 16.81},
			   {0.00, 8.38, 14.07, 18.49}};

Int_t dump_chi2=0;
Int_t DumpChi2 =0;//This is the one you set for chi2 dump
/* function prototypes */
void fcn(Int_t &, Double_t *, Double_t &, Double_t *, Int_t );
void read_data(int ExpErr, int PhiDep);
void minimize_org(int NPAR);
int minimize(int Mode, int QuietMode, int CnfLevel, int SingleExe, int PhiDep);
void initArgList(Double_t arglist[NUM_PARMS]);
void getFixedParameter(int Mode, Double_t arglist[NUM_PARMS], int &nFix);
double getDeltaChi2(int CnfLevel, int NFree);
void fit_data(int Mode, int ShowChiSQR, int ParOut,int QuietMode, int DataOut, 
	      int CnfLevel, int ExpErr, int RtnPcode, int Ch2Weight, int SingleExe, int PhiDep);
int getOutputKinema(int NDATA);
int ReturnPcodeAndExit(int Mode);

void ShowObsFlag();
void ListMode();
void ShowChi2(Double_t chisq, double chi[], int ctr[], int NCTG, int NPAR);
void PrintFit();
void getParameters(int Mode, int NPAR, int ParOut, int DataOut, int Ch2Weight);
int getPCODE(int Mode);
int getMODE(int Pcode);
double Error(double x, double dx, double y, double dy);
double Abs(double X);

/* global variables */
TComplex e0p[MAX_DATA], l0p[MAX_DATA], 
  e1p[MAX_DATA], m1p[MAX_DATA], l1p[MAX_DATA], 
  m1m[MAX_DATA], l1m[MAX_DATA], 
  m1m0[MAX_DATA], l1m0[MAX_DATA], 
  //m1m1[MAX_DATA], l1m1[MAX_DATA], 
  f1[MAX_DATA], f2[MAX_DATA], f3[MAX_DATA], f4[MAX_DATA], f5[MAX_DATA], f6[MAX_DATA],
  e0p3[MAX_DATA], l0p3[MAX_DATA],
  e1p3[MAX_DATA], m1p3[MAX_DATA], l1p3[MAX_DATA], 
  m1m3[MAX_DATA], l1m3[MAX_DATA],
  e0p1[MAX_DATA], l0p1[MAX_DATA],
  e1p1[MAX_DATA], m1p1[MAX_DATA], l1p1[MAX_DATA], 
  m1m1[MAX_DATA], l1m1[MAX_DATA];
  

int obs[MAX_DATA];
Double_t q2[MAX_DATA], w[MAX_DATA], t_pq[MAX_DATA], ph_pq[MAX_DATA], eps[MAX_DATA],
  data_exp[MAX_DATA], stat_err[MAX_DATA], inst_err[MAX_DATA], mod_err[MAX_DATA], 
  TotErr2[MAX_DATA],
  e1pr[MAX_DATA], e1pi[MAX_DATA], 
  m1pr[MAX_DATA],  m1pi[MAX_DATA],  l1pr[MAX_DATA],  l1pi[MAX_DATA], 
  e0pr[MAX_DATA],  e0pi[MAX_DATA],  l0pr[MAX_DATA],  l0pi[MAX_DATA],
  m1mr[MAX_DATA],  m1mi[MAX_DATA],  l1mr[MAX_DATA],  l1mi[MAX_DATA], 
  f1r[MAX_DATA], f1i[MAX_DATA], f2r[MAX_DATA], f2i[MAX_DATA], f3r[MAX_DATA], 
  f3i[MAX_DATA], f4r[MAX_DATA], f4i[MAX_DATA], f5r[MAX_DATA], f5i[MAX_DATA], 
  f6r[MAX_DATA], f6i[MAX_DATA], t_pi_rad[MAX_DATA], ph_pi_rad[MAX_DATA], omega_cm[MAX_DATA],
  e1p3r[MAX_DATA], e1p3i[MAX_DATA], m1p3r[MAX_DATA], m1p3i[MAX_DATA],
  l1p3r[MAX_DATA], l1p3i[MAX_DATA],  e0p3r[MAX_DATA], e0p3i[MAX_DATA],
  l0p3r[MAX_DATA], l0p3i[MAX_DATA],  m1m3r[MAX_DATA], m1m3i[MAX_DATA],
  l1m3r[MAX_DATA], l1m3i[MAX_DATA],  e0p1r[MAX_DATA], e0p1i[MAX_DATA],
  l0p1r[MAX_DATA], l0p1i[MAX_DATA],  e1p1r[MAX_DATA], e1p1i[MAX_DATA],
  m1p1r[MAX_DATA], m1p1i[MAX_DATA],  l1p1r[MAX_DATA], l1p1i[MAX_DATA],
  m1m1r[MAX_DATA], m1m1i[MAX_DATA],  l1m1r[MAX_DATA], l1m1i[MAX_DATA];

Double_t Theory[MAX_DATA], Delta[MAX_DATA];
  
Double_t err_mat[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
Double_t err_mat3[3][3]={{0,0,0},{0,0,0},{0,0,0}};

/* in Old Data Format */
Double_t  m1m0r[MAX_DATA], m1m0i[MAX_DATA], l1m0r[MAX_DATA], l1m0i[MAX_DATA]; 
//Double_t  m1m1r[MAX_DATA], m1m1i[MAX_DATA], l1m1r[MAX_DATA], l1m1i[MAX_DATA]; 



#endif /* FIT_DATA_H */

