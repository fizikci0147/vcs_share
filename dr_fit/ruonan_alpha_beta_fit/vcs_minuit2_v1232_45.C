#include <iostream>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TVector.h"
#include "TMath.h"
#include <fstream>
//#include "constant.h"
#include "dr_func.h"
//const double d2r =TMath::DegToRad();
//const double r2d =TMath::RadToDeg();
//const double Mp=0.93827231;//GeV
R__LOAD_LIBRARY(libVCS.dylib) 

    using namespace std;

//    //calling DR code
//    extern "C"
//{
//    double v_compton_diff_cross_driver_dr_(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *,double *, double *, double *);
//}
//read in data
void ImportInput(vector<double> &gamma,vector<double> &sig0,vector<double> &sig0E, vector<double> &sig180, vector<double> &sig180E,const char*filename )
{
    double igamma,isig0,isig0E,isig180,isig180E;
    ifstream infile(filename);
    if(infile.fail()){
        cout<<filename <<" input file doesn't exist!"<<endl;
        exit(1);
    }else{
        while(!infile.eof()){
            infile>>igamma>>isig0>>isig0E >>isig180>>isig180E;
            gamma.push_back(igamma);
            sig0.push_back(isig0);
            sig0E.push_back(isig0E);
            sig180.push_back(isig180);
            sig180E.push_back(isig180E);
        }
        infile.close();
        gamma.pop_back();
        sig0.pop_back();
        sig0E.pop_back();
        sig180.pop_back();
        sig180E.pop_back();
    }
}

////calc E' and th for the new W
//void Calc_ep_eth(double W, double Q2, double ebeam, double &eth, double &ep)
//{
//
//    ep=ebeam-(W*W+Q2-Mp*Mp)/(2.*Mp);
//    eth=2.*asin(sqrt(Q2/(4.*ebeam*ep)));
//
//}

//import xsec with different W
void ImportInputWithKinematics(vector<double> &ebeam, vector<double> &e_mom, vector<double> &eth, vector<double> &thgg,vector<double> &W, 
        vector<double> &Q2, vector<double> &phgg1,vector<double> &phgg2,vector<double> &sig0,vector<double> &sig0E, vector<double> &sig180, vector<double> &sig180E,const char*filename )
{
    double iebeam,ie_mom,ieth,ithgg,iW,iQ2,iphgg1,iphgg2,isig0,isig0E,isig180,isig180E;
    ifstream infile(filename);
    if(infile.fail()){
        cout<<filename <<" input file doesn't exist!"<<endl;
        exit(1);
    }else{
        while(!infile.eof()){
            infile>>iebeam >>ithgg>>iW >> iQ2 >>iphgg1>>iphgg2 >>isig0>>isig0E >>isig180>>isig180E;
            //calc ep and eth
            Calc_ep_eth(iW,iQ2,iebeam,ieth,ie_mom);
       
            ebeam.push_back(iebeam);
            e_mom.push_back(ie_mom);
            eth.push_back(ieth);
            thgg.push_back(ithgg);
            W.push_back(iW);
            Q2.push_back(iQ2);
            phgg1.push_back(iphgg1);
            phgg2.push_back(iphgg2);
            sig0.push_back(isig0);
            //sig0E.push_back(isig0E);
            sig0E.push_back(isig0*0.1);
            sig180.push_back(isig180);
            //sig180E.push_back(isig180E);
            sig180E.push_back(isig180*0.1);
            //sig180E.push_back(isig180E);
        }
        infile.close();
        ebeam.pop_back();
        e_mom.pop_back();
        eth.pop_back();
        thgg.pop_back();
        W.pop_back();
        Q2.pop_back();
        phgg1.pop_back();
        phgg2.pop_back();
        sig0.pop_back();
        sig0E.pop_back();
        sig180.pop_back();
        sig180E.pop_back();
    }
}


//calculate Assymetries from xs
void CalcAssym(double sig0, double sig0E, double sig180, double sig180E, double &Assym, double &AssymE)
{
    double denom,numer,denomE, numerE;
    numer =(sig0-sig180);
    denom =(sig0+sig180);
    numerE =sqrt(pow(sig0E,2)+pow(sig180E,2));
    denomE =numerE;
    Assym  =100.*numer/denom;
    AssymE =Assym*sqrt(pow(numerE,2)/pow(numer,2)+pow(denomE,2)/pow(denom,2));

}
//calculate Chi2 to minimize
double GetChi2(const double *xx)
{
    double la   = xx[0];
    double lb   = xx[1];
    double S1   = xx[2];
    double M1   = xx[3];
    double E1   = xx[4];
    double chi2=0.;

    int opt=1;//1 for reading different W and 0 for same W
    vector<double> kin,kout,kth,thgg,W,Q2,phgg1,phgg2,Datasig0,Datasig0E,Datasig180,Datasig180E;
    //vector<double> thgg,Datasig0,Datasig0E,Datasig180,Datasig180E;
    int model =5,ffact=4;
    //double phcm1=0.,phcm2=TMath::Pi();

    //ImportInput(thgg,Datasig0,Datasig0E,Datasig180,Datasig180E,Form("vcs_hallc_data.txt")); //s1=0,  a=.86,b=0.0035
    //ImportInputWithKinematics(kin,kout,kth,thgg,W,Q2,phgg1,phgg2,Datasig0,Datasig0E,Datasig180,Datasig180E,Form("input/prop_007.txt")); //s1=0,  a=.86,b=0.0035
    //ImportInputWithKinematics(kin,kout,kth,thgg,W,Q2,phgg1,phgg2,Datasig0,Datasig0E,Datasig180,Datasig180E,Form("input/prop_015.txt")); //s1=0,  a=.86,b=0.0035
    //ImportInputWithKinematics(kin,kout,kth,thgg,W,Q2,phgg1,phgg2,Datasig0,Datasig0E,Datasig180,Datasig180E,Form("input/prop_2_2_005.txt")); //s1=0,  a=.86,b=0.0035
    //ImportInputWithKinematics(kin,kout,kth,thgg,W,Q2,phgg1,phgg2,Datasig0,Datasig0E,Datasig180,Datasig180E,Form("input/prop_1_1_005_v2t.txt")); //s1=0,  a=.86,b=0.0035
    //ImportInputWithKinematics(kin,kout,kth,thgg,W,Q2,phgg1,phgg2,Datasig0,Datasig0E,Datasig180,Datasig180E,Form("input/prop_2_2_0_15_v2t.txt")); //s1=0,  a=.86,b=0.0035
    ImportInputWithKinematics(kin,kout,kth,thgg,W,Q2,phgg1,phgg2,Datasig0,Datasig0E,Datasig180,Datasig180E,Form("input/prop_2_2_045_t.txt")); //s1=0,  a=.86,b=0.0035
    const int Ndata=thgg.size();
    double DataAssym[Ndata],DataAssymE[Ndata];
//    double kin=4.55,kout =4.0344,kth=7.68*d2r,q2=0.33;

    for(int i=0;i<Ndata;i++){
        //cout<<kin[i]<<" "<<kout[i]<<" "<<kth[i]*r2d<<" " <<thgg[i]<<" "<<phgg1[i] << " "<<phgg2[i]<<" "<<Q2[i]<<endl;
        double thggr =thgg[i]*d2r;
        double phcm1 =phgg1[i]*d2r;
        double phcm2 =phgg2[i]*d2r;

        //kin[i]=4.55,kout[i] =4.0344,kth[i]=7.68*d2r,Q2[i]=0.33;

        CalcAssym(Datasig0[i],Datasig0E[i],Datasig180[i],Datasig180E[i],DataAssym[i],DataAssymE[i]);
        double Modsig0  =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin[i],&kout[i],&kth[i],&thggr,&phcm1,&Q2[i],&S1,&M1,&E1);
        double Modsig180 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin[i],&kout[i],&kth[i],&thggr,&phcm2,&Q2[i],&S1,&M1,&E1);
        double ModAssym = 100.*(Modsig0-Modsig180)/(Modsig0+Modsig180); 
   //     cout<<" la " <<la<<" lb " <<lb<<" Datxs0 "<< Datasig0[i]<<" Mod xs0 "<<Modsig0<<" Datx180 "<<Datasig180[i]<<" Mod xs180 "<<Modsig180 <<" DatAssy" <<DataAssym[i] <<" ModAssym "<<ModAssym  << endl;

        if(Datasig0[i]>0){
        chi2 += TMath::Power(((Datasig0[i]-Modsig0)/Datasig0E[i]),2);
        }
        if(Datasig180[i]>0){
        chi2 += TMath::Power(((Datasig180[i]-Modsig180)/Datasig180E[i]),2);
        }
        if(abs(phgg1[i]-phgg2[i])==180.&&Datasig0[i]>0.&&Datasig180[i]>0.){
        //chi2 += TMath::Power(((DataAssym[i]-ModAssym)/DataAssymE[i]),2);
        }

    }
    cout <<chi2<<" "<<la<<" "<<lb<<endl;
    return chi2;
}
//main function
int vcs_minuit2_v1232_45(const char * minName ="Minuit2",const char *algoName="")
{

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    // set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);  // for GSL
    min->SetTolerance(0.01);
    min->SetPrintLevel(1);

    // create funciton wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&GetChi2,2);
    min->SetFunction(f);

    // the number of parameters
    const int npar = 5;
    Int_t ierflg = 0;
    double par[npar] ={0.7,0.7,0.711,1.0,1.0};               // the start values
    //double par[npar] ={0.7,0.75,0.711,1.0,1.0};               // the start values
    double stepSize[npar]={0.001,0.001,0.001,0.001,0.001};          // step sizes 
    //double minVal[npar]={0.,0.72,0.,0.,0.};            // minimum bound on parameter 
    //double minVal[npar]={0.,0.5,0.,0.,0.};            // minimum bound on parameter 
    double minVal[npar]={0.,0.,0.,0.,0.};            // minimum bound on parameter 
    //double maxVal[npar]={1.,0.78,1.,1.,1.};            // maximum bound on parameter
    double maxVal[npar]={2.,1.0,1.,1.,1.};            // maximum bound on parameter
    //double maxVal[npar]={2.,1.,1.,1.,1.};            // maximum bound on parameter
    string parName[npar]={"la","lb","S1","M1","E1"};

    //set variables
    for (int i=0; i<npar; i++){
        min->SetVariable(i, parName[i].c_str(), par[i], stepSize[i]);
        //min->SetVariableLimits(i,minVal[i],maxVal[i]);
    }

    //limit lb to get reasonable uncertainty 
    min->SetVariableLimits(1,minVal[1],maxVal[1]);
    //S1,M1 and E1 are fixed.
    min->FixVariable(2);
    min->FixVariable(3);
    min->FixVariable(4);
    min->SetPrecision(0.00001);

    // do the minimization
    min->Minimize();

    const double *xs = min->X();
    std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "<< min->MinValue()  << std::endl;



    //double W=1.232,la=0.7,lb=0.5, kin=4.55,kout ,kth,q2=0.33,thgg=140.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    //double W=1.230,la=0.,lb=0., kin=2.2,kout ,kth,q2=0.07;//thgg=140.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    double W=1.230,la=0.,lb=0., kin=2.2,kout ,kth,q2=0.45;//thgg=140.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    int model =5,ffact=4;
    double gp, gplab,Q20,dla,dlb;
    Calckinvars(W,q2,kin,kth,kout,gp,gplab,Q20);
    const double *DE = min->Errors();
    dla =DE[0];
    dlb =DE[1];
    la=xs[0];
    lb=xs[1];
    //double dla =min->GetMinosError(0), dlb=min->GetMinosError(1);
    cout<<"dla: " <<dla <<" dlb: "<<dlb<<endl;
    double la_low = la-dla, la_high =la+dla, lb_low=lb-dlb, lb_high =lb+dlb;

    double alpha =alpha_driver_dr_(&la,&Q20);
    double alpha_low =alpha_driver_dr_(&la_low,&Q20);
    double alpha_high =alpha_driver_dr_(&la_high,&Q20);
    double beta  =beta_driver_dr_(&lb,&Q20);
    double beta_low  =beta_driver_dr_(&lb_low,&Q20);
    double beta_high  =beta_driver_dr_(&lb_high,&Q20);
    double pll =p_ll_tot_driver_(&q2,&la,&ffact);
    double plt =p_lt_tot_driver_(&q2,&lb,&ffact);
    double ptt =p_tt_driver_(&q2,&ffact);


    std::cout <<" alpha: "<<Form("%.10e",alpha) <<" alpha low: " <<Form("%.10e",alpha_low) <<" alpha high: "<<Form("%.10e",alpha_high)<<" "
              <<"\n beta:  "<<Form("%.10e",beta)  <<" beta low: "  <<Form("%.10e",beta_low)  <<" beta high: " <<Form("%.10e",beta_high)<<endl;

    return 0;
}

