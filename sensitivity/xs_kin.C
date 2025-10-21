#include <iostream>
#include "dr_func.h"
#include <iostream>
R__LOAD_LIBRARY(libVCS.dylib) 
using namespace std;

void xs_kin(){
/*
 * Kin I
 */
    double la=0.7,lb=0.7, kin=1.1,kout[5] ,kth[5],q2=0.05,thgg=0.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    //double la=0.7,lb=0.7, kin=2.2,kout[5] ,kth[5],q2=0.07,thgg=0.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
//    double la=0.7,lb=0.7, kin=2.2,kout[5] ,kth[5],q2=0.15,thgg=0.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    double la1=0.1,la2=1.3,lb1=0.1,lb2=1.3;
    int model =5,ffact=4;
    double W[5]={1.210,1.22,1.23,1.24,1.25};
    double gp, gplab,Q20;
   // kout=1.82562;
   // kth = 7.57*d2r;
    phcm0=0*d2r;
    thgg =0*d2r;
   // cout<<"result is : " <<alpha_driver_dr_(&la,&qq2[0]) <<alpha_driver_dr_(&la,&qq2[1])<<endl;
    double S1=0.711,M1=1.,E1=1.;

 /*   double W=1.230,la=0.7,lb=0.7, kin=2.2,kout ,kth,q2=0.15,thgg=0.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    double la1=0.1,la2=1.3,lb1=0.1,lb2=1.3;
    int model =5,ffact=4;
    double gp, gplab,Q20;
    Calckinvars(W,q2,kin,kth,kout,gp,gplab,Q20);
    kout=1.78299;
    kth = 11.22*d2r;
    phcm0=0*d2r;
    thgg =0*d2r;
    cout <<" W: "<<W<<" Q2: "<<q2<<" kin: "<<kin<<" kout: "<<kout<<" kth: "<<kth*r2d<< "thgg: "<<thgg*r2d<<endl;
   // cout<<"result is : " <<alpha_driver_dr_(&la,&qq2[0]) <<alpha_driver_dr_(&la,&qq2[1])<<endl;
    double S1=0.711,M1=1.,E1=1.;
*/

        TGraph *gW_0[5];
        TGraph *gW_pi[5];
        TGraphErrors *gp_0[5];
        TGraphErrors *gp_pi[5];



    for(int k =0;k<5;k++){

        double step =2.;
        double min =100;
        double max =180;
        const int N = int((max-min)/step);

        gW_0[k]=new TGraph();
        gW_pi[k]=new TGraph();
        gp_0[k]  =new TGraphErrors();
        gp_pi[k]=new TGraphErrors();

        gp_0[k]->SetMarkerStyle(20); 
        gp_pi[k]->SetMarkerStyle(20);
        gp_0[k]->SetMarkerColor(2);
        gp_pi[k]->SetMarkerColor(2);
        gp_0[k]->SetLineColor(2);
        gp_pi[k]->SetLineColor(2);



        for(int i =0;i<N;i++){
         //   cout <<" i: "<<i<<" thgg: "<<min<<endl;
    Calckinvars(W[k],q2,kin,kth[k],kout[k],gp,gplab,Q20);
    //cout <<" W: "<<W[k]<<" Q2: "<<q2<<" kin: "<<kin<<" kout: "<<kout[k]<<" kth: "<<kth[k]*r2d<< "thgg: "<<min<<endl;
    //kth[k]*=d2r;
            double thst =min*d2r;
       double xs0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin,&kout[k],&kth[k],&thst,&phcm0,&q2,&S1,&M1,&E1);
       double xspi =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin,&kout[k],&kth[k],&thst,&phcmpi,&q2,&S1,&M1,&E1);


       double asy =(xs0-xspi)/(xs0+xspi)*100.;

            gW_0[k]->SetPoint(i,min,xs0);
            gW_pi[k]->SetPoint(i,min,xspi);
    

            min +=step;
        }
            gp_0[k]->SetPoint(0,110,gW_0[k]->Eval(110.));
            gp_0[k]->SetPointError(0,0,gW_0[k]->Eval(110.)*0.03);
            gp_pi[k]->SetPoint(0,110,gW_pi[k]->Eval(110.));
            gp_pi[k]->SetPointError(0,0,gW_pi[k]->Eval(110.)*0.03);

            cout <<kin<<"\t"<<110. <<"\t"<<W[k]<<" "<<q2<< 0<<"\t"<<180<<"\t"<<gW_0[k]->Eval(110.) <<"\t"<<gW_0[k]->Eval(110.)*0.03<<"\t"<<gW_pi[k]->Eval(110.)<<"\t"<<gW_pi[k]->Eval(110.)*0.03<<endl; 
            cout <<kin<<"\t"<<135. <<"\t"<<W[k]<<" "<<q2<< 0<<"\t"<<180<<"\t"<<gW_0[k]->Eval(135.) <<"\t"<<gW_0[k]->Eval(135.)*0.03<<"\t"<<gW_pi[k]->Eval(135.)<<"\t"<<gW_pi[k]->Eval(135.)*0.03<<endl; 
            cout <<kin<<"\t"<<160. <<"\t"<<W[k]<<" "<<q2<< 0<<"\t"<<180<<"\t"<<gW_0[k]->Eval(160.) <<"\t"<<gW_0[k]->Eval(160.)*0.03<<"\t"<<gW_pi[k]->Eval(160.)<<"\t"<<gW_pi[k]->Eval(160.)*0.03<<endl; 

            gp_0[k]->SetPoint(1,135,gW_0[k]->Eval(135.));
            gp_0[k]->SetPointError(1,0,gW_0[k]->Eval(135.)*0.03);
            gp_pi[k]->SetPoint(1,135,gW_pi[k]->Eval(135.));
            gp_pi[k]->SetPointError(1,0,gW_pi[k]->Eval(135.)*0.03);


            gp_0[k]->SetPoint(2,160,gW_0[k]->Eval(160.));
            gp_0[k]->SetPointError(2,0,gW_0[k]->Eval(160.)*0.03);
            gp_pi[k]->SetPoint(2,160,gW_pi[k]->Eval(160.));
            gp_pi[k]->SetPointError(2,0,gW_pi[k]->Eval(160.)*0.03);


    }


        TLegend *t =new TLegend(0.6,0.6,0.8,0.8);
        t->AddEntry(gW_0[0]," ","l");
        t->SetBorderSize(0);


        TCanvas *c =new TCanvas("c","c");
        c->Divide(5,2);

        for(int i=0;i<5;i++){
        c->cd(1+i);
        gW_0[i]->Draw("al");
        gp_0[i]->Draw("psame");
        gW_0[i]->SetTitle("0");
        gW_0[i]->GetYaxis()->SetTitle("d^{5}#sigma (nb/GeV^{2}/sr^{2})");
        gW_0[i]->GetXaxis()->SetTitle("#theta_{#gamma^{*} #gamma} (Deg)");
        gW_0[i]->GetYaxis()->SetRangeUser(2000,10000);
        gW_0[i]->GetYaxis()->SetMaxDigits(3);
        gW_0[i]->GetXaxis()->SetNdivisions(505);

        c->cd(6+i);
        gW_pi[i]->Draw("al");
        gp_pi[i]->Draw("psame");
        gW_pi[i]->SetTitle("180");
        gW_pi[i]->GetYaxis()->SetTitle("d^{5}#sigma (nb/GeV^{2}/sr^{2})");
        gW_pi[i]->GetXaxis()->SetTitle("#theta_{#gamma^{*} #gamma} (Deg)");
        gW_pi[i]->GetYaxis()->SetRangeUser(2000,10000);
        gW_pi[i]->GetYaxis()->SetMaxDigits(3);
        gW_pi[i]->GetXaxis()->SetNdivisions(505);

        }


}
