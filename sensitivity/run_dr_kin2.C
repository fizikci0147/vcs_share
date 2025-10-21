#include <iostream>
#include "dr_func.h"
#include <iostream>
R__LOAD_LIBRARY(libVCS.dylib) 
using namespace std;

void run_dr_kin2(){
    double W=1.230,la=0.7,lb=0.7, kin=2.2,kout ,kth,q2=0.15,thgg=0.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
//    gSystem->Load("object_template/vcompton_vcsdr07.o");    
//    gSystem->Load("../libVCS.dylib");    

//    double la =0.859;
   // double qq2[2] ={0.057,0.1};
   // double q20[2] ={Q20(qq2[0],1.232),Q20(qq2[1],1.232)};

    //cout <<q20[0] <<" "<<q20[1]<<endl;

  
  //  double walpha[2] ={7.85e-4,6.6e-4};
  //  double wdalpha[2] ={1.06e-4,1.09e-4};
   // cout<<"result is : " <<alpha_driver_dr_(&la,&qq2[0]) <<alpha_driver_dr_(&la,&qq2[1])<<endl;

    //double W=1.220,la=0.7,lb=0.5, kin=4.55,kout =4.03444,kth=7.68791*d2r,q2=0.33,thgg=150.*d2r,phcm0=30.*d2r,phcmpi=215.*d2r;
    int model =5,ffact=4;
    double gp, gplab,Q20;
    Calckinvars(W,q2,kin,kth,kout,gp,gplab,Q20);
    kout=1.78299;
    kth = 11.22*d2r;
    phcm0=0*d2r;
    thgg =0*d2r;
    cout <<" W: "<<W<<" Q2: "<<q2<<" kin: "<<kin<<" kout: "<<kout<<" kth: "<<kth*r2d<< "thgg: "<<thgg*r2d<<endl;
  //  cout<<"result is : " <<alpha_driver_dr_(&la,&qq2[0]) <<alpha_driver_dr_(&la,&qq2[1])<<endl;
    double S1=0.711,M1=1.,E1=1.;

        TGraph *g0 =new TGraph();
        TGraph *g180 =new TGraph();
        double step =2.;
        double min =0;
        double max =180;
        const int N = int((max-min)/step);
    
        for(int i =0;i<N;i++){
            cout <<"thgg: "<<min<<endl;
            double thst =min*d2r;
       double xs0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin,&kout,&kth,&thst,&phcm0,&q2,&S1,&M1,&E1);
       double xspi =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin,&kout,&kth,&thst,&phcmpi,&q2,&S1,&M1,&E1);

       if(xs0!=xs0||xspi!=xspi) continue;

            g0->SetPoint(i,min,xs0);
            g180->SetPoint(i,min,xspi);
    
            min +=step;
    }
   
        TLegend *t =new TLegend(0.6,0.6,0.8,0.8);
        t->AddEntry(g0,"0","l");
        t->AddEntry(g180,"180","l");


    TCanvas *c =new TCanvas("c","c");
    g0->SetLineStyle(1);
    g180->SetLineStyle(7);
    g180->SetLineColor(2);
    g0->Draw("al");
    g0->SetTitle("");
    g0->GetYaxis()->SetTitle("d^{5}#sigma (nb/GeV^{2}/sr^{2})");
    g0->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    g0->GetYaxis()->SetRangeUser(0,60000);
    g180->Draw("same");
    t->Draw("same");
    
    
    /*  
    double alpha =alpha_driver_dr_(&la,&q2);
    //double alpha =alpha_driver_dr_(&la,&Q20);
    double beta  =beta_driver_dr_(&lb,&Q20);
    double pll =p_ll_tot_driver_(&q2,&la,&ffact);
    double plt =p_lt_tot_driver_(&q2,&lb,&ffact);
    double ptt =p_tt_driver_(&q2,&ffact);

    cout<<"Q2 "<<q2<<" Q20: "<<Q20 <<" Ep: " <<kout<<" kth: "<<kth*r2d<<endl;
    cout<<"xs0 "  <<" xsE: "<<" xs180: " <<" xs180E: "<<" \n" 
     <<xs0 <<"  "<<0.0144*xs0<<" " <<xspi<<"  "<<0.0144*xspi 
        <<"\n"<<" alpha: " <<Form("%.10e",alpha)<<" beta: "  <<Form("%.10e",beta)  <<" ";


    cout <<" PLL: "<<pll<<" PLT: "<<plt <<" PTT: "<<ptt<<endl;

*/
}
