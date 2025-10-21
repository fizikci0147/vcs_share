#include <iostream>
#include "dr_func.h"
#include <iostream>
R__LOAD_LIBRARY(libVCS.dylib) 
using namespace std;

void run_alpha(){
    double W=1.230,la=0.7,lb=0.7, kin=2.2,kout=1.8 ,kth=11.2*d2r,q2=0.15,thgg=10.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    int model =5,ffact=4;
    double S1=0.711,M1=1.,E1=1.;

       double xs0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin,&kout,&kth,&thgg,&phcm0,&q2,&S1,&M1,&E1);
        TGraph *g0 =new TGraph();
        TGraphErrors *gp =new TGraphErrors();
        double step =0.01;
        double min =0.001;
        double max =1.0;
        const int N = int((max-min)/step);
    
        for(int i =0;i<N;i++){
            double alph = alpha_driver_dr_(&la,&min);
            g0->SetPoint(i,min,alph);
            cout <<"i :"<<i <<" Q2: "<<min<<" alpha: "<<alph<<endl;
            min +=step;
    }
   
        gp->SetPoint(0,0.07,g0->Eval(0.07));
        gp->SetPoint(1,0.1,g0->Eval(0.1));
        gp->SetPoint(2,0.15,g0->Eval(0.15));

        gp->SetPointError(0,0.0,g0->Eval(0.07)*0.1);
        gp->SetPointError(1,0.,g0->Eval(0.07)*0.1);
        gp->SetPointError(2,0.,g0->Eval(0.07)*0.1);

        gp->SetMarkerColor(2);
        gp->SetLineColor(2);
        gp->SetMarkerStyle(20);



        TLegend *t =new TLegend(0.6,0.6,0.8,0.8);
        t->AddEntry(g0,"0","l");

    TCanvas *c =new TCanvas("c","c");
    g0->SetLineStyle(1);
    g0->Draw("al");
    g0->SetTitle("");
    g0->GetYaxis()->SetTitle("#alpha_{E} (fm^{3})");
    g0->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    gp->Draw("psame");
    g0->GetYaxis()->SetMaxDigits(2);
   // g0->GetYaxis()->SetRangeUser(0,1);
   // g180->Draw("same");
   // t->Draw("same");
    
    
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
