#include <iostream>
#include "dr_func.h"
#include <iostream>
R__LOAD_LIBRARY(libVCS.dylib) 
using namespace std;

void sensitivity_a_b_vs_W_nikos(){
/*
 * Kin I
 */
    double W=0,la=0.7,lb=0.7, kin=2.2,kout ,kth,q2=0.07,thgg=135.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
    double la1=0.1,la2=1.3,lb1=0.1,lb2=1.3;
    int model =5,ffact=4;



  //  phcm0=0*d2r;
  //  thgg =0*d2r;
    cout <<" W: "<<W<<" Q2: "<<q2<<" kin: "<<kin<<" kout: "<<kout<<" kth: "<<kth*r2d<< "thgg: "<<thgg*r2d<<endl;
   // cout<<"result is : " <<alpha_driver_dr_(&la,&qq2[0]) <<alpha_driver_dr_(&la,&qq2[1])<<endl;
    double S1=0.711,M1=1.,E1=1.;


    /*
     *KIn II
     *
    double W=1.230,la=0.7,lb=0.7, kin=2.2,kout ,kth,q2=0.15,thgg=0.*d2r,phcm0=0.*d2r,phcmpi=180.*d2r;
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
        TGraph *g0 =new TGraph();
        TGraph *g180 =new TGraph();

        TGraph *glb1_0 =new TGraph();
        TGraph *glb2_0 =new TGraph();
        TGraph *gla1_0 =new TGraph();
        TGraph *gla2_0 =new TGraph();

        TGraph *glb1_180 =new TGraph();
        TGraph *glb2_180 =new TGraph();
        TGraph *gla1_180 =new TGraph();
        TGraph *gla2_180 =new TGraph();
        gla2_0->SetLineStyle(7);
        glb2_0->SetLineStyle(7);
        gla2_180->SetLineStyle(7);
        glb2_180->SetLineStyle(7);

        TGraph *gasy =new TGraph();
        TGraph *gasya1 =new TGraph();
        TGraph *gasya2 =new TGraph();

        TGraph *gasyb1 =new TGraph();
        TGraph *gasyb2 =new TGraph();

        double step =0.001;
        double min =1.190;
        double max =1.260;
        const int N = int((max-min)/step);
    
        for(int i =0;i<N;i++){

    W=min;
    double gp, gplab,Q20;
    Calckinvars(W,q2,kin,kth,kout,gp,gplab,Q20);
    //kth *=d2r;



          //  cout <<" i: "<<i<<" thgg: "<<min<<endl;
            double thst =thgg;
    
            cout<<i <<" "<<kin<<" "<<kout<<" kth "<<kth <<" "<<thst <<endl;
       double xs0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin,&kout,&kth,&thst,&phcm0,&q2,&S1,&M1,&E1);
       double xspi =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb,&kin,&kout,&kth,&thst,&phcmpi,&q2,&S1,&M1,&E1);


       double xslb1_0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb1,&kin,&kout,&kth,&thst,&phcm0,&q2,&S1,&M1,&E1);
       double xslb2_0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb2,&kin,&kout,&kth,&thst,&phcm0,&q2,&S1,&M1,&E1);
       double xsla1_0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la1,&lb,&kin,&kout,&kth,&thst,&phcm0,&q2,&S1,&M1,&E1);
       double xsla2_0 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la2,&lb,&kin,&kout,&kth,&thst,&phcm0,&q2,&S1,&M1,&E1);


       double xslb1_180 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb1,&kin,&kout,&kth,&thst,&phcmpi,&q2,&S1,&M1,&E1);
       double xslb2_180 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la,&lb2,&kin,&kout,&kth,&thst,&phcmpi,&q2,&S1,&M1,&E1);
       double xsla1_180 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la1,&lb,&kin,&kout,&kth,&thst,&phcmpi,&q2,&S1,&M1,&E1);
       double xsla2_180 =v_compton_diff_cross_driver_dr_(&model,&ffact,&la2,&lb,&kin,&kout,&kth,&thst,&phcmpi,&q2,&S1,&M1,&E1);


       double asy =(xs0-xspi)/(xs0+xspi)*100.;
       double asya1 =(xsla1_0-xsla1_180)/(xsla1_0+xsla1_180)*100.;
       double asya2 =(xsla2_0-xsla2_180)/(xsla2_0+xsla2_180)*100.;
       double asyb1 =(xslb1_0-xslb1_180)/(xslb1_0+xslb1_180)*100.;
       double asyb2 =(xslb2_0-xslb2_180)/(xslb2_0+xslb2_180)*100.;



       if(xs0!=xs0||xspi!=xspi) continue;

            g0->SetPoint(i,min,xs0);
            g180->SetPoint(i,min,xspi);
    
            gla1_0->SetPoint(i,min,xsla1_0);
            gla2_0->SetPoint(i,min,xsla2_0);
            glb1_0->SetPoint(i,min,xslb1_0);
            glb2_0->SetPoint(i,min,xslb2_0);

            gla1_180->SetPoint(i,min,xsla1_180);
            gla2_180->SetPoint(i,min,xsla2_180);
            glb1_180->SetPoint(i,min,xslb1_180);
            glb2_180->SetPoint(i,min,xslb2_180);


            gasy->SetPoint(i,min,asy);
            gasya1->SetPoint(i,min,asya1);
            gasya2->SetPoint(i,min,asya2);
            gasyb1->SetPoint(i,min,asyb1);
            gasyb2->SetPoint(i,min,asyb2);




            min +=step;
        }



        gla1_0->SetLineColor(2);  
        gla2_0->SetLineColor(2);
        glb1_0->SetLineColor(4);
        glb2_0->SetLineColor(4);


        gla1_180->SetLineColor(2);  
        gla2_180->SetLineColor(2);
        glb1_180->SetLineColor(4);
        glb2_180->SetLineColor(4);

        gasya1->SetLineColor(2);
        gasya2->SetLineColor(2);
        gasyb1->SetLineColor(4);
        gasyb2->SetLineColor(4);

        gasya2->SetLineStyle(7);
        gasyb2->SetLineStyle(7);

        TLegend *t =new TLegend(0.6,0.6,0.8,0.8);
        //t->AddEntry(g0,"la=0.7,lb=0.7","l");
        t->AddEntry(gla1_0,"la=0.1,lb=0.7","l");
        t->AddEntry(gla2_0,"la=1.3,lb=0.7","l");
       // t->AddEntry(glb1_0,"la=0.7,lb=0.1","l");
       // t->AddEntry(glb2_0,"la=0.7,lb=1.3","l");
        t->SetBorderSize(0);


        TCanvas *c =new TCanvas("c","c");
        c->Divide(2,1);
        c->cd(1);
            g0->SetLineStyle(1);
        g180->SetLineStyle(1);
        //g180->SetLineColor(1);
        gla1_0->Draw("al");
        gla1_0->SetTitle("0");
        gla1_0->GetYaxis()->SetTitle("d^{5}#sigma (nb/GeV^{2}/sr^{2})");
        gla1_0->GetXaxis()->SetTitle("W (GeV)");
      //  g0->GetYaxis()->SetRangeUser(0,40000);
        t->Draw("same");
        gla1_0->Draw("lsame"); 
        gla2_0->Draw("lsame");
        //glb1_0->Draw("lsame");
        //glb2_0->Draw("lsame");
        gla1_0->GetYaxis()->SetMaxDigits(3);
        gla1_0->GetXaxis()->SetNdivisions(505);

        c->cd(2);
        gla1_180->Draw("al");
        gla1_180->SetTitle("180");
        gla1_180->GetYaxis()->SetTitle("d^{5}#sigma (nb/GeV^{2}/sr^{2})");
        gla1_180->GetXaxis()->SetTitle("W (GeV)");
       // g180->GetYaxis()->SetRangeUser(0,40000);
        gla1_180->GetXaxis()->SetNdivisions(505);

        gla1_180->Draw("lsame"); 
        gla2_180->Draw("lsame");
        //glb1_180->Draw("lsame");
        //glb2_180->Draw("lsame");
        gla1_180->GetYaxis()->SetMaxDigits(3);

       /* c->cd(3);
        gasy->Draw("al");
        gasy->SetTitle("");
        gasy->GetXaxis()->SetNdivisions(505);
        gasy->GetYaxis()->SetTitle("Assymetry %");
        gasy->GetXaxis()->SetTitle("#theta_{#gamma^{*} #gamma} (Deg)");
        gasya1->Draw("lsame");
        gasya2->Draw("lsame");
        gasyb1->Draw("lsame");
        gasyb2->Draw("lsame");
        */


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
