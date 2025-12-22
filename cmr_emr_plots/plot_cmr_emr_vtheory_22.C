#include "func.h"
//in this version removed the syst errors for Adams and propoposal points
//Sean's point is changed to published version
void plot_cmr_emr_vtheory_22(){

    vector<double> CQ2,cmr,CE1,CE2;
    vector<double> pCQ2,pcmr,pCE1,pCE2;
    vector<double> EQ2,emr,EE1,EE2;
    vector<double> pEQ2,pemr,pEE1,pEE2;
    vector<string> Cauth, Eauth;
    vector<string> pCauth, pEauth;
    //ImportData(CQ2,cmr,CE1,CE2,Cauth,Form("cmr_proposal.txt"));
    ImportData(CQ2,cmr,CE1,CE2,Cauth,Form("cmr.txt"));
    ImportData(pCQ2,pcmr,pCE1,pCE2,pCauth,Form("cmr_proposal_final_updated2.txt"));
    //ImportData(pCQ2,pcmr,pCE1,pCE2,pCauth,Form("cmr_proposal_final.txt"));
    ImportData(EQ2,emr,EE1,EE2,Eauth,Form("emr_all.txt"));
    ImportData(pEQ2,pemr,pEE1,pEE2,pEauth,Form("emr_proposal_updated2.txt"));

    TGraphErrors *gE =new TGraphErrors();
    TGraphErrors *gEp =new TGraphErrors();
    TGraphErrors *gC =new TGraphErrors();
    TGraphErrors *gCp =new TGraphErrors();

    gE->SetMarkerStyle(21);
    gEp->SetMarkerStyle(20);
    gEp->SetMarkerColor(2);
    gEp->SetLineColor(2);
    gC->SetMarkerStyle(21);
    gCp->SetMarkerStyle(20);
    gCp->SetMarkerColor(2);
    gCp->SetLineColor(2);

    TGraph *ghqm =new TGraph("theory_data/hqm.txt");
    TGraph *gcapstick =new TGraph("theory_data/capstick.txt");
    TGraph *gdmt =new TGraph("theory_data/dmt.txt");
    TGraph *gdsem =new TGraph("theory_data/dsem.txt");
    TGraph *ggh =new TGraph("theory_data/gh.txt");
    TGraph *glarge_ne =new TGraph("theory_data/large_ne.txt");
    TGraph *gmaid =new TGraph("theory_data/maid.txt");
    TGraph *gpv =new TGraph("theory_data/pv.txt");
    TGraph *gsaid =new TGraph("theory_data/said.txt");
    TGraph *gsato_lee =new TGraph("theory_data/sato_lee.txt");
    TGraph *gsato_lee_bare =new TGraph("theory_data/sato_lee_bare.txt");

    ghqm->SetLineColor(1);
    gcapstick->SetLineColor(kOrange-3);
    gdmt->SetLineColor(4);
    gdsem->SetLineColor(4);
    ggh->SetLineColor(kGray);
    glarge_ne->SetLineColor(kMagenta+3);
    gmaid->SetLineColor(1);
    gpv->SetLineColor(1);
    gsaid->SetLineColor(kGreen-2);
    gsato_lee->SetLineColor(2);
    gsato_lee_bare->SetLineColor(2);

    ghqm->SetLineStyle(3);
    gcapstick->SetLineStyle(8);
    gdmt->SetLineStyle(5); //maybe 10
    gdsem->SetLineStyle(1);
    ggh->SetLineStyle(2);
    glarge_ne->SetLineStyle(8);
    gmaid->SetLineStyle(7);
    gpv->SetLineStyle(4);
    gsaid->SetLineStyle(3);
    gsato_lee->SetLineStyle(1);
    gsato_lee_bare->SetLineStyle(7);

    TLegend *lg =new TLegend(0.3,0.3,0.8,0.8);
    lg->AddEntry(ghqm,"HQM","l");
    lg->AddEntry(gcapstick,"Capstick","l");
    lg->AddEntry(gdmt,"DMT","l");
    lg->AddEntry(gdsem,"DSEM","l");
    lg->AddEntry(ggh,"GH","l");
    lg->AddEntry(glarge_ne,"Large-Ne","l");
    lg->AddEntry(gmaid,"MAID","l");
    lg->AddEntry(gpv,"PV","l");
    lg->AddEntry(gsaid,"SAID","l");
    lg->AddEntry(gsato_lee,"Sato Lee","l");
    lg->AddEntry(gsato_lee_bare,"Sato Lee (bare)","l");

    lg->SetBorderSize(0);



    for(int i=0;i<emr.size();i++){
        if(i<3) EE2[i]=0.;
        double EE = sqrt(EE1[i]*EE1[i]+EE2[i]*EE2[i]);
        //double EE = sqrt(EE1[i]*EE1[i]+EE2[i]*EE2[i]);
        gE->SetPoint(i,EQ2[i],emr[i]);
        gE->SetPointError(i,0.,EE);
    }

    for(int i=0;i<cmr.size();i++){
        if(i<3) CE2[i]=0.;
        double CE = sqrt(CE1[i]*CE1[i]+CE2[i]*CE2[i]);
        //double CE = sqrt(CE1[i]*CE1[i]+CE2[i]*CE2[i]);
        gC->SetPoint(i,CQ2[i],-cmr[i]);
        gC->SetPointError(i,0.,CE);
    }

    for(int i=0;i<pemr.size();i++){
        //double pEE = sqrt(pEE1[i]*pEE1[i]+pEE2[i]*pEE2[i]);
        double pEE =pEE1[i];// sqrt(pEE1[i]*pEE1[i]+pEE2[i]*pEE2[i]);
        //double pEE = sqrt(pEE1[i]*pEE1[i]+pEE2[i]*pEE2[i]);
        gEp->SetPoint(i,pEQ2[i],pemr[i]);
        gEp->SetPointError(i,0.,pEE*0.92);//adding off plane xs decreased uncertainties by 8% 
    }

    for(int i=0;i<pcmr.size();i++){
        //double pCE =sqrt(pCE1[i]*pCE1[i]+pCE2[i]*pCE2[i]);
        double pCE =pCE1[i];//sqrt(pCE1[i]*pCE1[i]+pCE2[i]*pCE2[i]);
        //double pCE = sqrt(pCE1[i]*pCE1[i]+pCE2[i]*pCE2[i]);
        gCp->SetPoint(i,pCQ2[i],pcmr[i]);
        gCp->SetPointError(i,0.,pCE*0.9);//adding off plane xs decreased uncertainties by 10% 
    }
    gE->SetLineColor(4);
    gE->SetMarkerColor(4);
    gC->SetLineColor(4);
    gC->SetMarkerColor(4);


    auto *gE2 =(TGraphErrors*)gE->Clone();
    auto *gC2 =(TGraphErrors*)gC->Clone();

    TCanvas * cE= new TCanvas("cE","cE");

    gE->Draw("ap");
    gEp->Draw("psame");
    gE->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gE->GetYaxis()->SetTitle("EMR (%)");
    gE->SetTitle(" ");

    TCanvas * cC= new TCanvas("cC","cC");
    gC->Draw("ap");
    gCp->Draw("psame");
    gC->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gC->GetYaxis()->SetTitle("CMR (%)");
    gC->SetTitle(" ");
/*
    TCanvas * cb= new TCanvas("cb","cb");
    cb->Divide(2,2);
    cb->cd(1);
    gC->Draw("ap");
    gC->GetXaxis()->SetLimits(0,0.93);
    gC->GetYaxis()->SetRangeUser(-9.3,0);
    gCp->Draw("psame");
    gC->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gC->GetYaxis()->SetTitle("CMR (%)");
    gC->SetTitle(" ");

    cb->cd(3);

    gE->Draw("ap");
    gE->GetXaxis()->SetLimits(-0.01,0.93);
    gE->GetYaxis()->SetRangeUser(-6,0);
    gEp->Draw("psame");
    gE->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gE->GetYaxis()->SetTitle("EMR (%)");
    gE->SetTitle(" ");



    cb->cd(2);
    gC2->Draw("ap");
    gC2->GetXaxis()->SetLimits(0,0.26);
    gC2->GetYaxis()->SetRangeUser(-9.3,0);
    gCp->Draw("psame");
    gC2->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gC2->GetYaxis()->SetTitle("CMR (%)");
    gC2->SetTitle(" ");

    cb->cd(4);

    gE2->Draw("ap");
    gE2->GetXaxis()->SetLimits(-0.01,0.26);
    gE2->GetYaxis()->SetRangeUser(-6,0);
    gEp->Draw("psame");
    gE2->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gE2->GetYaxis()->SetTitle("EMR (%)");
    gE2->SetTitle(" ");
*/


    TCanvas * cb= new TCanvas("cb","cb");
    cb->Divide(1,1);
    cb->cd(1);
    gC->Draw("ap");
    gC->GetXaxis()->SetLimits(0,1.02);
    gC->GetYaxis()->SetRangeUser(-9.3,0.8);
    gC->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gC->GetYaxis()->SetTitle("CMR (%)");
    gC->SetTitle(" ");

    ghqm->Draw("lsame");
    gcapstick->Draw("lsame");
    gdmt->Draw("lsame");
    gdsem->Draw("lsame");
    ggh->Draw("lsame");
    glarge_ne->Draw("lsame");
    gmaid->Draw("lsame");
    gpv->Draw("lsame");
    gsaid->Draw("lsame");
    gsato_lee->Draw("lsame");
    gsato_lee_bare->Draw("lsame");

    gCp->Draw("psame");
    gC->Draw("psame");

    //lg->Draw("same");

  /*  cb->cd(2);

    gE->Draw("ap");
    gE->GetXaxis()->SetLimits(-0.01,1.02);
    gE->GetYaxis()->SetRangeUser(-6,0);
    gEp->Draw("psame");
    gE->GetXaxis()->SetTitle("Q^{2} (GeV/c)");
    gE->GetYaxis()->SetTitle("EMR (%)");
    gE->SetTitle(" ");
*/



}
