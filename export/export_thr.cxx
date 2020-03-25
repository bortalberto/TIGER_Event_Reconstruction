#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <string>
void export_thr_L1(int run_num)

{
    TString iname = "thr_run_";
    iname += Form("%i", run_num);
    iname += ".root";
    auto datafile = new TFile(iname);
    auto datatree = (TTree*)datafile->Get("tree");
    datatree->SetMarkerColor(2);
    datatree->SetMarkerStyle(20);
    datatree->SetMarkerSize(0.5);
    TImage *img = TImage::Create();


    TCanvas * c3 = new TCanvas("c3","Thresholds layer 1, strip X",1000,600);
    gStyle->SetOptStat(0);
    TH2F *h3 = new TH2F("h3", "Thresholds Layer 1 (strip X) ", 9000, 0.0, 900.0, 500, -2.0, 25.0);
    TH2F *h32 = new TH2F("h32", "Thresholds Layer 1 (strip X) ", 9000, 0.0, 900.0, 500, -2.0, 25.0);

    datatree->Draw("thr_t_fC:strip_x>>h3", "layer_id==1 && strip_x>0 ");
    datatree->Draw("thr_e_fC:strip_x>>h32", "layer_id==1 && strip_x>0 ","same");

    h3->GetXaxis()->SetTitle("Strip X");
    h3->GetYaxis()->SetTitle("Thresholds level [fC]");
    h3->SetMarkerColor(1);
    h3->SetMarkerStyle(20);
    h3->SetMarkerSize(0.5);
    h32->SetMarkerColor(2);
    h32->SetMarkerStyle(20);
    h32->SetMarkerSize(0.5);
    img->FromPad(c3);
    img->WriteImage("L1_X.png");

    TCanvas * c4 = new TCanvas("c4","Thresholds Layer 1, strip V",1000,600);
    gStyle->SetOptStat(0);
    TH2F *h4 = new TH2F("h4", "Thresholds Layer 1 (strip V) ", 12000, 0.0, 1200.0, 500, -2.0, 25.0);
    TH2F *h42 = new TH2F("h42", "Thresholds Layer 1 (strip V) ", 12000, 0.0, 1200.0, 500, -2.0, 25.0);

    datatree->Draw("thr_t_fC:strip_v>>h4", "layer_id==1 && strip_v>0 " );
    datatree->Draw("thr_e_fC:strip_v>>h42", "layer_id==1 && strip_v>0 ","same");

    h4->GetXaxis()->SetTitle("Strip V");
    h4->GetYaxis()->SetTitle("Thresholds level [fC]");
    h4->SetMarkerColor(1);
    h4->SetMarkerStyle(20);
    h4->SetMarkerSize(0.5);
    h42->SetMarkerColor(2);
    h42->SetMarkerStyle(20);
    h42->SetMarkerSize(0.5);
    img->FromPad(c4);
    img->WriteImage("L1_V.png");

    TCanvas * c5 = new TCanvas("c5","Thresholds layer 2 (strip V), hist",800,600);
    gStyle->SetOptStat(0);
    TH1F *h5 = new TH1F("h5", "Thresholds Layer 1 (strip V) ", 54, -2.0, 25.0);
    TH1F *h52 = new TH1F("h52", "Thresholds Layer 1 (strip V) ", 54, -2.0, 25.0);
    h5->SetLineColorAlpha(1, 0.35);
    h5->SetLineStyle(8);
    h52->SetLineColorAlpha(2, 0.35);
    h52->SetLineStyle(8);
    datatree->Draw("thr_t_fC>>h5", "strip_v>0 && layer_id==1");
    datatree->Draw("thr_e_fC>>h52", "strip_v>0 && layer_id==1");
    h5->Draw();
    h52->Draw("same");
    int max;
    if (h5->GetMaximum() > h52->GetMaximum()){
        max=h5->GetMaximum();
        }
    else{
        max=h52->GetMaximum();
        }

    h5->GetXaxis()->SetTitle("Thresholds level [fC]");
    h5->GetYaxis()->SetTitle("N");

    float m1,m2,rms1,rms2;
    auto legend8 = new TLegend(0.5,0.7,0.9,0.9);
    m1= h5->GetMean();  rms1= h5->GetRMS();
    m2= h52->GetMean();  rms2= h52->GetRMS();

    legend8->AddEntry(h5,Form("T branch, mu= %.2f, rms=%.2f",m1,rms1),"f");
    legend8->AddEntry(h52,Form("E branch, mu= %.2f, rms=%.2f",m2,rms2),"f");
    legend8->Draw();


    img->FromPad(c5);
    img->WriteImage("L1_V_hist.png");



    TCanvas * c6 = new TCanvas("c6","Thresholds Layer 1, strip X, hist",800,600);
    gStyle->SetOptStat(0);
    TH1F *h6 = new TH1F("h6", "Thresholds Layer 1 (strip X) ", 54, -2.0, 25.0);
    TH1F *h62 = new TH1F("h62", "Thresholds Layer 1 (strip X) ", 54, -2.0, 25.0);


    h6->SetLineColorAlpha(1, 0.35);
    h6->SetLineStyle(8);
    h62->SetLineColorAlpha(2, 0.35);
    h62->SetLineStyle(8);



    h6->GetXaxis()->SetTitle("Thresholds level [fC]");
    h6->GetYaxis()->SetTitle("N");
    datatree->Draw("thr_t_fC>>h6", "strip_x>0 && layer_id==1");
    datatree->Draw("thr_e_fC>>h62", "strip_x>0 && layer_id==1");
    h6->Draw();
    h62->Draw("same");


    auto legend9 = new TLegend(0.5,0.7,0.9,0.9);
    m1= h6->GetMean();  rms1= h6->GetRMS();
    m2= h62->GetMean();  rms2= h62->GetRMS();

    legend9->AddEntry(h6,Form("T branch, mu= %.2f, rms=%.2f",m1,rms1),"f");
    legend9->AddEntry(h62,Form("E branch, mu= %.2f, rms=%.2f",m2,rms2),"f");
    legend9->Draw();

    img->FromPad(c6);
    img->WriteImage("L1_X_hist.png");



}

void export_thr_L2(int run_num)

{
    TString iname = "thr_run_";
    iname += Form("%i", run_num);
    iname += ".root";
    auto datafile = new TFile(iname);
    auto datatree = (TTree*)datafile->Get("tree");
    datatree->SetMarkerColor(2);
    datatree->SetMarkerStyle(20);
    datatree->SetMarkerSize(0.5);
    TImage *img = TImage::Create();


    TCanvas * c3 = new TCanvas("c3","Thresholds layer 2, strip X",1000,600);
    gStyle->SetOptStat(0);
    TH2F *h3 = new TH2F("h3", "Thresholds Layer 2 (strip X) ", 13000, 0.0, 1300.0, 500, -2.0, 25.0);
    TH2F *h32 = new TH2F("h32", "Thresholds Layer 2 (strip X)", 13000, 0.0, 1300.0, 500, -2.0, 25.0);

    datatree->Draw("thr_t_fC:strip_x>>h3", "layer_id==2 && strip_x>0 ");
    datatree->Draw("thr_e_fC:strip_x>>h32", "layer_id==2 && strip_x>0 ","same");

    h3->GetXaxis()->SetTitle("Strip X");
    h3->GetYaxis()->SetTitle("Thresholds level [fC]");
    h3->SetMarkerColor(1);
    h3->SetMarkerStyle(20);
    h3->SetMarkerSize(0.5);
    h32->SetMarkerColor(2);
    h32->SetMarkerStyle(20);
    h32->SetMarkerSize(0.5);
    img->FromPad(c3);
    img->WriteImage("L2_X.png");

    TCanvas * c4 = new TCanvas("c4","Thresholds layer 2, strip V",1000,600);
    gStyle->SetOptStat(0);
    TH2F *h4 = new TH2F("h4", "Thresholds Layer 2 (strip V) ", 22000, 0.0, 2200.0, 500, -2.0, 25.0);
    TH2F *h42 = new TH2F("h42", "Thresholds Layer 2 (strip V) ", 22000, 0.0, 2200.0, 500, -2.0, 25.0);

    datatree->Draw("thr_t_fC:strip_v>>h4", "layer_id==2 && strip_v>0 " );
    datatree->Draw("thr_e_fC:strip_v>>h42", "layer_id==2 && strip_v>0 ","same");

    h4->GetXaxis()->SetTitle("Strip V");
    h4->GetYaxis()->SetTitle("Thresholds level [fC]");
    h4->SetMarkerColor(1);
    h4->SetMarkerStyle(20);
    h4->SetMarkerSize(0.5);
    h42->SetMarkerColor(2);
    h42->SetMarkerStyle(20);
    h42->SetMarkerSize(0.5);
    img->FromPad(c4);
    img->WriteImage("L2_V.png");

    TCanvas * c5 = new TCanvas("c5","Thresholds layer 2 (strip V), hist",800,600);
    gStyle->SetOptStat(0);
    TH1F *h5 = new TH1F("h5", "Thresholds Layer 2 (strip V) ", 54, -2.0, 25.0);
    TH1F *h52 = new TH1F("h52", "Thresholds Layer 2 (strip V) ", 54, -2.0, 25.0);
    h5->SetLineColorAlpha(1, 0.35);
    h5->SetLineStyle(8);
    h52->SetLineColorAlpha(2, 0.35);
    h52->SetLineStyle(8);
    datatree->Draw("thr_t_fC>>h5", "strip_v>0 && layer_id==2");
    datatree->Draw("thr_e_fC>>h52", "strip_v>0 && layer_id==2");
    h5->Draw();
    h52->Draw("same");
    int max;
    if (h5->GetMaximum() > h52->GetMaximum()){
        max=h5->GetMaximum();
        }
    else{
        max=h52->GetMaximum();
        }

    h5->GetXaxis()->SetTitle("Thresholds level [fC]");
    h5->GetYaxis()->SetTitle("N");

    float m1,m2,rms1,rms2;
    auto legend8 = new TLegend(0.5,0.7,0.9,0.9);
    m1= h5->GetMean();  rms1= h5->GetRMS();
    m2= h52->GetMean();  rms2= h52->GetRMS();

    legend8->AddEntry(h5,Form("T branch, mu= %.2f, rms=%.2f",m1,rms1),"f");
    legend8->AddEntry(h52,Form("E branch, mu= %.2f, rms=%.2f",m2,rms2),"f");
    legend8->Draw();


    img->FromPad(c5);
    img->WriteImage("L2_V_hist.png");



    TCanvas * c6 = new TCanvas("c6","Thresholds Layer 2, strip X, hist",800,600);
    gStyle->SetOptStat(0);
    TH1F *h6 = new TH1F("h6", "Thresholds layer 2 (strip X) ", 54, -2.0, 25.0);
    TH1F *h62 = new TH1F("h62", "Thresholds layer 2 (strip X)", 54, -2.0, 25.0);


    h6->SetLineColorAlpha(1, 0.35);
    h6->SetLineStyle(8);
    h62->SetLineColorAlpha(2, 0.35);
    h62->SetLineStyle(8);



    h6->GetXaxis()->SetTitle("Thresholds level [fC]");
    h6->GetYaxis()->SetTitle("N");
    datatree->Draw("thr_t_fC>>h6", "strip_x>0 && layer_id==2");
    datatree->Draw("thr_e_fC>>h62", "strip_x>0 && layer_id==2");
    h6->Draw();
    h62->Draw("same");


    auto legend9 = new TLegend(0.5,0.7,0.9,0.9);
    m1= h6->GetMean();  rms1= h6->GetRMS();
    m2= h62->GetMean();  rms2= h62->GetRMS();

    legend9->AddEntry(h6,Form("T branch, mu= %.2f, rms=%.2f",m1,rms1),"f");
    legend9->AddEntry(h62,Form("E branch, mu= %.2f, rms=%.2f",m2,rms2),"f");
    legend9->Draw();

    img->FromPad(c6);
    img->WriteImage("L2_X_hist.png");



}
