#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"

void check_L1(TString iname="event.root", TString oname="check_L1.pdf"){
    auto file = new TFile(iname);
    auto tree = (TTree*)file->Get("tree");
    //auto newfile = new TFile(oname,"recreate");

    //TH2F* h_charge_sh_channel[16];
    TH1F* h_charge_sh_chip[16][2];
    TH1F* h_time_chip[16][2];


    auto h_FEB_label_channel = new TH2F("FEB_lable_channel","Hits in FEB_label:channel",65,0,64,16,0,16);
    h_FEB_label_channel->GetXaxis()->SetTitle("channel");
    h_FEB_label_channel->GetYaxis()->SetTitle("FEB_label");
    tree->Draw("FEB_label:channel>>FEB_lable_channel");


    auto h_charge_time = new TH2F("charge_time","Charge_SH:t_min_ttrigg",80,-80,80,100,-10,40);
    h_charge_time->GetXaxis()->SetTitle("t_min_ttrigg");
    h_charge_time->GetYaxis()->SetTitle("Charge_SH");
    tree->Draw("charge_SH:t_min_ttrigg>>charge_time","charge_SH>-10&&charge_SH<40&&t_min_ttrigg>-80&&t_min_ttrigg<80&&channel!=20&&layer==1");

    //for(int i=0; i<16; i++){
    //    TString t_name1 = Form("charge_sh_channel_%d",i);
    //    TString t_name2 = Form("FEB_label_%d",i);
    //    h_charge_sh_channel[i] = new TH2F(t_name1,t_name2,65,0,64,100,-25,25);
    //    TString name0 = Form("charge_SH:channel>>charge_sh_channel_%d",i);
    //    TString name1 = Form("FEB_label==%d&&charge_SH>-25&&charge_SH<25",i);
    //    h_charge_sh_channel[i]->GetXaxis()->SetTitle("channel");
    //    h_charge_sh_channel[i]->GetYaxis()->SetTitle("Charge_SH");
    //    tree->Draw(name0,name1);
    //}

    
    for(int i=0; i<16; i++){
	for(int j=0; j<2; j++){
	    TString t_name = Form("FEB_label_%d_chip_%d",i,j+1);
	    h_charge_sh_chip[i][j] = new TH1F(t_name,t_name,100,-10,40);
	    TString name0 = Form("charge_SH>>FEB_label_%d_chip_%d",i,j+1);
	    TString name1 = Form("FEB_label==%d&&chip==%d&&charge_SH>-10&&charge_SH<40&&t_min_ttrigg>-80&&t_min_ttrigg<80&&channel!=20&&layer==1",i,j+1);
	    h_charge_sh_chip[i][j]->GetXaxis()->SetTitle("Charge_SH");
	    h_charge_sh_chip[i][j]->GetYaxis()->SetTitle("Entries");
	    tree->Draw(name0,name1);
	}
    }


    for(int i=0; i<16; i++){
	for(int j=0; j<2; j++){
	    TString t_name = Form("time_FEB_label_%d_chip_%d",i,j+1);
	    h_time_chip[i][j] = new TH1F(t_name,t_name,80,-80,80);
	    TString name0 = Form("t_min_ttrigg>>time_FEB_label_%d_chip_%d",i,j+1);
	    TString name1 = Form("FEB_label==%d&&chip==%d&&charge_SH>-10&&charge_SH<40&&t_min_ttrigg>-80&&t_min_ttrigg<80&&channel!=20&&layer==1",i,j+1);
	    h_time_chip[i][j]->GetXaxis()->SetTitle("t_min_ttrigg");
	    h_time_chip[i][j]->GetYaxis()->SetTitle("Entries");
	    tree->Draw(name0,name1);
	}
    }

    gStyle->SetOptStat(0);
    auto can_feb_label_channel = new TCanvas("FEB_label_channel","FEB_label_channel",1200,800);
    can_feb_label_channel->cd();
    h_FEB_label_channel->Draw("colz");

    can_feb_label_channel->Print(oname+"(");
    //can_feb_label_channel->Write();


    auto can_charge_time = new TCanvas("charge_time","charge_time",1200,800);
    can_charge_time->cd();
    h_charge_time->Draw("colz");

    can_charge_time->Print(oname);
    //can_charge_time->Write();

    //auto can_charge_sh_channel = new TCanvas("Charge_SH_channel","Charge_SH_channel",1600,800);
    //can_charge_sh_channel->Divide(4,7);
    //for(int i=0; i<16; i++){
    //    can_charge_sh_channel->cd(i+1);
    //    h_charge_sh_channel[i]->Draw();
    //}
    //can_charge_sh_channel->Print("check_L2.pdf");
    //can_charge_sh_channel->Write();


    gStyle->SetOptStat(1);


    auto can_charge = new TCanvas("charge","charge",1200,800);
    can_charge->cd();
    auto h_charge = new TH1F("h_charge","charge_SH", 100, -10, 40);
    h_charge->GetXaxis()->SetTitle("Charge_SH");
    h_charge->GetYaxis()->SetTitle("Hits");
    tree->Draw("charge_SH>>h_charge","charge_SH>-10&&charge_SH<40&&t_min_ttrigg>-80&&t_min_ttrigg<80&&channel!=20&&layer==1");
    //tree->Draw("charge_SH>>h_charge","charge_SH>-10&&charge_SH<40&&channel!=20");

    can_charge->Print(oname);
    //can_charge->Write();

    auto can_time = new TCanvas("time","time",1200,800);
    can_time->cd();
    auto h_time = new TH1F("h_time","h_time", 80, -80, 80);
    h_time->GetXaxis()->SetTitle("t_min_ttrigg");
    h_time->GetYaxis()->SetTitle("Hits");
    tree->Draw("t_min_ttrigg>>h_time","charge_SH>-10&&charge_SH<40&&t_min_ttrigg>-80&&t_min_ttrigg<80&&channel!=20&&layer==1");
    //tree->Draw("charge_SH>>h_charge","charge_SH>-10&&charge_SH<40&&channel!=20");

    can_time->Print(oname);
    //can_time->Write();


    auto can_charge_sh_chip1 = new TCanvas("FEB_label_16-29_chip","FEB_label_0-7_chip",1600,800);
    can_charge_sh_chip1->Divide(4,4);
    for(int i=0; i<16; i++){
	can_charge_sh_chip1->cd(i+1);
	h_charge_sh_chip[i/2][i%2]->Draw();
    }
    can_charge_sh_chip1->Print(oname);
    //can_charge_sh_chip1->Write();

    auto can_charge_sh_chip2 = new TCanvas("FEB_label_30-43_chip","FEB_label_8-15_chip",1600,800);
    can_charge_sh_chip2->Divide(4,4);
    for(int i=16; i<32; i++){
	can_charge_sh_chip2->cd(i-15);
	h_charge_sh_chip[i/2][i%2]->Draw();
    }
    can_charge_sh_chip2->Print(oname);
    //can_charge_sh_chip2->Write();


    auto can_time_chip1 = new TCanvas("time_FEB_label_16-29_chip","time_FEB_label_0-7_chip",1600,800);
    can_time_chip1->Divide(4,4);
    for(int i=0; i<16; i++){
	can_time_chip1->cd(i+1);
	h_time_chip[i/2][i%2]->Draw();
    }
    can_time_chip1->Print(oname);
    //can_time_chip1->Write();

    auto can_time_chip2 = new TCanvas("time_FEB_label_30-43_chip","time_FEB_label_8-15_chip",1600,800);
    can_time_chip2->Divide(4,4);
    for(int i=16; i<32; i++){
	can_time_chip2->cd(i-15);
	h_time_chip[i/2][i%2]->Draw();
    }
    can_time_chip2->Print(oname+")");
    //can_time_chip2->Write();


    //newfile->Write();

}
