void GEMROC_Analyzer() {
  
  const int MAX_CHANNELS = 5;
  int run;
  
  std::cout << "\n***** Which RUN is this? *****" << std::endl;
  std::cin >> run;

  TCanvas *c1 = new TCanvas("c1", "multipads", 900, 700);
  gStyle->SetOptStat(0);

  TFile *f0 = TFile::Open("event.root","r");
  TTree *t0 = (TTree *)f0->Get("tree");
  
  TH1F *Diff_Tiger = new TH1F("Diff_Tiger", " ", 100, -1000, 1000);
  
  double avg_0_DT[88], sigma_0_DT[88], avg_1_DT[88], sigma_1_DT[88], avg_2_DT[88], sigma_2_DT[88], avg_3_DT[88], sigma_3_DT[88], avg_4_DT[88], sigma_4_DT[88], TID[88];

  TString what_to_plot, how_to_plot;

  for(int chID = 0; chID < MAX_CHANNELS; chID++) { // Looping over the i-th channel
    for(int i = 0; i < 88;i++) {                   // Looping over the i-th TIGER
      
      what_to_plot = Form("TP_diff[%d][%d]>>Diff_Tiger" , i, chID);
      if(i < 32) how_to_plot  = Form("TP_diff[0][%d] != -999"        ,    chID);
      else       how_to_plot  = Form("TP_diff[32][%d] != -999"       ,    chID);

      t0->Draw(what_to_plot, how_to_plot);

      if(chID == 0) {
	avg_0_DT[i]   = Diff_Tiger->GetMean();
	sigma_0_DT[i] = Diff_Tiger->GetStdDev();
      }
      else if(chID == 1) {
	avg_1_DT[i]   = Diff_Tiger->GetMean();
	sigma_1_DT[i] = Diff_Tiger->GetStdDev();
      }
      else if(chID == 2) {
	avg_2_DT[i]   = Diff_Tiger->GetMean();
	sigma_2_DT[i] = Diff_Tiger->GetStdDev();
      }
      else if(chID == 3) {
	avg_3_DT[i]   = Diff_Tiger->GetMean();
	sigma_3_DT[i] = Diff_Tiger->GetStdDev();
      }
      else if(chID == 4) {
	avg_4_DT[i]   = Diff_Tiger->GetMean();
	sigma_4_DT[i] = Diff_Tiger->GetStdDev();
      }

      if(chID == MAX_CHANNELS-1) TID[i] = i;

      Diff_Tiger->Reset();
    }
  }
  
  c1->Divide(1,2);

  c1->cd(1);
  TGraph* mean_f = new TGraph(88, TID, avg_3_DT);
  mean_f->SetMarkerStyle(28);
  mean_f->SetTitle("i-th TIGER #mu_{#Delta T} per channel");
  mean_f->GetXaxis()->SetTitle("TIGER ID");
  mean_f->GetXaxis()->SetRangeUser(0, 91);
  mean_f->GetYaxis()->SetTitle("#mu_{#Delta T}  [Clock Count]");
  mean_f->Draw("ap");

  TLegend* mean_leg = new TLegend(0.1, 0.36, 0.33, 0.56 );
  mean_leg->AddEntry(mean_f , "Channel 20", "p");
 
  if(run > 58) {   
    TGraph* mean_s = new TGraph(88, TID, avg_1_DT);
    mean_s->SetMarkerStyle(3);
    mean_s->SetMarkerColor(kRed);
    mean_s->Draw("p");

    mean_leg->AddEntry(mean_s , "Channel 10", "p");
  }

  if(run > 59) {
    TGraph* mean_t = new TGraph(88, TID, avg_2_DT);
    mean_t->SetMarkerStyle(22);
    mean_t->SetMarkerColor(kMagenta);
    mean_t->SetMarkerSize(0.9);
    mean_t->Draw("p");

    mean_leg->AddEntry(mean_t , "Channel 15", "p");    
  }

  if(run > 60) {
    TGraph* mean_fo = new TGraph(88, TID, avg_0_DT);
    mean_fo->SetMarkerStyle(4);
    mean_fo->SetMarkerColor(kBlue);
    mean_fo->SetMarkerSize(0.7);
    mean_fo->Draw("p");
    
    mean_leg->AddEntry(mean_fo, "Channel 5" , "p");
  }

  if(run == 62) {
    TGraph* mean_fi = new TGraph(88, TID, avg_4_DT);
    mean_fi->SetMarkerStyle(30);
    mean_fi->Draw("p");

    mean_leg->AddEntry(mean_fi, "Channel 25", "p");
  }
  
  mean_leg->Draw();

  c1->cd(2);
  TGraph* StdDev_f = new TGraph(88, TID, sigma_3_DT);
  StdDev_f->SetMarkerStyle(28);
  StdDev_f->SetTitle("i-th TIGER #sigma_{#Delta T} per Channel");
  StdDev_f->GetXaxis()->SetTitle("TIGER ID");
  StdDev_f->GetXaxis()->SetRangeUser(0, 91);
  StdDev_f->GetYaxis()->SetTitle("#sigma_{#Delta T} [Clock Count]");
  StdDev_f->Draw("ap");

  TLegend* StdDev_leg = new TLegend(0.1, 0.7, 0.48, 0.9 );
  StdDev_leg->AddEntry(StdDev_f , "Channel 20", "p");
  
  if(run > 58) {
    TGraph* StdDev_s = new TGraph(88, TID, sigma_1_DT);
    StdDev_s->SetMarkerStyle(3);
    StdDev_s->SetMarkerColor(kRed);
    StdDev_s->Draw("p");

    StdDev_leg->AddEntry(StdDev_s , "Channel 10", "p");
  }

  if(run > 59) {
    TGraph* StdDev_t = new TGraph(88, TID, sigma_2_DT);
    StdDev_t->SetMarkerStyle(22);
    StdDev_t->SetMarkerColor(kMagenta);
    StdDev_t->SetMarkerSize(0.9);
    StdDev_t->Draw("p");

    StdDev_leg->AddEntry(StdDev_t , "Channel 15", "p");    
  }

  if(run > 60) {
  TGraph* StdDev_fo = new TGraph(88, TID, sigma_0_DT);
  StdDev_fo->SetMarkerStyle(4);
  StdDev_fo->SetMarkerColor(kBlue);
  StdDev_fo->SetMarkerSize(0.7);
  StdDev_fo->Draw("p");

  StdDev_leg->AddEntry(StdDev_fo, "Channel 5" , "p");  
  }

  if(run == 62) {
    TGraph* StdDev_fi = new TGraph(88, TID, sigma_4_DT);
    StdDev_fi->SetMarkerStyle(30);
    StdDev_fi->Draw("p");

    StdDev_leg->AddEntry(StdDev_fi, "Channel 25", "p");
  }

  StdDev_leg->Draw();

}


