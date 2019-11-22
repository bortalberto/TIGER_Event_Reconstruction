#include "daq.h"

//PARAMETERS FOR THE EVENT/HITS SELECTION
const int INIT_FEB = 16; //16;
const int MAX_FEB  = 43; //43;  
const int N_FEB = MAX_FEB-INIT_FEB+1;
const int N_TIGER = 2*N_FEB;
const int N_CHIP = 2;
const int signal_time_cut_min = -32;
const int signal_time_cut_max =  16;
const int noise_time_cut_min_left  = -80;
const int noise_time_cut_max_left  = -40;
const int noise_time_cut_min_right =  40;
const int noise_time_cut_max_right =  80;
const int charge_cut = 0;
const int TIME_BIN   = 6.25; //ns
const int NCHANNEL   = 64;
const double time_window_noise=(noise_time_cut_max_left-noise_time_cut_min_left + noise_time_cut_max_right-noise_time_cut_min_right)*6.25e-9; //seconds                                                                                                                             



//PARAMETERS FOR THE SELECTION OF BAD CHIPS
const double max_rate = 10000;
const double max_noise_charge = 10;
const double max_threshold = 5.0;
const double min_average_efficiency = 0.90;
const double min_efficiency = 0.8;
const double failure_rate_efficienct = 0.1;
const double min_saturation_value = 40;

//VARIABLE OF THE QAQC
TString inFile, inFileTP, outFile;
TChain ch("tree");
TChain chTP("t1");
TCanvas *c1 = new TCanvas("c1","c1",800,600);
int error_counter[N_FEB][N_CHIP]; //0->feb_id 1->chip_id
std::fstream outStream;
double expected_noise_rate = 2000; //Hz
double expected_noise_per_chip = time_window_noise*expected_noise_rate*NCHANNEL;

void daq(int run){
  cout<<"ciao"<<endl;
  daq(run,-1);
}

void daq(int run, int subrun){
  cout<<"++ QAQC analysis and plot"<<endl;
  if(init(run,subrun)){
    cout<<"++ Initialization completed"<<endl;
    for(int FEB_i=INIT_FEB;FEB_i<=MAX_FEB;FEB_i++){
      for(int chip_i=1; chip_i<3; chip_i++){
	noise_check(FEB_i,chip_i);
	charge_check(FEB_i,chip_i);
	comunication_efficiency(FEB_i,chip_i);
	print_error(FEB_i,chip_i);
      }
    }
  }
  outStream.close();
  cout<<"++ Terminated QAQC"<<endl;
  return;
}

bool init(int run, int subrun){
  //open the root file
  if(subrun<0){
    inFile = DataDir+to_string(run)+"/event.root";
    inFileTP = DataDir+to_string(run)+"/TP_event.root";
    outFile = DataDir+to_string(run)+"/QAQC_log.txt";
  }
  else{
    inFile = DataDir+to_string(run)+Form("/Sub_RUN_event_%d.root",subrun);
    inFileTP = DataDir+to_string(run)+Form("/Sub_RUN_TP_event_%d.root",subrun);
    outFile = DataDir+to_string(run)+Form("/QAQC_log_%d.txt",subrun);
  }
  cout<<"*******************************************************************"<<endl;
  cout<<"**  The file analyzed is: "<<inFile<<" and "<<inFileTP<<endl;
  cout<<"**  The FEB range goes from "<<INIT_FEB<<" to "<<MAX_FEB<<endl;
  cout<<"**  Noise is evaluated from "<<noise_time_cut_min_left*TIME_BIN<<" to "<<noise_time_cut_max_left*TIME_BIN<<" ns and from "<<noise_time_cut_min_right*TIME_BIN<<" to "<<noise_time_cut_max_right*TIME_BIN<<" ns"<<endl;
  cout<<"**  Signal is evalueted from "<<signal_time_cut_min*TIME_BIN<<" to "<<signal_time_cut_max*TIME_BIN<<" ns"<<endl;
  cout<<"**  Charge cut used is Q(fC)>"<<charge_cut<<endl;
  cout<<"!!  The maximum mean noise rate per channel for good chips is "<<max_rate<<" Hz"<<endl;
  cout<<"!!  The maximum mean charge noise for good chips is "<<max_noise_charge<<" fC"<<endl;
  cout<<"!!  The maximum threshold for good chips is "<<max_threshold<<" fC"<<endl;
  cout<<"!!  The limit for the mean efficiency for good chips is "<<min_average_efficiency*100<<"\%"<<endl;
  cout<<"!!  The maximum failure (eff<"<<min_efficiency*100<<"%) rate for good chips is "<<failure_rate_efficienct*100<<"\%"<<endl;
  cout<<"!!  The minimum saturation for good chips is "<<min_saturation_value<<" fC"<<endl;
  cout<<"*******************************************************************"<<endl;
  std::ifstream inStream(inFile, std::ios::binary);
  if (!inStream) {
    std::cerr << "File " << inFile << " not found" << std::endl;
    return false;
  }
  std::ifstream inStreamTP(inFileTP, std::ios::binary);
  if (!inStreamTP) {
    std::cerr << "File " << inFileTP << " not found" << std::endl;
    return false;
  }
  outStream.open(outFile,std::ios_base::out);
  if (!outStream) {
    std::cerr << "File " << outFile << " not found" << std::endl;
    return false;
  }
  ch.Add(inFile);
  chTP.Add(inFileTP);
  expected_noise_per_chip*=ch.GetEntries();
  outStream<<"FEB CHIP N_ERROR"<<endl;
  for(int ifeb=0;ifeb<=N_FEB;ifeb++) for(int ichip=0;ichip<N_CHIP;ichip++) error_counter[ifeb][ichip]=0;
  return true;
}

TString Get_Cut(int caso, int var){
  return Get_Cut(caso,var,-1);
}

TString Get_Cut(int caso, int var1, int var2){
  TString command = "ciao";
  switch(caso){
  case 0:
    //Noise check
    command = Form("channel!=62&&FEB_label==%d&&chip==%d&&((t_min_ttrigg>%d && t_min_ttrigg<%d)|| (t_min_ttrigg>%d&&t_min_ttrigg<%d))&&charge_SH>%d",var1,var2,noise_time_cut_min_left,noise_time_cut_max_left,noise_time_cut_min_right,noise_time_cut_max_right,charge_cut);
    return command;
  case 1:
    //Threshold check
    command = Form("channel!=62&&FEB_label==%d&&chip==%d&&t_min_ttrigg>%d && t_min_ttrigg<%d&&charge_SH>%d",var1,var2,signal_time_cut_min,signal_time_cut_max,charge_cut);
    return command;
  case 2:
    //TP comunication check
    command = Form("eff_Tpm1[%d] < %f",var1,min_efficiency);
    return command;
  default:
    cout<<"No cut has been selected"<<endl;
    return command;
  }
}

TString Get_Command(int caso, TString h_name){
  return Get_Command(caso, h_name, -1);  
}

TString Get_Command(int caso, TString h_name, int var){
  TString command;
  switch(caso){
  case 0:
    //Charge distribution
    command = "charge_SH>>"+h_name;
    return command;
  case 1:
    //TP comunication efficiency
    command = Form("eff_Tpm1[%d]",var) + (string)">>" + h_name;
    return command;
  default:
    cout<<"No command has been selected"<<endl;
    return command;
  }
}

void noise_check(int FEB_i, int chip_i){
  bool print_here = false;
  int n_events = ch.GetEntries();
  TString h_name="h_noise";
  TH1D *h_noise = new TH1D(h_name,h_name,100,0,60);
  if(print_here) cout<<ch.GetEntries()<<endl;
  if(print_here) cout<<Get_Cut(0,FEB_i,chip_i)<<endl;
  if(print_here) cout<<ch.GetEntries(Get_Cut(0,FEB_i,chip_i))<<endl;
  if(print_here) cout<<Get_Command(0,h_name)<<endl;
  ch.Draw(Get_Command(0,h_name),Get_Cut(0,FEB_i,chip_i));
  if(h_noise->GetEntries()==0 && NCHANNEL<expected_noise_per_chip){
    cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has NO NOISE HITS --> CONTROLL IT"<<endl;
    error_counter[FEB_i-INIT_FEB][chip_i-1]++;
  }
  else{
    int hit_of_noise = h_noise->GetEntries();
    double charge_noise = h_noise->GetMean();
    double rate = (double)hit_of_noise/time_window_noise/n_events/NCHANNEL;
    if(print_here) cout<<FEB_i<<" "<<chip_i<<" "<<hit_of_noise<<" "<<charge_noise<<" "<<rate<<endl;
    if(rate>=max_rate) {
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a RATE of "<<rate<<" Hz --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
    }
    if(charge_noise>=max_noise_charge) {
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a MEAN CHARGE of "<<charge_noise<<" fC --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
    }
  }
  delete h_noise;
  return;
}

void charge_check(int FEB_i, int chip_i){
  bool print_here = false;
  int min_charge_hist=0;
  int max_charge_hist=70;
  double charge_resolution = 0.5; //fC
  int n_bin = (max_charge_hist-min_charge_hist)/charge_resolution;
  TString h_name="h_charge";
  TH1D *h_charge = new TH1D(h_name,h_name,n_bin,min_charge_hist,max_charge_hist);
  if(print_here) cout<<ch.GetEntries()<<endl; 
  if(print_here) cout<<Get_Cut(1,FEB_i,chip_i)<<endl; 
  if(print_here) cout<<ch.GetEntries(Get_Cut(1,FEB_i,chip_i))<<endl; 
  if(print_here) cout<<Get_Command(0,h_name)<<endl;
  ch.Draw(Get_Command(0,h_name),Get_Cut(1,FEB_i,chip_i));
  if(h_charge->GetEntries()==0 && NCHANNEL<expected_noise_per_chip){
    cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has NO HITS --> CONTROLL IT"<<endl;
    error_counter[FEB_i-INIT_FEB][chip_i-1]++;
  }
  else{
    int max_bin_thr = -1;
    int min_bin_thr = -1;
    if(print_here) for(int ibin=1;ibin<=n_bin;ibin++) cout<<ibin<<" "<<h_charge->GetBinContent(ibin)<<endl;
    for(int ibin=1;ibin<=n_bin;ibin++) if(h_charge->GetBinContent(ibin)>h_charge->GetBinContent(ibin-1) && h_charge->GetBinContent(ibin)>h_charge->GetBinContent(ibin+1) && max_bin_thr==-1){max_bin_thr=ibin; break;}
    for(int ibin=1;ibin<=max_bin_thr;ibin++) if(h_charge->GetBinContent(ibin)>(double)0.5*h_charge->GetMaximum() && min_bin_thr==-1) {min_bin_thr=ibin; break;}
    double thresh = max_bin_thr*charge_resolution;
    double mean_thr = (max_bin_thr+min_bin_thr)*charge_resolution/2.;
    double sigma_thr = (max_bin_thr-min_bin_thr)*charge_resolution/2.;
    if(print_here) cout<<"Bins for thr: "<<min_bin_thr<<" "<<max_bin_thr<<" -< THR: "<<mean_thr<<" +- "<<sigma_thr<<endl;
    if(thresh >= max_threshold){
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a THRESHOLD of "<<mean_thr<<" +- "<<sigma_thr <<"fC --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
    }
    double min_sat_bin = -1;
    double max_sat_bin = -1;
    for(int ibin=n_bin-1;ibin>0;ibin--) if(h_charge->GetBinContent(ibin)>0 && h_charge->GetBinContent(ibin-1)>0) {max_sat_bin=ibin; break;}
    for(int ibin=max_sat_bin;ibin>0;ibin--) if(h_charge->GetBinContent(ibin)<1.2*h_charge->GetBinContent(ibin+1) && h_charge->GetBinContent(ibin)<1.2*h_charge->GetBinContent(ibin-1)) {min_sat_bin=ibin; break;}
    if(print_here) cout<<min_sat_bin<<" "<<max_sat_bin<<" "<<min_sat_bin*charge_resolution<<" "<<max_sat_bin*charge_resolution<<" "<<0.5*(max_sat_bin+min_sat_bin)*charge_resolution<<" +- "<<0.5*(max_sat_bin-min_sat_bin)*charge_resolution<<endl;
    double mean_sat = 0.5*(max_sat_bin+min_sat_bin)*charge_resolution;
    double sigma_sat = 0.5*(max_sat_bin-min_sat_bin)*charge_resolution;
    if(mean_sat<min_saturation_value){
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a SATURATION level of "<<mean_sat<<" +- "<<sigma_sat<<" --> TOO LOW"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
    }
  }
  delete h_charge;
  return;
}

void comunication_efficiency(int FEB_i, int chip_i){
  bool print_here = false;
  TString h_name="h_TP_eff";
  TH1D *h_TP_eff = new TH1D(h_name,h_name,1000,0,1);
  int chipID = FEB_i*2+chip_i-1;
  if(print_here) cout<<chTP.GetEntries()<<endl;
  if(print_here) cout<<Get_Command(1, h_name, chipID)<<endl;
  if(print_here) cout<<Get_Cut(2,chipID)<<endl; 
  chTP.Draw(Get_Command(1,h_name,chipID));
  if(h_TP_eff->GetEntries()==0){
    cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has NO TP communication informations --> CONTROLL IT"<<endl;
    error_counter[FEB_i-INIT_FEB][chip_i-1]++;
  }
  else{
    double mean_eff_TP = h_TP_eff->GetMean();
    double n_eff_below_min = chTP.GetEntries(Get_Cut(2,chipID));
    double ratio_eff_below = n_eff_below_min/chTP.GetEntries();
    if(print_here) cout<<chipID<<" "<<mean_eff_TP<<" "<<n_eff_below_min<<" "<<ratio_eff_below<<endl;
    if(mean_eff_TP<min_average_efficiency){
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a mean TP comunication EFFICIENCY of "<<mean_eff_TP*100<<"\% --> TOO LOW"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
    }
    if(ratio_eff_below>failure_rate_efficienct){
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a TP comunication failure rate (eff<"<<min_efficiency*100<<"\%) of "<<ratio_eff_below*100<<"\% --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
    }
  }
  delete h_TP_eff;
  return;
}

void print_error(int FEB_i, int chip_i){
  TString out = Form("%d %d %d",FEB_i,chip_i,error_counter[FEB_i-INIT_FEB][chip_i-1]);
  outStream<<out<<endl;
  return;
}
