#include "daq.h"

//ON and OFF
bool extraction          = 0;
bool general_check       = 0;
bool check_thr           = general_check * true;
bool check_noise         = general_check * true;
bool check_communication = general_check * true;
//Boolean detector
bool check_L1 = true;
bool check_L2 = true;
bool check_planar = false;
bool no_TP = true;

//PARAMETERS FOR THE EVENT/HITS SELECTION
const int INIT_FEB =  0; //0;
const int MAX_FEB  = 55; //55;  
const int N_FEB = MAX_FEB-INIT_FEB+1;
const int N_TIGER = 2*N_FEB;
const int N_CHIP = 2;
const int signal_time_cut_min = -32;
const int signal_time_cut_max =  16;
const int noise_time_cut_min_left  = -180;
const int noise_time_cut_max_left  =  -70;
const int noise_time_cut_min_right =   40;
const int noise_time_cut_max_right =   80;
const int signal_time_cut_min_noTP      = -1420;
const int signal_time_cut_max_noTP      = -1370;
const int noise_time_cut_min_left_noTP  = -1570;
const int noise_time_cut_max_left_noTP  = -1460;
const int noise_time_cut_min_right_noTP = -1240;
const int noise_time_cut_max_right_noTP = -1200;
const int charge_cut = 1;
const int TIME_BIN   = 6.25; //ns
const int NCHANNEL   = 1;//64
const double time_window_noise=(noise_time_cut_max_left-noise_time_cut_min_left + noise_time_cut_max_right-noise_time_cut_min_right)*6.25e-9; //seconds                                                                   
const double time_window_noise_noTP=(noise_time_cut_max_left_noTP-noise_time_cut_min_left_noTP + noise_time_cut_max_right_noTP-noise_time_cut_min_right_noTP)*6.25e-9; //seconds                                                                   
const int thr_ch = 62;
                                                          
//PARAMETERS FOR THE SELECTION OF BAD CHIPS
const double max_rate = 10000;
const double max_noise_charge = 8;
const double max_threshold = 3.;
const double min_average_efficiency = 0.90;
const double min_efficiency = 0.80;
const double failure_rate_efficienct = 0.1;
const double min_saturation_value = 40;

//VARIABLE OF THE QAQC
TString inFile, inFileTP, outFile;
TChain ch("tree");
TChain chTP("t1");
TCanvas *c1 = new TCanvas("c1","c1",800,600);
int error_counter[N_FEB][N_CHIP]; //0->feb_id 1->chip_id
bool error_counter_type[N_FEB][N_CHIP][3]; //0->feb_id 1->chip_i 2->noise/thr/eff
int err_noise(0), err_charge(1), err_comunication(2);
std::fstream outStream;
double expected_noise_rate = 2000; //Hz
double expected_noise_per_chip = time_window_noise*expected_noise_rate*NCHANNEL;

//VARIABLE FOR EXTRACTION
TString ofile_name;
TTree *otree;
TFile *ofile;
int t_feb, t_chip, t_ch, t_strip_x, t_strip_v;
double t_noise, t_thr, t_thr_wid;

//VARIABLE FOR MAPPING
TString map_file="mapping_IHEP_L2_2planari_penta.root";
TFile *mapfile = new TFile(map_file);
TTree *maptree = (TTree*)mapfile->Get("tree");
int channel_id, pos_x, pos_v, chip_id, FEB_label_id;
int mx[N_FEB][N_CHIP][NCHANNEL];
int mv[N_FEB][N_CHIP][NCHANNEL];
int mchip_id[N_FEB][N_CHIP][NCHANNEL];
int mFEB_label_id[N_FEB][N_CHIP][NCHANNEL];

inline double fitfunctionFD(double *x, double *par){ return (par[0]+ par[1]/(1+TMath::Exp(-(x[0]-par[2])/par[3])));}

void daq(int run){
  cout<<"ciao"<<endl;
  daq(run,-1);
}

void daq(int run, int subrun){
  cout<<"++ QAQC analysis and plot"<<endl;
  if(init(run,subrun)){
    cout<<"++ Initialization completed"<<endl;
    if(check_L1)     daq_i( 0,15); //L1
    if(check_L2)     daq_i(16,43); //L2
    if(check_planar) daq_i(44,55); //planari
  }
  outStream.close();
  cout<<"++ Terminated QAQC"<<endl;
  return;
}

void daq_i(int start_FEB, int end_FEB){
  for(int ifeb=0;ifeb<=N_FEB;ifeb++) for(int ichip=0;ichip<N_CHIP;ichip++) {error_counter[ifeb][ichip]=0;for(int i_type=0;i_type<3;i_type++)error_counter_type[ifeb][ichip][i_type]=0;}
  for(int FEB_i=start_FEB;FEB_i<=end_FEB;FEB_i++){
    if(FEB_i==44 || FEB_i==45 || FEB_i==46 || FEB_i==47 || FEB_i==48 || FEB_i==50 || FEB_i==51) continue;
    for(int chip_i=1; chip_i<3; chip_i++){
      if(check_noise)         noise_check(FEB_i,chip_i);
      if(check_thr)           charge_check(FEB_i,chip_i);
      if(check_communication) comunication_efficiency(FEB_i,chip_i);
      if(extraction)          extract_noise_and_thr(FEB_i,chip_i);
      print_error(FEB_i,chip_i);
    }
  }
  count_error(start_FEB, end_FEB);  
  if(extraction) ofile->Write();
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
    //inFile = DataDir+to_string(run)+Form("/Sub_RUN_event_%d.root",subrun);
    //inFileTP = DataDir+to_string(run)+Form("/Sub_RUN_TP_event_%d.root",subrun);
    inFile = DataDir+to_string(run)+Form("/event_till_%d.root",subrun); 
    inFileTP = DataDir+to_string(run)+Form("/TP_event_till_%d.root",subrun);
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
  for(int ifeb=0;ifeb<=N_FEB;ifeb++) for(int ichip=0;ichip<N_CHIP;ichip++) {error_counter[ifeb][ichip]=0;for(int i_type=0;i_type<3;i_type++)error_counter_type[ifeb][ichip][i_type]=0;}
  //Define the extraction tree
  if(extraction){
    ofile_name=Form("/home/ihep_data/data/raw_daq/extracted_noise_thr_%d.root",run);
    ofile = new TFile(ofile_name,"RECREATE");
    otree = new TTree("tree","tree");
    otree->Branch("FEB",&t_feb,"FEB/I");
    otree->Branch("chip",&t_chip,"chip/I");
    otree->Branch("channel",&t_ch,"channel/I");
    otree->Branch("noise",&t_noise,"noise/D");
    otree->Branch("threshold",&t_thr,"threshold/D");
    otree->Branch("thr_width",&t_thr_wid,"thr_width/D");
    otree->Branch("strip_x",&t_strip_x,"strip_x/I");
    otree->Branch("strip_v",&t_strip_v,"strip_v/I");
  }
  //Define the mapping
  maptree->SetBranchAddress("channel_id", &channel_id);
  maptree->SetBranchAddress("pos_x", &pos_x);
  maptree->SetBranchAddress("pos_v", &pos_v);
  maptree->SetBranchAddress("chip_id", &chip_id);
  maptree->SetBranchAddress("FEB_label", &FEB_label_id);
  memset(mx, -1, sizeof(mx));
  memset(mv, -1, sizeof(mv));
  for (int i = 0; i < maptree->GetEntries(); i++) {
    maptree->GetEntry(i);
    if(pos_x*pos_v<0 && (mx[FEB_label_id][chip_id-1][channel_id])==-1 && mv[FEB_label_id][chip_id-1][channel_id]==-1 ){ 
      mx[FEB_label_id][chip_id-1][channel_id] = pos_x;
      mv[FEB_label_id][chip_id-1][channel_id] = pos_v;
    }
  }
  return true;
}

TString Get_Cut(int caso, int var){
  return Get_Cut(caso,var,-1);
}

TString Get_Cut(int caso, int var1, int var2){
  return Get_Cut(caso,var1,var2,-1);
}

TString Get_Cut(int caso, int var1, int var2, int var3){
  TString cut = "ciao";
  switch(caso){
  case 0:
    //Noise check
    cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&((t_min_ttrigg>%d && t_min_ttrigg<%d)|| (t_min_ttrigg>%d&&t_min_ttrigg<%d))&&charge_SH>%d",thr_ch,var1,var2,noise_time_cut_min_left,noise_time_cut_max_left,noise_time_cut_min_right,noise_time_cut_max_right,charge_cut);
    return cut;
  case 1:
    //Threshold check
    if(!no_TP) cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&t_min_ttrigg>%d && t_min_ttrigg<%d&&charge_SH>%d",thr_ch,var1,var2,signal_time_cut_min,signal_time_cut_max,charge_cut);
    else cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&tcoarse_min_ts>%d && tcoarse_min_ts<%d&&charge_SH>%d",thr_ch,var1,var2,signal_time_cut_min_noTP,signal_time_cut_max_noTP,charge_cut);
    return cut;
  case 2:
    //TP comunication check
    //cut = Form("eff_Tpm1[%d] < %f",var1,min_efficiency);
    cut = Form("FEB==%d && chip==%d", var1, var2);
    return cut;
  case 3:
    //Noise check per channel
    if(!no_TP) cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&((t_min_ttrigg>%d && t_min_ttrigg<%d)|| (t_min_ttrigg>%d&&t_min_ttrigg<%d))&&charge_SH>%d && channel==%d",thr_ch,var1,var2,noise_time_cut_min_left,noise_time_cut_max_left,noise_time_cut_min_right,noise_time_cut_max_right,charge_cut,var3);
    else cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&((tcoarse_min_ts>%d && tcoarse_min_ts<%d)|| (tcoarse_min_ts>%d&&tcoarse_min_ts<%d))&&charge_SH>%d && channel==%d",thr_ch,var1,var2,noise_time_cut_min_left_noTP,noise_time_cut_max_left_noTP,noise_time_cut_min_right_noTP,noise_time_cut_max_right_noTP,charge_cut,var3);
    return cut;
  case 4:
    //Noise check per channel 
    cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&charge_SH>%d && channel==%d",thr_ch,var1,var2,charge_cut,var3);
    cout<<cut<<endl;
    return cut;
  case 5:
    cut = Form("TP_eff<%f && FEB==%d && chip==%d", min_efficiency, var1, var2);
    return cut;
  default:
    cout<<"No cut has been selected"<<endl;
    return cut;
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
    //command = Form("eff_Tpm1[%d]",var) + (string)">>" + h_name;
    command = "TP_eff>>"+h_name;
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
    error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_noise]++;
  }
  else{
    int hit_of_noise = h_noise->GetEntries();
    double charge_noise = h_noise->GetMean();
    double rate = (double)hit_of_noise/time_window_noise/n_events/NCHANNEL;
    if(print_here) cout<<FEB_i<<" "<<chip_i<<" "<<hit_of_noise<<" "<<charge_noise<<" "<<rate<<endl;
    if(rate>=max_rate) {
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a RATE of "<<rate<<" Hz --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
      error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_noise]++;
    }
    if(charge_noise>=max_noise_charge) {
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a MEAN CHARGE of "<<charge_noise<<" fC --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
      error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_noise]++;
    }
  }
  delete h_noise;
  return;
}

void charge_check(int FEB_i, int chip_i){
  bool print_here = false;
  int min_charge_hist=0;
  int max_charge_hist=70;
  double charge_resolution = 0.25; //fC
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
    error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_charge]++;
  }
  else{
    int max_bin_thr = -1;
    int min_bin_thr = -1;
    if(print_here) for(int ibin=1;ibin<=n_bin;ibin++) cout<<ibin<<" "<<h_charge->GetBinContent(ibin)<<endl;
    for(int ibin=1;ibin<=n_bin;ibin++) if(h_charge->GetBinContent(ibin)>h_charge->GetBinContent(ibin-1) && h_charge->GetBinContent(ibin)>h_charge->GetBinContent(ibin+1) && h_charge->GetBinContent(ibin)>h_charge->GetBinContent(ibin+5) && h_charge->GetBinContent(ibin)>h_charge->GetBinContent(ibin+6) && max_bin_thr==-1){max_bin_thr=ibin; break;}
    for(int ibin=1;ibin<=max_bin_thr;ibin++) if(h_charge->GetBinContent(ibin)>(double)0.1*h_charge->GetMaximum() && min_bin_thr==-1) {min_bin_thr=ibin; break;}
    double thresh = max_bin_thr*charge_resolution;
    double mean_thr = (max_bin_thr+min_bin_thr)*charge_resolution/2.;
    double sigma_thr = (max_bin_thr-min_bin_thr)*charge_resolution/2.;
    if(sigma_thr<=charge_resolution) sigma_thr=2*charge_resolution;
    TF1 *f_thr = new TF1("f_thr",fitfunctionFD,min_bin_thr-5,max_bin_thr+3,4);
    f_thr->SetParameters(0,h_charge->GetBinContent(max_bin_thr),mean_thr,sigma_thr);
    f_thr->SetParLimits(0,-0.1*h_charge->GetBinContent(max_bin_thr),0.1*h_charge->GetBinContent(max_bin_thr));
    f_thr->SetParLimits(1,h_charge->GetBinContent(max_bin_thr)*0.5,h_charge->GetBinContent(max_bin_thr)*2);
    f_thr->SetParLimits(2,mean_thr-sigma_thr,mean_thr+sigma_thr);
    f_thr->SetParLimits(3,0.1*sigma_thr,1.2*sigma_thr);
    h_charge->Fit("f_thr","Q","",0,max_bin_thr*charge_resolution+2);
    double thr_fit = f_thr->GetParameter(2);
    double err_thr_fit = f_thr->GetParError(2);
    if(abs(thr_fit-mean_thr)<2*sigma_thr) {mean_thr=thr_fit; sigma_thr=err_thr_fit;}
    delete f_thr;
    if(print_here) cout<<"Bins for thr: "<<min_bin_thr<<" "<<max_bin_thr<<" -< THR: "<<mean_thr<<" +- "<<sigma_thr<<endl;
    if(mean_thr >= max_threshold){
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a THRESHOLD of "<<mean_thr<<" +- "<<sigma_thr <<"fC --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
      error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_charge]++;
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
      error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_charge]++;
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
  //if(print_here) cout<<Get_Cut(2,chipID)<<endl; 
  if(print_here) cout<<Get_Cut(2,FEB_i,chip_i)<<endl; 
  chTP.Draw(Get_Command(1,h_name,chipID),Get_Cut(2,FEB_i,chip_i));
  if(h_TP_eff->GetEntries()==0){
    cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has NO TP communication informations --> CONTROLL IT"<<endl;
    error_counter[FEB_i-INIT_FEB][chip_i-1]++;
    error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_comunication]++;
   }
  else{
    double mean_eff_TP = h_TP_eff->GetMean();
    double n_eff_below_min = chTP.GetEntries(Get_Cut(5,chipID));
    double ratio_eff_below = n_eff_below_min/ h_TP_eff->GetEntries();
    if(print_here) cout<<chipID<<" "<<mean_eff_TP<<" "<<n_eff_below_min<<" "<<ratio_eff_below<<endl;
    if(mean_eff_TP<min_average_efficiency){
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a mean TP comunication EFFICIENCY of "<<mean_eff_TP*100<<"\% --> TOO LOW"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
      error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_comunication]++;
    }
    if(ratio_eff_below>failure_rate_efficienct){
      cout<<"The FEB "<<FEB_i<<" chip "<<chip_i<<" has a TP comunication failure rate (eff<"<<min_efficiency*100<<"\%) of "<<ratio_eff_below*100<<"\% --> TOO HIGH"<<endl;
      error_counter[FEB_i-INIT_FEB][chip_i-1]++;
      error_counter_type[FEB_i-INIT_FEB][chip_i-1][err_comunication]++;
     }
  }
  delete h_TP_eff;
  return;
}

void extract_noise_and_thr(int FEB_i, int chip_i){
  bool print_here = false;
  for(int ch_i=0;ch_i<NCHANNEL;ch_i++){
    //noise
    double nhitsofnoise = ch.GetEntries(Get_Cut(3,FEB_i,chip_i,ch_i));
    double rate = 0;
    if(!no_TP) rate = nhitsofnoise/time_window_noise/ch.GetEntries();
    else       rate = nhitsofnoise/time_window_noise_noTP/ch.GetEntries();
    //thr
    if(print_here) cout<<"rate: "<<rate<<endl;
    int min_charge_hist=0;
    int max_charge_hist=30;
    double charge_resolution = 0.25; //fC                                                                                                                                                                                 
    int n_bin = (max_charge_hist-min_charge_hist)/charge_resolution;
    TString h_name="h_charge";
    TH1D *h_charge = new TH1D(h_name,h_name,n_bin,min_charge_hist,max_charge_hist);
    TF1 *f_thr = new TF1("f_thr","[0]/(1+exp(-(x-[1])/[2]))");
    if(print_here) cout<<Get_Command(0,h_name)<<endl;
    if(print_here) cout<<Get_Cut(4,FEB_i,chip_i,ch_i)<<endl;
    ch.Draw(Get_Command(0,h_name),Get_Cut(4,FEB_i,chip_i,ch_i));
    double threshold = 0;
    double thr_wid = 0;
    if(h_charge->GetEntries()){
      f_thr->SetParameters(h_charge->GetMaximum(),0.5*h_charge->GetMaximumBin()*charge_resolution,0.2);
      f_thr->SetParLimits(0,0,1.2*h_charge->GetMaximum());
      f_thr->SetParLimits(1,0,h_charge->GetMaximumBin()*charge_resolution+1);
      f_thr->SetParLimits(2,0,1);
      if(print_here){
	cout<<0<<" "<<0<<" "<<1.5*h_charge->GetMaximum()<<endl;
	cout<<1<<" "<<0<<" "<<h_charge->GetMaximumBin()*charge_resolution<<endl;
	cout<<2<<" "<<0<<" "<<2<<endl;
	cout<<h_charge->GetMaximumBin()<<" "<<h_charge->GetMaximumBin()*charge_resolution<<endl;
	for(int i=0;i<n_bin;i++)cout<<h_charge->GetBinContent(i+1)<<" ";
	cout<<endl;
      }
      h_charge->Fit("f_thr","Q","",h_charge->GetMaximumBin()*charge_resolution-2,h_charge->GetMaximumBin()*charge_resolution+1);
      threshold=f_thr->GetParameter(1);
      thr_wid  =f_thr->GetParameter(2);
    }
    cout<<FEB_i<<" "<<chip_i<<" "<<ch_i<<" "<<rate<<" "<<threshold<<" "<<mx[FEB_i][chip_i-1][ch_i]<<" "<<mv[FEB_i][chip_i-1][ch_i]<<endl;
    cout<<"*****************************"<<endl;
    //Fill the ntuple
    t_feb=FEB_i;
    t_chip=chip_i;
    t_ch=ch_i;
    t_noise=rate;
    t_thr=threshold;
    t_thr_wid=thr_wid;
    t_strip_x=mx[FEB_i][chip_i-1][ch_i];
    t_strip_v=mv[FEB_i][chip_i-1][ch_i];
    otree->Fill();
    delete h_charge, f_thr;
  }
  return;
}

void print_error(int FEB_i, int chip_i){
  TString out = Form("%d %d %d",FEB_i,chip_i,error_counter[FEB_i-INIT_FEB][chip_i-1]);
  outStream<<out<<endl;
  return;
}

void count_error(int f_start, int f_end){
  int nerr=0;
  int nerr_type[3]={0,0,0};
  for(int FEB_i=f_start;FEB_i<=f_end;FEB_i++){
    for(int chip_i=1; chip_i<3; chip_i++){
      nerr+=error_counter[FEB_i-INIT_FEB][chip_i-1];
      for(int type_i=0; type_i<3; type_i++){
	nerr_type[type_i]+=error_counter_type[FEB_i-INIT_FEB][chip_i-1][type_i];
      }
    }
  }  
  cout<<"***************************"<<endl;
  cout<<"Number of total errors: "<<nerr<<endl;
  for(int type_i=0;type_i<3; type_i++){
    cout<<"Number of chip with error type "<<error_type(type_i)<<": "<<nerr_type[type_i]<<endl;
  }
  cout<<"***************************"<<endl;
  return;
}

string error_type(int type){
  if(type==0) return "noise";
  if(type==1) return "threshold";
  if(type==2) return "communication";
  return "dunno";
}
