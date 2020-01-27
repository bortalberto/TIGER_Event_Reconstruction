#include "ext.h"

//ON and OFF
bool extraction          = 1;
bool general_check       = 0;
bool check_thr           = general_check * true;
bool check_noise         = general_check * true;
bool check_communication = general_check * true;
//Boolean detector
bool check_L1 = true;
bool check_L2 = true;
bool check_planar = false;
bool no_TP = false;

//PARAMETERS FOR THE EVENT/HITS SELECTION
const int INIT_FEB =  0; //0;
const int MAX_FEB  = 55; //55;  
const int N_FEB = MAX_FEB-INIT_FEB+1;
const int N_TIGER = 2*N_FEB;
const int N_CHIP = 2;
const int signal_time_cut_min = -32;
const int signal_time_cut_max =  16;
const int noise_time_cut_min_left  = -120;
const int noise_time_cut_max_left  =  -70;
const int noise_time_cut_min_right =   80;//40;
const int noise_time_cut_max_right =   80;//80;
const int signal_time_cut_min_noTP      = -1420;
const int signal_time_cut_max_noTP      = -1380;
const int noise_time_cut_min_left_noTP  = -1500;
const int noise_time_cut_max_left_noTP  = -1450;
const int noise_time_cut_min_right_noTP = -1300;//-1240;
const int noise_time_cut_max_right_noTP = -1300;
const int charge_cut = -5;
const int TIME_BIN   = 6.25; //ns
const int NCHANNEL   = 64;
const double time_window_noise=(noise_time_cut_max_left-noise_time_cut_min_left + noise_time_cut_max_right-noise_time_cut_min_right)*6.25e-9; //seconds                                                                   
const double time_window_noise_noTP=(noise_time_cut_max_left_noTP-noise_time_cut_min_left_noTP + noise_time_cut_max_right_noTP-noise_time_cut_min_right_noTP)*6.25e-9; //seconds                                                                   
int thr_ch = 62;
                                                          
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

//VARIABLE FOR COMMANDS
int n_runs=0;

inline double fitfunctionFD(double *x, double *par){ return (par[0]+ par[1]/(1+TMath::Exp(-(x[0]-par[2])/par[3])));}

void ext(int run, int FEB_i, int chip_i, int ch_i){
  bool print_here = false;
  inFile = DataDir+to_string(run)+"/event.root";
  ch.Add(inFile);
  //noise                  
  int n_bin_noise=100;
  TH1D *h_noise;
  TString name_noise = "h_noise";
  if(no_TP){
    n_bin_noise= noise_time_cut_max_right_noTP - noise_time_cut_min_left_noTP;
    h_noise = new TH1D("h_noise","h_noise",n_bin_noise,noise_time_cut_min_left_noTP,noise_time_cut_max_right_noTP);
    //ch.Draw("tcoarse_min_ts>>h_noise",Get_Cut(3,FEB_i,chip_i,ch_i),"goff");
    //cout<<Get_Cut(3,FEB_i,chip_i,ch_i)<<endl;
  }
  if(!no_TP) {
    n_bin_noise= noise_time_cut_max_right - noise_time_cut_min_left;
    h_noise = new TH1D("h_noise","h_noise",n_bin_noise,noise_time_cut_max_right,noise_time_cut_min_left);
    //ch.Draw("t_min_ttrigg>>h_noise",Get_Cut(3,FEB_i,chip_i,ch_i),"goff");
    //cout<<Get_Cut(3,FEB_i,chip_i,ch_i)<<endl;
  }
  if(print_here)cout<<Get_Command(2,name_noise)<<endl;
  if(print_here)cout<<Get_Cut(3,FEB_i,chip_i,ch_i)<<endl;
  ch.Draw(Get_Command(2,name_noise),Get_Cut(3,FEB_i,chip_i,ch_i),"goff");
  //double nhitsofnoise = ch.GetEntries(Get_Cut(3,FEB_i,chip_i,ch_i));
  double nhitsofnoise = h_noise->GetEntries();
  double rate = 0;
  if(!no_TP) rate = nhitsofnoise/time_window_noise/ch.GetEntries();
  else       rate = nhitsofnoise/time_window_noise_noTP/ch.GetEntries();
  //thr                                                                                                                                                                                                                  
  int min_charge_hist=0;
  int max_charge_hist=30;
  double charge_resolution = 0.25; //fC                                                                                                                                                                                 
  int n_bin = (max_charge_hist-min_charge_hist)/charge_resolution;
  TString h_name="h_charge";
  TH1D *h_charge = new TH1D(h_name,h_name,n_bin,min_charge_hist,max_charge_hist);
  TF1 *f_thr = new TF1("f_thr","[0]/(1+exp(-(x-[1])/[2]))");
  ch.Draw(Get_Command(0,h_name),Get_Cut(4,FEB_i,chip_i,ch_i),"goff");
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
  //mapping
  maptree->SetBranchAddress("channel_id", &channel_id);
  maptree->SetBranchAddress("pos_x", &pos_x);
  maptree->SetBranchAddress("pos_v", &pos_v);
  maptree->SetBranchAddress("chip_id", &chip_id);
  maptree->SetBranchAddress("FEB_label", &FEB_label_id);
  int mx(-1),mv(-1);
  for (int i = 0; i < maptree->GetEntries(); i++) {
    maptree->GetEntry(i);
    if(pos_x*pos_v<0 && mx==-1 && mv==-1 && FEB_label_id==FEB_i && chip_id==chip_i && channel_id==ch_i){
      mx = pos_x;
      mv = pos_v;
      break;
    }
  }
  if(print_here)cout<<"FEB: "<<FEB_i<<" Chip: "<<chip_i<<" Channel: "<<ch_i<<" Rate: "<<rate<<" Thr: "<<threshold<<" StripX: "<<mx<<" StripV: "<<mv<<endl;
  outFile = DataDir+to_string(run)+"/extraction.txt";
  outStream.open(outFile,std::ios::out |std::ios::app);
  outStream<<FEB_i<<" "<<chip_i<<" "<<ch_i<<" "<<rate<<" "<<threshold<<" "<<thr_wid<<" "<<mx<<" "<<mv<<endl;
  outStream.close();
  return;
}

int count_line(int run){
  outFile =  DataDir+to_string(run)+"/extraction.txt";
  int numLines = 0;                                                                                                                                                                                                         ifstream in(outFile);                                                                                                                                                                                                     std::string unused;                                                                                                                                                                                                       while ( std::getline(in, unused) ) ++numLines; 
  return numLines;
}

void ext(int run){
  //cout<<"++ QAQC analysis and plot"<<endl;
  if(init(run)){
    cout<<"++ Initialization completed"<<endl;
    if(check_L1)     ext_i(run, 0,15); //L1
    if(check_L2)     ext_i(run,16,43); //L2
    if(check_planar) ext_i(run,44,55); //planari
  }
  //cout<<"++ Terminated QAQC"<<endl;
  cout<<"Number of running channels: "<<n_runs<<endl;
  while(count_line(run)<n_runs) {
    cout<<"Number of analyzed channels: "<<count_line(run)<<endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(10000));
  }
  outFile = DataDir+to_string(run)+"/extraction.txt";
  ifstream in(outFile); 
  while(in>>t_feb>>t_chip>>t_ch>>t_noise>>t_thr>>t_thr_wid>>t_strip_x>>t_strip_v) otree->Fill();
  otree->Write();
  return;
}

void ext_i(int run, int start_FEB, int end_FEB){
  for(int FEB_i=start_FEB;FEB_i<=end_FEB;FEB_i++){
    if(FEB_i==44 || FEB_i==45 || FEB_i==46 || FEB_i==47 || FEB_i==48 || FEB_i==50 || FEB_i==51) continue;
    for(int chip_i=1; chip_i<3; chip_i++){
      for(int ch_i=0;ch_i<NCHANNEL;ch_i++){
	TString bash_command = Form("ts $exe_ter -X %d %d %d %d",run,FEB_i,chip_i,ch_i);
	//cout<<bash_command<<endl;
	gSystem->Exec(bash_command);
	n_runs++;
      }
    }
  }
  return;
}

bool init(int run){
  //open the root file
  inFile = DataDir+to_string(run)+"/event.root";
  outFile = DataDir+to_string(run)+"/QAQC_log.txt";
  std::ifstream inStream(inFile, std::ios::binary);
  if (!inStream) {
    std::cerr << "File " << inFile << " not found" << std::endl;
    return false;
  }
  outStream.open(outFile,std::ios_base::out);
  if (!outStream) {
    std::cerr << "File " << outFile << " not found" << std::endl;
    return false;
  }
  ch.Add(inFile);
  expected_noise_per_chip*=ch.GetEntries();
  //Define the extraction tree
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
  //Clean the txt file
  outFile = DataDir+to_string(run)+"/extraction.txt";
  TString bash_command = "rm "+outFile;
  gSystem->Exec(bash_command);
  //Set ts number of thread
  gSystem->Exec("ts -S 40");

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
  if(no_TP) thr_ch = -1;
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
  case 2: 
    if(no_TP)  command = "tcoarse_min_ts>>"+h_name;
    if(!no_TP) command = "t_min_ttrigg>>"+h_name;
    return command;
  default:
    cout<<"No command has been selected"<<endl;
    return command;
  }
}
