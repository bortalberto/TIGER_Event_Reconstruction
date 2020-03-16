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
bool no_TP = true;

//PARAMETERS FOR THE EVENT/HITS SELECTION
const int INIT_FEB =  0;
const int MAX_FEB  = 55;  
const int N_FEB = MAX_FEB-INIT_FEB+1;
const int N_TIGER = 2*N_FEB;
const int N_CHIP = 2;
const int NCHANNEL   = 64;
const int charge_cut = -5; //fC
const int charge_max = 60; //fC
const int TIME_BIN   = 6.25; //ns   
const int signal_time_cut_min = -32;
const int signal_time_cut_max =  18;
const int noise_time_cut_min_left  = -120;
const int noise_time_cut_max_left  =  -70;
const int noise_time_cut_min_right =   80;//40;
const int noise_time_cut_max_right =   80;//80;
const int signal_time_cut_min_noTP      = -8937;//-1430*TIME_BIN;
const int signal_time_cut_max_noTP      = -8437;//-1350*TIME_BIN;
const int noise_time_cut_min_left_noTP  = -9375;//-1500*TIME_BIN;
const int noise_time_cut_max_left_noTP  = -9062;//-1450*TIME_BIN;
const int noise_time_cut_min_right_noTP = -8125;//-1300*TIME_BIN;//-1240;
const int noise_time_cut_max_right_noTP = -8125;//-1300*TIME_BIN;
const double time_window_noise=(noise_time_cut_max_left-noise_time_cut_min_left + noise_time_cut_max_right-noise_time_cut_min_right)*6.25e-9; //seconds                                                                   
const double time_window_noise_noTP=(noise_time_cut_max_left_noTP-noise_time_cut_min_left_noTP + noise_time_cut_max_right_noTP-noise_time_cut_min_right_noTP)*1e-9; //seconds                                                                   
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
TString inFile;
TString  inFileTP, outFile;
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
double t_noise_Hz, t_noise_Q, t_thr, t_thr_wid;
double t_time_start1,t_time_start2,t_time_start3,t_time_stop1,t_time_stop2,t_time_stop3,t_time_sigma1,t_time_sigma2,t_time_sigma3; 
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

void ext(vector<int> runs, int FEB_i, int chip_i, int ch_i){
  gStyle->SetOptFit(1111);
  bool print_here = false;
  for(int i=0;i<runs.size();i++){
    inFile=DataDir+to_string(runs.at(i))+"/event.root";
    ch.Add(inFile);
  }
  
  //
  //
  //noise                  
  int n_bin_noise=100;
  int n_bin_charge=0.5*(charge_max-charge_cut);
  TH1D *h_noise, *h_noise_charge;
  TString name_noise  = "h_noise";
  TString name_noise_charge = "h_noise_charge";
  //time distribution and noise rate
  if(no_TP){
    n_bin_noise= noise_time_cut_max_right_noTP - noise_time_cut_min_left_noTP;
    h_noise = new TH1D(name_noise,name_noise,n_bin_noise,noise_time_cut_min_left_noTP,noise_time_cut_max_right_noTP);
  }
  if(!no_TP) {
    n_bin_noise= noise_time_cut_max_right - noise_time_cut_min_left;
    h_noise = new TH1D(name_noise,name_noise,n_bin_noise,noise_time_cut_max_right,noise_time_cut_min_left);
  }
  if(print_here)cout<<Get_Command(2,name_noise)<<endl;
  if(print_here)cout<<Get_Cut(3,FEB_i,chip_i,ch_i)<<endl;
  ch.Draw(Get_Command(2,name_noise),Get_Cut(3,FEB_i,chip_i,ch_i),"goff");
  double nhitsofnoise = h_noise->GetEntries();
  double rate = 0;
  if(!no_TP) rate = nhitsofnoise/time_window_noise/ch.GetEntries();
  else       rate = nhitsofnoise/time_window_noise_noTP/ch.GetEntries();
  //charge distribution and noise charge
  h_noise_charge= new TH1D(name_noise_charge,name_noise_charge,n_bin_charge,charge_cut,charge_max);
  if(print_here)cout<<Get_Command(0,name_noise_charge)<<endl;
  if(print_here)cout<<Get_Cut(3,FEB_i,chip_i,ch_i)<<endl;
  ch.Draw(Get_Command(0,name_noise_charge),Get_Cut(3,FEB_i,chip_i,ch_i),"goff");
  double noise_charge = h_noise_charge->GetMean();

  //   
  //   
  //signal threshold
  int min_charge_hist=0;
  int max_charge_hist=30;
  double charge_resolution = 0.25; //fC                                                                                                                                                                                 
  int n_bin = (max_charge_hist-min_charge_hist)/charge_resolution;
  TString h_name="h_charge";
  TH1D *h_charge = new TH1D(h_name,h_name,n_bin,min_charge_hist,max_charge_hist);
  TF1 *f_thr = new TF1("f_thr","[0]/(1+exp(-(x-[1])/[2]))");
  if(print_here)cout<<Get_Command(0,h_name)<<endl;
  if(print_here)cout<<Get_Cut(4,FEB_i,chip_i,ch_i)<<endl;
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

  //
  //
  //signal time distribution
  TH1D *h_signal_time;
  int n_bin_time=100;
  TString h_name_time="h_signal_time";
  TF1 *f1 = new TF1("ciao1","[0] + [1]/(1+TMath::Exp(-(x-[2])/[3])) +  [1]/(1+TMath::Exp(-(x-[4])/(-[3]))) - [1]");
  TF1 *f2 = new TF1("ciao2","[0] + [1]*TMath::Exp(-[2]*(x-[3]))/(1+TMath::Exp(-(x-[4])/[5]))");
  TF1 *f3 = new TF1("ciao3","[0] + [1]/(1+TMath::Exp(-(x-[2])/(-[3])))");
  TF1 *f4 = new TF1("ciao4","[0] + [1]/(1+TMath::Exp(-(x-[2])/[3])) +  [1]/(1+TMath::Exp(-(x-[4])/(-[3]))) - [1] + gaus(5)");
  double t_start1,t_start2,t_start3,t_stop1,t_stop2,t_stop3,t_sigma1,t_sigma2,t_sigma3;
  t_start1=t_start2=t_start3=t_stop1=t_stop2=t_stop3=t_sigma1=t_sigma2=t_sigma3=0;
  if(no_TP){
    n_bin_time=signal_time_cut_max_noTP-signal_time_cut_min_noTP;
    h_signal_time = new TH1D(h_name_time,h_name_time,n_bin_time,signal_time_cut_min_noTP,signal_time_cut_max_noTP);
  }
  if(!no_TP){
    n_bin_time=signal_time_cut_max-signal_time_cut_min;
    h_signal_time = new TH1D(h_name_time,h_name_time,n_bin_time,signal_time_cut_min,signal_time_cut_max);
  }
  if(print_here)cout<<Get_Command(2,h_name_time)<<endl;
  if(print_here)cout<<Get_Cut(6,FEB_i,chip_i,ch_i)<<endl;
  ch.Draw(Get_Command(2,h_name_time),Get_Cut(6,FEB_i,chip_i,ch_i),"goff");
  if(print_here)cout<<"ENTRIES: "<<h_signal_time->GetEntries()<<endl;
  if(h_signal_time->GetEntries()){
    //fit time - function 1
    double par0=0;
    int i_count=0;
    for(int i=1;i<7;i++) if(h_signal_time->GetBinContent(i)){par0+=h_signal_time->GetBinContent(i);i_count++;}
    par0/=i_count;  
    if(no_TP) {
      f1->SetParameters(par0,0.5*h_signal_time->GetMaximum(),signal_time_cut_min_noTP+12,5,signal_time_cut_max_noTP-15);
      f1->SetParLimits(0,0,h_signal_time->GetMaximum());
      f1->SetParLimits(1,0,2*h_signal_time->GetMaximum());
      f1->SetParLimits(2,signal_time_cut_min_noTP,signal_time_cut_max_noTP);
      f1->SetParLimits(3,0,100);
      f1->SetParLimits(4,signal_time_cut_min_noTP,signal_time_cut_max_noTP);
    }
    if(!no_TP){
      f1->SetParameters(par0,0.5*h_signal_time->GetMaximum(),signal_time_cut_min+12,5,signal_time_cut_max-15);
      f1->SetParLimits(0,0,h_signal_time->GetMaximum());
      f1->SetParLimits(1,0,2*h_signal_time->GetMaximum());
      f1->SetParLimits(2,signal_time_cut_min,signal_time_cut_max);
      f1->SetParLimits(3,0,100);
      f1->SetParLimits(4,signal_time_cut_min,signal_time_cut_max);
    }
    f1->SetParLimits(0,0,0.2*h_signal_time->GetMaximum());
    f1->SetParLimits(1,0.2*h_signal_time->GetMaximum(),h_signal_time->GetMaximum());
    h_signal_time->Fit(f1,"Q");
    if(((gMinuit->fCstatu)=="CONVERGED ") || ((gMinuit->fCstatu)=="SUCCESSFUL") || ((gMinuit->fCstatu)=="OK ")|| ((gMinuit->fCstatu)=="CALL LIMIT"))        {
      t_start1 = f1->GetParameter(2);
      t_stop1  = f1->GetParameter(4);
      t_sigma1 = f1->GetParameter(3);
    }
    //fit time - function 2
    f2->SetParameters(f1->GetParameter(0),f1->GetParameter(1),0.05,f1->GetParameter(2),f1->GetParameter(2),f1->GetParameter(3));
    f2->SetParLimits(0,0,3*f1->GetParameter(0));
    f2->SetParLimits(1,0,2*h_signal_time->GetMaximum());
    f2->SetParLimits(2,0.,1);
    f2->SetParLimits(3,f1->GetParameter(2),f1->GetParameter(2)+7*f1->GetParameter(3));
    f2->SetParLimits(4,f1->GetParameter(2)-5*f1->GetParameter(3),f1->GetParameter(2)+5*f1->GetParameter(3));
    f2->SetParLimits(5,0,5*f1->GetParameter(2));
    h_signal_time->Fit(f2,"Q","",signal_time_cut_min_noTP,0.5*(f1->GetParameter(2)+f1->GetParameter(4)));
    if(((gMinuit->fCstatu)=="CONVERGED ") || ((gMinuit->fCstatu)=="SUCCESSFUL") || ((gMinuit->fCstatu)=="OK ")|| ((gMinuit->fCstatu)=="CALL LIMIT"))        {
      t_start2 = f2->GetParameter(4);
      t_sigma2 = f2->GetParameter(5);
    }
    //fit time - function 3
    f3->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(4),f1->GetParameter(3));
    f3->SetParLimits(0,0,3*f1->GetParameter(0));
    f3->SetParLimits(1,0,2*h_signal_time->GetMaximum());
    f3->SetParLimits(2,f1->GetParameter(4)-5*f1->GetParameter(3),f1->GetParameter(4)+5*f1->GetParameter(3));
    f3->SetParLimits(3,0,5*f1->GetParameter(3));
    h_signal_time->Fit(f3,"Q","",0.5*(f1->GetParameter(2)+f1->GetParameter(4)),signal_time_cut_max_noTP);
    if(((gMinuit->fCstatu)=="CONVERGED ") || ((gMinuit->fCstatu)=="SUCCESSFUL") || ((gMinuit->fCstatu)=="OK ")|| ((gMinuit->fCstatu)=="CALL LIMIT"))        {
      t_stop2  = f3->GetParameter(2);
    }
    //fit time - function 4
    f4->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(3),f1->GetParameter(4),0.1*f1->GetParameter(0),f1->GetParameter(2)+3*f1->GetParameter(3),f1->GetParameter(3));
    f4->SetParLimits(0,0.9*f1->GetParameter(0),1.1*f1->GetParameter(0));
    f4->SetParLimits(1,0.9*f1->GetParameter(1),1.1*f1->GetParameter(1));
    f4->SetParLimits(2,0.9*f1->GetParameter(2),1.1*f1->GetParameter(2));
    f4->SetParLimits(3,0.9*f1->GetParameter(3),1.1*f1->GetParameter(3));
    f4->SetParLimits(4,0.9*f1->GetParameter(4),1.1*f1->GetParameter(4));
    f4->SetParLimits(5,0,2*f1->GetParameter(0));
    f4->SetParLimits(6,f1->GetParameter(2)-3*f1->GetParameter(3),f1->GetParameter(2)+3*f1->GetParameter(3));
    f4->SetParLimits(7,0,10*f1->GetParameter(3));
    h_signal_time->Fit(f4,"Q");
    if(((gMinuit->fCstatu)=="CONVERGED ") || ((gMinuit->fCstatu)=="SUCCESSFUL") || ((gMinuit->fCstatu)=="OK ")|| ((gMinuit->fCstatu)=="CALL LIMIT"))        {
      t_start3 = f4->GetParameter(2);
      t_stop3  = f4->GetParameter(4);
      t_sigma3 = f4->GetParameter(5);
    }
  }
  //   
  //   
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
  //
  //
  //Write the output of the single strip on the file
  if(print_here)cout<<"FEB: "<<FEB_i<<" Chip: "<<chip_i<<" Channel: "<<ch_i<<" Rate: "<<rate<<" Thr: "<<threshold<<" StripX: "<<mx<<" StripV: "<<mv<<endl;
  outFile = DataDir+to_string(runs.at(0))+"/extraction.txt";
  outStream.open(outFile,std::ios::out |std::ios::app);
  outStream<<FEB_i<<" "<<chip_i<<" "<<ch_i<<" "<<rate<<" "<<noise_charge<<" "<<threshold<<" "<<thr_wid<<" "<<t_start1<<" "<<t_start2<<" "<<t_start3<<" "<<t_stop1<<" "<<t_stop2<<" "<<t_stop3<<" "<<t_sigma1<<" "<<t_sigma2<<" "<<t_sigma3<<" "<<mx<<" "<<mv<<endl;
  if(1 || print_here) cout<<FEB_i<<" "<<chip_i<<" "<<ch_i<<" "<<rate<<" "<<noise_charge<<" "<<threshold<<" "<<thr_wid<<" "<<t_start1<<" "<<t_start2<<" "<<t_start3<<" "<<t_stop1<<" "<<t_stop2<<" "<<t_stop3<<" "<<t_sigma1<<" "<<t_sigma2<<" "<<t_sigma3<<" "<<mx<<" "<<mv<<endl;
  outStream.close();
  return;
}

int count_line(int run){
  outFile =  DataDir+to_string(run)+"/extraction.txt";
  int numLines = 0;
  ifstream in(outFile);                                                                                          
  std::string unused;                                                                                                                                                                     
  while ( std::getline(in, unused) ) ++numLines; 
  return numLines;
}

void ext(int run){
  vector<int> runs;
  runs.push_back(run);
  ext(runs);
  return;
}

void ext(vector<int> runs){
  //cout<<"++ QAQC analysis and plot"<<endl;
  if(init(runs)){
    cout<<"++ Initialization completed"<<endl;
    if(check_L1)     ext_i(runs, 0,15); //L1
    if(check_L2)     ext_i(runs,16,43); //L2
    if(check_planar) ext_i(runs,44,55); //planari
  }
  else return;
  //cout<<"++ Terminated QAQC"<<endl;
  cout<<"Number of running channels: "<<n_runs<<endl;
  while(count_line(runs.at(0))<n_runs) {
    cout<<"Number of analyzed channels: "<<count_line(runs.at(0))<<endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(10000));
  }
  //Read the output from each channel and write it on a file
  outFile = DataDir+to_string(runs.at(0))+"/extraction.txt";
  ifstream in(outFile); 
  while(in>>t_feb>>t_chip>>t_ch>>t_noise_Hz>>t_noise_Q>>t_thr>>t_thr_wid>>t_time_start1>>t_time_start2>>t_time_start3>>t_time_stop1>>t_time_stop2>>t_time_stop3>>t_time_sigma1>>t_time_sigma2>>t_time_sigma3>>t_strip_x>>t_strip_v) {
    if(t_noise_Hz==0) t_thr=t_thr_wid=t_time_start1=t_time_start2=t_time_start3=t_time_stop1=t_time_stop2=t_time_stop3=t_time_sigma1=t_time_sigma2=t_time_sigma3=0; 
    otree->Fill();
  }
  otree->Write(); 
  ofile->Close();
  //Copy ttree if other runs
  if(runs.size()>1){
    TString source_name=Form("/home/ihep_data/data/raw_daq/extracted_noise_thr_%d.root",runs.at(0));
    for(int i=1;i<runs.size();i++){
      TString destin_name=Form("/home/ihep_data/data/raw_daq/extracted_noise_thr_%d.root",runs.at(1));
      TString bash_command = "cp " + source_name + " " + destin_name;
      cout<<bash_command<<endl;
      gSystem->Exec(bash_command);
    }
  }
  //Clean
  gSystem->Exec("rm -f /tmp/ihep_data/ts*");
  return;
}

void ext_i(vector<int> runs, int start_FEB, int end_FEB){
  for(int FEB_i=start_FEB;FEB_i<=end_FEB;FEB_i++){
    if(FEB_i==44 || FEB_i==45 || FEB_i==46 || FEB_i==47 || FEB_i==48 || FEB_i==50 || FEB_i==51) continue;
    for(int chip_i=1; chip_i<3; chip_i++){
      for(int ch_i=0;ch_i<NCHANNEL;ch_i++){
	TString bash_command;
	if(runs.size()==1) bash_command=Form("ts -n $exe_ter -X %d %d %d %d %d",runs.size(),runs.at(0),FEB_i,chip_i,ch_i);
	if(runs.size()==2) bash_command=Form("ts -n $exe_ter -X %d %d %d %d %d %d",runs.size(),runs.at(0),runs.at(1),FEB_i,chip_i,ch_i);
        if(runs.size()==3) bash_command=Form("ts -n $exe_ter -X %d %d %d %d %d %d %d",runs.size(),runs.at(0),runs.at(1),runs.at(2),FEB_i,chip_i,ch_i);
        if(runs.size()==4) bash_command=Form("ts -n $exe_ter -X %d %d %d %d %d %d %d %d",runs.size(),runs.at(0),runs.at(1),runs.at(2),runs.at(3),FEB_i,chip_i,ch_i);
	//cout<<bash_command<<endl;
	gSystem->Exec(bash_command);
	n_runs++;
      }
    }
  }
  return;
}

bool init(vector<int> runs){
  bool isok=runs.size();
  for(int i=0;i<runs.size();i++){
    cout<<"test: "<<i+1<<" --> "<<runs.at(i)<<endl;
    isok*=init(runs.at(i),runs.at(0));
  }
  return isok;
}

bool init(int run, int run0){
  //open the root file
  inFile = DataDir+to_string(run)+"/event.root";
  std::ifstream inStream(inFile, std::ios::binary);
  if (!inStream) {
    std::cerr << "File " << inFile << " not found" << std::endl;
    return false;
  }
  TFile *file1 = new TFile((TString)inFile);
  TTree *tree1 = (TTree*)file1->Get("tree");
  bool b1(0),b2(0),b3(0);
  b1=file1->IsOpen();
  b2=(bool)tree1;
  if(b2) b3=tree1->IsFolder();
  if(b1*b2*b3 == 0) {
    cout<<"ERROR in file "<<inFile<<" "<<b1<<" "<<b2<<" "<<b3<<endl;
    return false;
  }
  ch.Add((TString)inFile);
  if(ch.GetEntries()==0){
    std::cerr << "No entries in the event file RUN "<< run<<endl;
  }
  expected_noise_per_chip*=ch.GetEntries();
  //Define the extraction tree
  ofile_name=Form("/home/ihep_data/data/raw_daq/extracted_noise_thr_%d.root",run0);
  ofile = new TFile(ofile_name,"RECREATE");
  otree = new TTree("tree","tree");
  //General info
  otree->Branch("FEB",&t_feb,"FEB/I");
  otree->Branch("chip",&t_chip,"chip/I");
  otree->Branch("channel",&t_ch,"channel/I");
  //Mapping
  otree->Branch("strip_x",&t_strip_x,"strip_x/I");
  otree->Branch("strip_v",&t_strip_v,"strip_v/I");
  //Noise
  otree->Branch("noise_Hz",&t_noise_Hz,"noise_Hz/D");
  otree->Branch("noise_Q",&t_noise_Q,"noise_Q/D");
  //Threshold
  otree->Branch("threshold",&t_thr,"threshold/D");
  otree->Branch("thr_width",&t_thr_wid,"thr_width/D");
  //Time
  otree->Branch("time_start1",&t_time_start1,"time_start1_ns/D");
  otree->Branch("time_start2",&t_time_start2,"time_start2_ns/D");
  otree->Branch("time_start3",&t_time_start3,"time_start3_ns/D");
  otree->Branch("time_stop1" ,&t_time_stop1, "time_stop1_ns/D");
  otree->Branch("time_stop2" ,&t_time_stop2, "time_stop2_ns/D");
  otree->Branch("time_stop3" ,&t_time_stop3, "time_stop3_ns/D");
  otree->Branch("time_sigma1",&t_time_sigma1,"time_sigma1_ns/D");
  otree->Branch("time_sigma2",&t_time_sigma2,"time_sigma2_ns/D");
  otree->Branch("time_sigma3",&t_time_sigma3,"time_sigma3_ns/D");


  //Clean the txt file
  outFile = DataDir+to_string(run0)+"/extraction.txt";
  TString bash_command = "rm -f "+outFile;
  gSystem->Exec(bash_command);
  //Set ts number of thread
  gSystem->Exec("export TS_MAXCONN=100");
  gSystem->Exec("TS_MAXCONN=100");
  gSystem->Exec("ts -S 40");
  gSystem->Exec("sleep 1");
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
    //Noise out of time selection - feb - chip
    if(!no_TP) cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&((t_min_ttrigg>=%d && t_min_ttrigg<%d)|| (t_min_ttrigg>=%d&&t_min_ttrigg<%d))&&charge_SH>%d",thr_ch,var1,var2,noise_time_cut_min_left,noise_time_cut_max_left,noise_time_cut_min_right,noise_time_cut_max_right,charge_cut);
    //if(no_TP) cut = Form("FEB_label==%d&&chip==%d&&((t_min_ttrigg>=%d && t_min_ttrigg<%d)|| (t_min_ttrigg>=%d&&t_min_ttrigg<%d))&&charge_SH>%d",var1,var2,noise_time_cut_min_left,noise_time_cut_max_left,noise_time_cut_min_right,noise_time_cut_max_right,charge_cut);
    if(no_TP) cut = Form("FEB_label==%d&&chip==%d&&((time_ns>=%d && time_ns<%d)|| (time_ns>=%d&&time_ns<%d))&&charge_SH>%d",var1,var2,noise_time_cut_min_left_noTP,noise_time_cut_max_left_noTP,noise_time_cut_min_right_noTP,noise_time_cut_max_right_noTP,charge_cut);
    return cut;
  case 1:
    //Signal in time selection - feb - chip
    if(!no_TP) cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&t_min_ttrigg>=%d && t_min_ttrigg<%d&&charge_SH>%d",thr_ch,var1,var2,signal_time_cut_min,signal_time_cut_max,charge_cut);
    if(no_TP)  cut = Form("FEB_label==%d&&chip==%d&&time_ns>=%d && time_ns<%d && charge_SH>%d",var1,var2,signal_time_cut_min_noTP,signal_time_cut_max_noTP,charge_cut);
    return cut;
  case 2:
    //TP comunication selection - feb - chip
    cut = Form("FEB==%d && chip==%d", var1, var2);
    return cut;
  case 3:
    //Noise out of time selection on one channel - feb - chip - channel
    if(!no_TP) cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&((t_min_ttrigg>=%d && t_min_ttrigg<%d)|| (t_min_ttrigg>=%d&&t_min_ttrigg<%d))&&charge_SH>%d && channel==%d",thr_ch,var1,var2,noise_time_cut_min_left,noise_time_cut_max_left,noise_time_cut_min_right,noise_time_cut_max_right,charge_cut,var3);
    if(no_TP)  cut = Form("FEB_label==%d&&chip==%d&&((time_ns>=%d && time_ns<%d)|| (time_ns>=%d&&time_ns<%d))&&charge_SH>%d && channel==%d",var1,var2,noise_time_cut_min_left_noTP,noise_time_cut_max_left_noTP,noise_time_cut_min_right_noTP,noise_time_cut_max_right_noTP,charge_cut,var3);
    return cut;
  case 4:
    //Channel selection - feb - chip - channel
    if(no_TP) cut = Form("FEB_label==%d&&chip==%d&&charge_SH>%d && channel==%d",var1,var2,charge_cut,var3);
    if(!no_TP) cut = Form("channel!=%d && FEB_label==%d && chip==%d && charge_SH>%d && channel==%d",thr_ch,var1,var2,charge_cut,var3);
    return cut;
  case 5:
    //TP comunication selection - min_eff - feb - chip
    cut = Form("TP_eff<%f && FEB==%d && chip==%d", min_efficiency, var1, var2);
    return cut;
  case 6:
    //Signal in time selection - feb - chip - channel
    if(!no_TP) cut = Form("channel!=%d&&FEB_label==%d&&chip==%d&&t_min_ttrigg>=%d && t_min_ttrigg<%d&&charge_SH>%d&&channel==%d",thr_ch,var1,var2,signal_time_cut_min,signal_time_cut_max,charge_cut,var3);
    if(no_TP)  cut = Form("FEB_label==%d&&chip==%d&&time_ns>=%d && time_ns<%d && charge_SH>%d&&channel==%d",var1,var2,signal_time_cut_min_noTP,signal_time_cut_max_noTP,charge_cut,var3);
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
  //Charge distribution
  case 0:
    command = "charge_SH>>"+h_name;
    return command;
  //TP comunication efficiency
  case 1:
    command = "TP_eff>>"+h_name;
    return command;
  //Time distribution
  case 2: 
    //if(no_TP)  command = "tcoarse_min_ts-tfine>>"+h_name;
    if(no_TP)  command = "time_ns>>"+h_name;
    if(!no_TP) command = "t_min_ttrigg-tfine>>"+h_name;
    return command;
  //None
  default:
    cout<<"No command has been selected"<<endl;
    return command;
  }
}
