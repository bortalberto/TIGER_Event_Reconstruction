#include "event.h"

const bool save_TP             = true ;
const bool at_least_two_roc    = false;
const bool test_chip_channel   = false;
const bool test_ROC_efficiency = true;
const bool DEBUG               = false;
const bool print_hits          = false;

const int  MAX_SIZE   = 50000;
const int  N_TIGER    =  150;
const int  MAX_EVENT  = 30000;
const int  wrong_TP   = -9999;

const int  count_min_cout = 332;
const int  count_max_cout = count_min_cout + MAX_EVENT;

int trigg_channel = 62;
//if(run>=118) trigg_channel=62;

//In file variables
int dchannel, dgemroc, dFEB, dcount, dtimestamp, dstrip_x, dstrip_v, dl1ts_min_tcoarse, dlasttigerframenum, dchip, dFEB_label, drunNo, dlayer, dtrigg_flag, dtcoarse_min_ts;
float dcharge_SH, dcharge_TOT, dpos_phi, dtcoarse, decoarse, dtfine, define, dttrigg, dtrigg_tcoarse, dconstant, dslope, dqmax, dtime, dradius, ddelta_coarse; 

//TP test
int TP_count[N_TIGER];
int TP_value[N_TIGER];
int TP_diff [N_TIGER];
double eff[N_TIGER], eff_T[N_TIGER], eff_Tpm1[N_TIGER];
int nhit_TIGER[N_TIGER];

TTree *otree, *ootree;

int count_unique(std::vector<int> v){
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
    return v.size();
}

int count_diff(std::vector<int> v){
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());

    return fabs(v.back() - v.front());
}

int moda(vector<int> v){

  if(print_hits) for(int i=0;i<v.size();i++) cout<<v[i]<<endl;
  int mode =wrong_TP;
  int count = 1;
  if(v.size()==1) return v[0];
  for(int j=0;j<v.size();j++){
    int number=v[j];
    for (int i = 1, countMode = 1; i < v.size(); ++i) {
      if (v[i] == number)

        ++countMode;
      if (countMode > count) {
        count = countMode;
        mode = number;
      }
    }
  }
  if(print_hits)cout<<"mode: "<<mode<<endl;  
  return mode;
}


void event(int run, int subrun){
  std::string iname=ANADIR;
  iname=iname+std::to_string(run)+"/Sub_RUN_ana_"  +std::to_string(subrun)+".root";
  std::string oname=ANADIR;
  oname=oname+std::to_string(run)+"/Sub_RUN_event_"+std::to_string(subrun)+".root";

  std::ifstream inStream(iname, std::ios::binary);
  if (!inStream) {
    std::cerr << "File " << iname << " not found" << std::endl;
    return;
  }



  //cout << "************* RUN " << run << endl;
  //cout << "************* TRIGGER CHANNEL " << trigg_channel << endl;
  if(run < 118) trigg_channel = 20;
  //cout << "************* post TRIGGER CHANNEL " << trigg_channel << endl;

  bool toohits = false;
  auto file = new TFile(iname.c_str());
  auto tree = (TTree*)file->Get("tree");
  
  int dchannel, dgemroc, dFEB, dcount, dtimestamp, dstrip_x, dstrip_v, dl1ts_min_tcoarse, dlasttigerframenum, dchip, dFEB_label, drunNo, dlayer, dtrigg_flag, dtac;
  float dcharge_SH, dhcarge_TOT, dpos_phi, dtcoarse, decoarse, dtfine, define, dttrigg, dtrigg_tcoarse, dconstant, dslope, dqmax, dtime, dradius, ddelta_coarse; 
  bool dsaturated;

  tree->SetBranchAddress("runNo",&drunNo);
  tree->SetBranchAddress("layer",&dlayer);
  tree->SetBranchAddress("channel",&dchannel);
  tree->SetBranchAddress("gemroc",&dgemroc);
  tree->SetBranchAddress("FEB",&dFEB);
  tree->SetBranchAddress("charge_SH",&dcharge_SH);
  tree->SetBranchAddress("charge_TOT",&dcharge_TOT);
  tree->SetBranchAddress("count",&dcount);
  tree->SetBranchAddress("timestamp",&dtimestamp);
  tree->SetBranchAddress("pos_phi",&dpos_phi);
  tree->SetBranchAddress("radius",&dradius);
  tree->SetBranchAddress("strip_v",&dstrip_v);
  tree->SetBranchAddress("strip_x",&dstrip_x);
  tree->SetBranchAddress("tcoarse",&dtcoarse);
  tree->SetBranchAddress("ecoarse",&decoarse);
  tree->SetBranchAddress("tfine",&dtfine);
  tree->SetBranchAddress("efine",&define);
  tree->SetBranchAddress("l1ts_min_tcoarse",&dl1ts_min_tcoarse);
  tree->SetBranchAddress("tcoarse_min_ts",&dtcoarse_min_ts);
  tree->SetBranchAddress("trigg_tcoarse",&dtrigg_tcoarse);
  tree->SetBranchAddress("lasttigerframenum",&dlasttigerframenum);
  tree->SetBranchAddress("chip",&dchip);
  tree->SetBranchAddress("FEB_label",&dFEB_label);
  tree->SetBranchAddress("QDCcali_constant",&dconstant);
  tree->SetBranchAddress("QDCcali_slope",&dslope);
  tree->SetBranchAddress("QDCcali_qmax",&dqmax);
  tree->SetBranchAddress("time",&dtime);
  tree->SetBranchAddress("trigg_flag",&dtrigg_flag);
  tree->SetBranchAddress("delta_coarse",&ddelta_coarse);
  tree->SetBranchAddress("saturated",&dsaturated);
  tree->SetBranchAddress("tac",&dtac);
  
  std::multimap<int,int> mevt;
  
  if(print_hits) {
    cout<<"*************************"<<endl;
    cout<<"*************************"<<endl;
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      if(dcount>=count_min_cout && dcount< count_max_cout) cout<<"index: "<<i<<" + ROC: "<<dgemroc<<" channel: "<<dchannel<<" count: "<<dcount<<endl;
    }
    
    cout<<"*************************"<<endl;
    cout<<"*************************"<<endl;
  }
  
  
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    mevt.insert(std::pair<int,int>(dcount,i));
  }
  
  
  if(test_chip_channel){
    //Channel analysis
    int max_channel = 64;
    int max_chip = N_TIGER;
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TH2I *channel_chip = new TH2I ("channel_chip","channel_chip",max_channel,0,max_channel,max_chip,0,max_chip);
    int mean =999999;
    int std = 999999;
    int count_cut = 0;
    int count_post = 0;
    double mean_post = 0;
    double std_post = 0;
    int bad_channel[max_chip][max_channel];
    int n_bad_channel[max_chip];
    int n_std=3; // it was 3
    float charge_min = 2.5;
    bool first_time=true;
    TCut bad_ch = "";
    TCut charge_cut = Form("charge_SH> %f && charge_SH<20 && channel!= %d && strip_x>0",charge_min,trigg_channel);
    for(int i=0;i<max_chip;i++) n_bad_channel[i]=0;
    tree->Draw("FEB_label*2+chip-1:channel>>channel_chip",charge_cut,"zcol");
    c1->SaveAs("channel_chip.pdf(","pdf");    
    while(1){
      //cout<<"hello"<<endl;
      for(int i=0;i<max_channel;i++){
	for(int j=0;j<max_chip;j++){
	  int tmp = channel_chip->GetBinContent(i+1,j+1);
	  //cout<<i<<" "<<j<<" "<<tmp<<endl;
	  if(tmp>mean+n_std*std) {
	    //cout<<tmp<<" "<<mean+std<<" "<<i<<" "<<j<<endl;
	    channel_chip->SetBinContent(i+1,j+1,0);
	    cout<<"Removed chip "<<j<<" channel "<<i<<endl;
	    bad_channel[j][n_bad_channel[j]]=i;
	    n_bad_channel[j]++;
	    //bad_ch += Form("(FEB_label+chip-1!=%i)||(FEB_label+chip-1==%i&&channel!=%i)",j,j,i);
	    //cout<<bad_ch<<endl;
	    count_cut++;
	  }
	  if(tmp && tmp<mean+n_std*std){
	    mean_post+=tmp;
	    count_post++;
	  }   
	}
      }
      mean_post/=count_post;
      for(int i=0;i<max_channel;i++){
	for(int j=0;j<max_chip;j++){
	  int tmp = channel_chip->GetBinContent(i+1,j+1);
	  //cout<<i<<" "<<j<<" "<<tmp<<endl;
	  if(tmp>0 && tmp<mean+n_std*std) {
	    std_post+=pow(mean_post-tmp,2);
	  }   
	}
      }
      std_post/=count_post;
      std_post=sqrt(std_post);
      //cout<<count_cut<<endl;
      cout<<count_post<<" "<<mean_post<<" "<<std_post<<endl;
      if(!count_cut && !first_time) break;
      first_time=false;
      count_cut=0;
      mean=mean_post;
      std=std_post;
      count_post=0;
      mean_post=0;
      std_post=0;
    }
    
    
    for(int j=0;j<max_chip;j++){
      cout<<"chip "<<j<<" n_bad_ch "<<n_bad_channel[j]<<endl;
      if(!n_bad_channel[j]) continue;
      TString tmp_cut = Form("(FEB_label+chip-1!=%i)||(FEB_label+chip-1==%i",j,j);
      for(int i=0;i<n_bad_channel[j];i++) tmp_cut +=Form("&&channel!=%i",bad_channel[j][i]);
      tmp_cut +=")";
      bad_ch += tmp_cut;
    }
    
    
    TH1D *h = new TH1D ("h","h",100,0,20);
    c1->SaveAs("channel_chip.pdf","pdf");
    //cout<<charge_cut<<" "<<charge_min<<endl;
    tree->Draw("charge_SH>>h",charge_cut);
    cout<<"Mean charge beofre cut: "<<h->GetMean()<<endl;
    c1->SaveAs("channel_chip.pdf","pdf");
    charge_cut += bad_ch;
    //cout<<charge_cut<<endl;
    tree->Draw("charge_SH>>h",charge_cut);
    cout<<"Mean charge after cut: "<<h->GetMean()<<endl;
    c1->SaveAs("channel_chip.pdf)","pdf");
  } 
 //return;
  
  
  
  auto ofile = new TFile(oname.c_str(),"RECREATE");
  otree = new TTree("tree","tree");
  //ootree = new TTree("tee","tee");

  int evtNo, nhits, ngemrocs, ntimestamp, runNo, ntcoarse_L1_TP, ntcoarse_L2_TP;
  float trigg_tcoarse;
  
  std::vector<int> tcount, tchannel, tgemroc, tFEB, ttimestamp, tstrip_x, tstrip_v, tl1ts_min_tcoarse, tchip, tFEB_label, tquality, tlayer, ttac, ttcoarse_min_ts;
  std::vector<float> tcharge_SH, tcharge_TOT, tpos_phi, ttcoarse, tecoarse, ttfine, tefine, t_min_ttrigg, tconstant, tslope, tqmax, ttime, tradius, ttrigg, delta_coarse;
  std::vector<bool> tsaturated;
  
  if(DEBUG) std::cout << "DEBUG::Just created all the vector variable" << std::endl;

  otree->Branch("runNo",&runNo,"runNo/i"); // the run index
  otree->Branch("evtNo",&evtNo,"evtNo/i"); // the event number index in each run
  otree->Branch("nhits",&nhits,"nhits/i"); // the number of hits in each event
  otree->Branch("ngemrocs",&ngemrocs,"ngemrocs/i"); // the number of gemrocs fired in each event
  otree->Branch("ntimestamp",&ntimestamp,"ntimestamp/i"); // the number of different GEMROC LOCAL_L1_TIMESTAMP values in each event
  otree->Branch("tcoarse_L1_TP_diff",&ntcoarse_L1_TP,"ntcoarse_L1_TP/i"); // the number of different TCOARSEs for test pulse in L1
  otree->Branch("tcoarse_L2_TP_diff",&ntcoarse_L2_TP,"ntcoarse_L2_TP/i"); // the number of different TCOARSEs for test pulse in L2
  otree->Branch("trigg_tcoarse",&trigg_tcoarse,"trigg_tcoarse/f"); // tcoarse value for the trigger channel in each event
  
  otree->Branch("local_l1_count"    , "vector<int>"  , &tcount           ); 
  otree->Branch("layer"             , "vector<int>"  , &tlayer           ); // Layer No. for each hit
  otree->Branch("channel"           , "vector<int>"  , &tchannel         ); // channel ID for each hit
  otree->Branch("gemroc"            , "vector<int>"  , &tgemroc          ); // GEMROC ID for each hit
  otree->Branch("FEB_SW"            , "vector<int>"  , &tFEB             ); // FEB ID for each hit
  otree->Branch("local_l1_timestamp", "vector<int>"  , &ttimestamp       ); // GEMROC LOCAL_L1_TIMESTAMP for each hit
  otree->Branch("strip_v"           , "vector<int>"  , &tstrip_v         ); // strip V number of each hit
  otree->Branch("strip_x"           , "vector<int>"  , &tstrip_x         ); // strip X number of each hit
  otree->Branch("chip"              , "vector<int>"  , &tchip            ); // FEB chip ID for each hit
  otree->Branch("FEB_label"         , "vector<int>"  , &tFEB_label       ); // FEB label number for each hit
  otree->Branch("l1ts_min_tcoarse"  , "vector<int>"  , &tl1ts_min_tcoarse); // GEMROC LOCAL_L1_TIMESTAMP - tcoarse for each hit
  otree->Branch("tcoarse_min_ts"    , "vector<int>"  , &ttcoarse_min_ts  ); // tcoarse - TIMESTAMP 
  otree->Branch("quality"           , "vector<int>"  , &tquality         ); // tags of good event for each hit
  otree->Branch("charge_SH"         , "vector<float>", &tcharge_SH       ); // charge in S&H mode, with QDC calibration, for each hit
  otree->Branch("charge_TOT"        , "vector<float>", &tcharge_TOT      ); 
  otree->Branch("radius"            , "vector<float>", &tradius          ); // radius for each hit, L1: 90.223 mm; L2: 129.8 mm
  otree->Branch("pos_phi"           , "vector<float>", &tpos_phi         ); // phi positon (strip_x) for each hit
  otree->Branch("tcoarse"           , "vector<float>", &ttcoarse         ); // tcoarse value for each hit
  otree->Branch("ecoarse"           , "vector<float>", &tecoarse         ); // ecoarse value for each hit
  otree->Branch("tfine"             , "vector<float>", &ttfine           ); // tfine value for each hit
  otree->Branch("efine"             , "vector<float>", &tefine           ); // efine value for each hit
  otree->Branch("ttrigg"            , "vector<float>", &ttrigg           ); // tcoarse of each trigger channel 
  otree->Branch("t_min_ttrigg"      , "vector<float>", &t_min_ttrigg     ); // tcoarse of each hit - tcoarse from trigger channel in each event
  otree->Branch("QDCcali_constant"  , "vector<float>", &tconstant        ); // constant value from DQC calibration curve for each hit
  otree->Branch("QDCcali_slope"     , "vector<float>", &tslope           ); // slope value from DQC calibration curve for each hit
  otree->Branch("QDCcali_qmax"      , "vector<float>", &tqmax            ); // Q max value from DQC calibration curve for each hit
  otree->Branch("time"              , "vector<float>", &ttime            );
  otree->Branch("delta_coarse"      , "vector<float>", &delta_coarse     );
  otree->Branch("count"             , "vector<int>"  , &tcount           ); // count number from the ROC
  otree->Branch("saturated"         , "vector<bool>" , &tsaturated       );
  otree->Branch("tac"               , "vector<int>"  , &ttac             );
  
  if(DEBUG) std::cout << "DEBUG::Just Saved all the vector branches" << std::endl;

  otree->Branch("TP_count"      , &TP_count   , "TP_count[200]/I");
  otree->Branch("TP_value"      , &TP_value   , "TP_value[200]/I");
  otree->Branch("TP_diff"       , &TP_diff    , "TP_diff[200]/I" );
  //ootree->Branch("eff"     , &eff     , "eff[200]/D");
  //ootree->Branch("eff_T"   , &eff_T   , "eff_T[200]/D");
  //ootree->Branch("eff_Tpm1", &eff_Tpm1, "eff_Tpm1[200]/D");

  std::multimap<int, int>::iterator it = mevt.begin();
  int itcount = it->first;
  nhits = 0;
  evtNo = 0;
  ngemrocs = 0;
  runNo = 0;
  std::vector<int> vgemrocs;
  std::vector<int> vtimestamp;
  std::vector<float> vttrigg;
  std::vector<int> vlasttigerframenum;
  std::vector<int> vtrigg_gemroc;
  
  std::vector<int> vl1tsmintcoarse_L1_TP;
  std::vector<int> vl1tsmintcoarse_L2_TP;
  for(int i=0; i<N_TIGER;i++) {
    TP_count[i]=0;
    nhit_TIGER[i]=0;
  }
  if(print_hits) {
    cout<<"Number of entries: "<<tree->GetEntries()<<endl;
    cout<<"+++++++++++++++++"<<endl;
    cout<<"+++++++++++++++++"<<endl;
    for (std::pair<int, int> elem : mevt){
      if(elem.first>= count_min_cout && elem.first < count_max_cout) cout<<elem.first<<" "<<elem.second<<endl;
      if(elem.first == count_max_cout) break;
    }
    cout<<"+++++++++++++++++"<<endl;
    cout<<"+++++++++++++++++"<<endl;
  }

  //  for (std::pair<int, int> elem : mevt){
  for(std::multimap<int, int>::iterator elem = mevt.begin(); elem != mevt.end(); elem++){
    //if(elem.second==0&&evtNo>2) break; COMMENTED
    toohits=false; 
    if(print_hits) {
      cout<<"elem: "<<elem->first<<" "<<elem->second<<" --> nhits: "<<nhits<<endl;
      cout<<"Bool itcount=elem.first = "<< (itcount==elem->first) <<endl;
      cout<<"Bool nhits==0           = "<< (nhits==0) << endl;
    }
    if(itcount==elem->first||nhits==0){ //an event collects all hits which have a same GEMROC LOCAL_COUNT number
      //	if(toohits)continue;
      
      itcount = elem->first;
      tree->GetEntry(elem->second);
      if(print_hits){
	cout<<"  itcount: "<<itcount<<endl;
	cout<<"    count: "<<dcount<<endl;
	cout<<"channel  : "<<dchannel<<endl;
	cout<<"ROC      : "<<dgemroc<<endl;
	cout<<"index    : "<<elem->second<<endl;
	cout<<"tcoarse  : "<<dtcoarse<<endl;
	cout<<"timestamp: "<<dtimestamp<<endl;
	cout<<"ts-tcoars: "<<dl1ts_min_tcoarse<<endl;
	cout<<endl;
      }
      if(!runNo) runNo = drunNo;
      
      //	  ttrigg[nhits] = -999;
      if(dtrigg_flag){ // trigger channgel
	trigg_tcoarse = dtrigg_tcoarse;
	vttrigg.push_back(dtrigg_tcoarse);
	vtrigg_gemroc.push_back(dgemroc);
	vlasttigerframenum.push_back(dlasttigerframenum);
	ttrigg.push_back(dtrigg_tcoarse);
	
      }
      else ttrigg.push_back(-999);
      
      //                                                                                                                                                                                                              
      //TEST ON TP

      TIGER_count(dFEB_label, dchip);
      TP_fill(dchannel, dFEB_label, dchip, dl1ts_min_tcoarse);
      
      //check test pulse                                                                                                                                                                                              
      if(dlayer==1){
	if(dchannel==trigg_channel) vl1tsmintcoarse_L1_TP.push_back(dl1ts_min_tcoarse);
      }
      if(dlayer==2){
	if(dchannel==trigg_channel) vl1tsmintcoarse_L2_TP.push_back(dl1ts_min_tcoarse);
      }
      
      vgemrocs.push_back(dgemroc);
      vtimestamp.push_back(dtimestamp);
      
      //To speed up the code the TP is not saved
      if(!save_TP && dchannel == trigg_channel) continue;
      
      //Check if there are BAD event due to ROC problem --> S/H only !!!
      if(ddelta_coarse!=25 && ddelta_coarse!=26) continue;
      
      tcount           .push_back(dcount           );
      tlayer           .push_back(dlayer           );
      tchannel         .push_back(dchannel         );
      tgemroc          .push_back(dgemroc          );
      tFEB             .push_back(dFEB             ); 
      ttimestamp       .push_back(dtimestamp       );
      tcharge_SH       .push_back(dcharge_SH       );
      tcharge_TOT      .push_back(dcharge_TOT      );
      tpos_phi         .push_back(dpos_phi         );
      tradius          .push_back(dradius          );
      tstrip_x         .push_back(dstrip_x         );
      tstrip_v         .push_back(dstrip_v         );
      ttcoarse         .push_back(dtcoarse         );
      tecoarse         .push_back(decoarse         );
      ttfine           .push_back(dtfine           );
      tchip            .push_back(dchip            );
      tFEB_label       .push_back(dFEB_label       );
      tefine           .push_back(define           );
      tl1ts_min_tcoarse.push_back(dl1ts_min_tcoarse);
      ttcoarse_min_ts  .push_back(dtcoarse_min_ts  );
      tconstant        .push_back(dconstant        );
      tslope           .push_back(dslope           );
      tqmax            .push_back(dqmax            );
      tquality         .push_back(0                );
      ttime            .push_back(dtime            ); 
      delta_coarse     .push_back(ddelta_coarse    );
      tsaturated       .push_back(dsaturated       );
      ttac             .push_back(dtac             );
      
      //if(DEBUG) std::cout << "DEBUG::Finished to push_back" << std::endl;
      nhits++;
      //	  if(nhits>MAX_SIZE && DEBUG) {std::cout<<"the number of hits is out of ARRAY range"<<std::endl; nhits=0; toohits=true; continue;}
    }
    else{
      ngemrocs = count_unique(vgemrocs);
      ntimestamp = count_unique(vtimestamp);
      
      if(!vl1tsmintcoarse_L1_TP.empty()) ntcoarse_L1_TP = count_diff(vl1tsmintcoarse_L1_TP);
      else ntcoarse_L1_TP = 9999;
      
      if(!vl1tsmintcoarse_L2_TP.empty()) ntcoarse_L2_TP = count_diff(vl1tsmintcoarse_L2_TP);
	else ntcoarse_L2_TP = 9999;
      
      if(at_least_two_roc && (ngemrocs<2 || nhits<2)){ //at least two gemrocs should be fired with two hits
	nhits = 0;
	vgemrocs.clear();
	vtimestamp.clear();
	vttrigg.clear();
	vl1tsmintcoarse_L1_TP.clear();
	vl1tsmintcoarse_L2_TP.clear();
	for(int i=0; i<N_TIGER;i++) {
	  TP_count[i]=0;
	  //nhit_TIGER[i]=0;
	}
	continue;
      }
      else{
	int good = 0;
	for(int i=0; i<nhits; i++){
	  good++; //no pos_phi requirement
	}
	if(good==0){
	  nhits = 0;
	  trigg_tcoarse = -9999;
	  vgemrocs.clear();
	  vtimestamp.clear();
	  vttrigg.clear();
	  
	  for(int i=0; i<N_TIGER;i++) {
	    TP_count[i]=0;
	    //nhit_TIGER[i]=0;
	  }
	  continue;
	}
	
	for(int i=0; i<nhits; i++){

	  if(tlayer.at(i)==0){
	    int TP_L0 = moda(vl1tsmintcoarse_L2_TP);
	    t_min_ttrigg.push_back(ttcoarse.at(i) - TP_L0);
	  }
	  else if(tlayer.at(i)==1){
	    int TP_L1 = moda(vl1tsmintcoarse_L1_TP);
	    if(TP_L1 == wrong_TP) TP_L1 = moda(vl1tsmintcoarse_L2_TP);
	    t_min_ttrigg.push_back(ttcoarse.at(i) - TP_L1);
	  }
	  else if(tlayer.at(i)==2){
	    int TP_L2 = moda(vl1tsmintcoarse_L2_TP);
	    if(TP_L2 == wrong_TP) TP_L2 = moda(vl1tsmintcoarse_L1_TP);
	    t_min_ttrigg.push_back(ttcoarse.at(i) - TP_L2);
	  } //if there are two trigger time in one event, only the first one is used here.
	  else {
	    int TP_L2 = moda(vl1tsmintcoarse_L2_TP);
	    t_min_ttrigg.push_back(ttcoarse.at(i) - TP_L2);
	  }
	}
      } 
      TP_test(vl1tsmintcoarse_L1_TP,vl1tsmintcoarse_L2_TP);    
      if(print_hits) {
	for(int i=0;i<nhits;i++){
	  cout<<"hits id: "<<i<<" ROC: "<<tgemroc[i]<<" channel: "<<tchannel[i]<<" time: "<<t_min_ttrigg[i]<<endl;
	}
	cout<<"L1 TP: "<<moda(vl1tsmintcoarse_L1_TP)<<endl;
	cout<<"L2 TP: "<<moda(vl1tsmintcoarse_L2_TP)<<endl;
	cout<<"Event end:   "<<evtNo<<endl;
	cout<<"---------------"<<endl;
      }
      evtNo++;
      otree->Fill();
      if(evtNo%100==0) std::cout << "Evt. No.\t" << evtNo << " \t & nHits \t" << nhits << std::endl;
      
      for(int i=0; i<N_TIGER;i++) {
	TP_count[i]=0;
	//nhit_TIGER[i]=0;
      }

      nhits=0;
      trigg_tcoarse = -9999;
      
      vgemrocs         .clear();
      vtimestamp       .clear();
      vttrigg          .clear();
      vl1tsmintcoarse_L1_TP   .clear();
      vl1tsmintcoarse_L2_TP   .clear();
      
      tcount           .clear();
      tlayer           .clear();
      tchannel         .clear();
      tgemroc          .clear();
      tFEB             .clear();
      ttimestamp       .clear();
      tstrip_v         .clear();
      tstrip_x         .clear();
      tchip            .clear();
      tFEB_label       .clear();
      tl1ts_min_tcoarse.clear();
      ttcoarse_min_ts  .clear();
      tquality         .clear();
      tcharge_SH       .clear();
      tcharge_TOT      .clear();
      tradius          .clear();
      tpos_phi         .clear();
      ttcoarse         .clear();
      tecoarse         .clear();
      ttfine           .clear();
      tefine           .clear();
      ttrigg           .clear();
      t_min_ttrigg     .clear();
      tconstant        .clear();
      tslope           .clear();
      tqmax            .clear();
      ttime            .clear();
      delta_coarse     .clear();
      tsaturated       .clear();
      ttac             .clear();
     
      if(evtNo >= MAX_EVENT) {
	std::cout << "***** ATTENTION THE MAXIMUM NUMBER OF EVENTS HAS BEEN REACHED, GOING TO END THE RUN. *****" << std::endl;
	std::cout << "***** \t MAX_EVENT \t" << MAX_EVENT << " \t *****"                                          << std::endl;
	break;
      }
      elem--;
    }//end of else not the same event
  }
  if(test_ROC_efficiency)TP_cout();    

  file->Close();
  otree->Write();
  ofile->Close();

  std::string ooname=ANADIR;  
  ooname=ooname+std::to_string(run)+"/Sub_RUN_TP_event_"+std::to_string(subrun)+".root";
  TFile *oofile = new TFile(ooname.c_str(),"RECREATE");
  ootree = new TTree("t1","t1");
  //ootree->Branch("eff"     , &eff     , "eff[200]/D");
  //ootree->Branch("eff_T"   , &eff_T   , "eff_T[200]/D");
  //ootree->Branch("eff_Tpm1", &eff_Tpm1, "eff_Tpm1[200]/D");
  //ootree->Fill();
  double TP_eff=0;
  int ichip=0;
  int iFEB=0;
  int ievent=0;
  int lsubRUN=subrun;
  int ihit=0;
  ootree->Branch("TP_eff"    , &TP_eff , "TP_eff/D");
  ootree->Branch("FEB_label" , &iFEB   , "FEB/I");
  ootree->Branch("chip"      , &ichip  , "chip/I");
  ootree->Branch("event"     , &ievent , "event/I");
  ootree->Branch("subRUN"    , &lsubRUN, "subRUN/I");
  ootree->Branch("nhit"      , &ihit   , "nhit/I");
  for(int iichip=0;iichip<N_TIGER;iichip++){
    TP_eff=eff_Tpm1[iichip];
    iFEB=iichip/2;
    ichip=iichip-(iichip/2*2)+1;
    ievent=evtNo;
    ihit=nhit_TIGER[iichip];
    ootree->Fill();
  }
  ootree->Write();
  oofile->Close();

  //file->Close();
  //ofile->Close();
}

void TIGER_count(int feb, int ip){
  if(feb==-1) return;
  int TIGER_ID = 2*feb;
  if(ip==2) TIGER_ID++;
  if(TIGER_ID<0 || TIGER_ID>N_TIGER){
    cout<<"Error in the TP test with TIGER ID: "<<TIGER_ID<<endl;
    return;
  }
  nhit_TIGER[TIGER_ID]++;
}

void TP_fill(int ch, int feb, int ip, double t){
  if(feb==-1) return;
  int TIGER_ID = 2*feb;
  if(ip==2) TIGER_ID++;
  if(TIGER_ID<0 || TIGER_ID>N_TIGER){
    cout<<"Error in the TP test with TIGER ID: "<<TIGER_ID<<endl;
    return;
  }
  if(ch == trigg_channel) {
    TP_count[TIGER_ID]++;
    TP_value[TIGER_ID]=t;
    if(DEBUG) std::cout << "\t TIGER_ID \t" << TIGER_ID << "\t TP_COUNT \t" << TP_count[TIGER_ID] << "\t TP_value \t" << TP_value[TIGER_ID] << std::endl;
  }
}

void TP_test(std::vector<int> v1, std::vector<int> v2){

  int TP_L1 = moda(v1);
  int TP_L2 = moda(v2);
  for(int i=0;i<32;i++){
    TP_diff[i]=TP_value[i]-TP_L1;
  }
  for(int i=32;i<N_TIGER;i++){
    TP_diff[i]=TP_value[i]-TP_L2;
  } 
}

void TP_cout(){

  TString cut, cot, cat;
  int n_TP_count_cut = 2;
  cout<<(double)(otree->GetEntries())<<endl;
  for(int i=0;i<N_TIGER;i++){
    eff[i]=eff_T[i]=eff_Tpm1[i]=-1;
  }
  for(int i=0;i<N_TIGER;i++){
    cot = Form("TP_count[%d]==1", i);
    cut = Form("TP_diff[%d]==0 && TP_count[%d]==1" , i, i);
    cat = Form("abs(TP_diff[%d])<%i && TP_count[%d]==1" , i, n_TP_count_cut, i);
    double num_eff      = (double)(otree->GetEntries(cot));
    double num_eff_T    = (double)(otree->GetEntries(cut));
    double num_eff_Tpm1 = (double)(otree->GetEntries(cat));
    double den          = (double)(otree->GetEntries(""));
    eff[i]      = num_eff/den;
    eff_T[i]    = num_eff_T/den;
    eff_Tpm1[i] = num_eff_Tpm1/den;
    cout<<"TIGER chip \t"<<i<<" efficiency: \t"<<eff[i]<<"\t and time efficiency: \t"<<eff_T[i]<<" and time efficiency +- "<<n_TP_count_cut-1<<" clock: \t"<<eff_Tpm1[i]<<endl;
  }
  return;
}
