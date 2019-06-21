#include "ana.h"
void ana(int run, int subrun){
  std::string iname=ANADIR;
  iname=iname+std::to_string(run)+"/Sub_RUN_dec_"+std::to_string(subrun)+".root";
  std::string oname=ANADIR;
  oname=oname+std::to_string(run)+"/Sub_RUN_ana_"+std::to_string(subrun)+".root";

  int trigg_channel=20;
  if(run>=118) trigg_channel=62;
  bool save_TP = true;
  //int trigg_FEB=4;
  //int trigg_gemroc=4;
  auto mapfile = new TFile("mapping_IHEP.root");
  auto maptree = (TTree*)mapfile->Get("tree");
  
  int channel_id, gemroc_id, SW_FEB_id, pos_x, pos_v, chip_id, FEB_label_id;
  float phi;
  
  maptree->SetBranchAddress("channel_id", &channel_id);
  maptree->SetBranchAddress("gemroc_id", &gemroc_id);
  maptree->SetBranchAddress("SW_FEB_id", &SW_FEB_id);
  maptree->SetBranchAddress("pos_x", &pos_x);
  maptree->SetBranchAddress("pos_v", &pos_v);
  maptree->SetBranchAddress("phi", &phi);
  maptree->SetBranchAddress("chip_id", &chip_id);
  maptree->SetBranchAddress("FEB_label", &FEB_label_id);
  
  int mx[11][8][64];
  float mphi[11][8][64];
  int mv[11][8][64];
  int mchip_id[11][8][64];
  int mFEB_label_id[11][8][64];
  memset(mx, -1, sizeof(mx));
  memset(mv, -1, sizeof(mv));
  memset(mphi, -1, sizeof(mphi));
  memset(mchip_id, -1, sizeof(mchip_id));
  memset(mFEB_label_id, -1, sizeof(mFEB_label_id));
  
  for (int i = 0; i < maptree->GetEntries(); i++) {
    maptree->GetEntry(i);
    mx[gemroc_id][SW_FEB_id][channel_id] = pos_x;
    mv[gemroc_id][SW_FEB_id][channel_id] = pos_v;
    mphi[gemroc_id][SW_FEB_id][channel_id] = phi;
    mchip_id[gemroc_id][SW_FEB_id][channel_id] = chip_id;
    mFEB_label_id[gemroc_id][SW_FEB_id][channel_id] = FEB_label_id;
  }
  
  auto consfile = new TFile("QDCcalib.root");
  auto constree = (TTree*)consfile->Get("tree");
  
  int cons_channel_id, cons_gemroc_id, cons_SW_FEB_id;
  float cons_constant, cons_slope, cons_qmax;
  
  constree->SetBranchAddress("channel_id", &cons_channel_id);
  constree->SetBranchAddress("gemroc_id", &cons_gemroc_id);
  constree->SetBranchAddress("SW_FEB_id", &cons_SW_FEB_id);
  constree->SetBranchAddress("constant", &cons_constant);
  constree->SetBranchAddress("slope", &cons_slope);
  constree->SetBranchAddress("qmax", &cons_qmax);
  
  float cons[11][8][64];
  float mslope[11][8][64];
  float mqmax[11][8][64];
  memset(cons, -1, sizeof(cons));
  memset(mslope, -1, sizeof(mslope));
  memset(mqmax, -1, sizeof(mqmax));
  
  for (int i = 0; i < constree->GetEntries(); i++) {
    constree->GetEntry(i);
    if(cons_gemroc_id<0||cons_SW_FEB_id<-1) continue;
    cons[cons_gemroc_id][cons_SW_FEB_id][cons_channel_id] = cons_constant;
    mslope[cons_gemroc_id][cons_SW_FEB_id][cons_channel_id] = cons_slope;
    mqmax[cons_gemroc_id][cons_SW_FEB_id][cons_channel_id] = cons_qmax;
  }
  
  
  auto TDCconsfile = new TFile("TDCcalib.root");
  auto TDCconstree = (TTree*)TDCconsfile->Get("tree");
  
  int TDCcons_channel_id, TDCcons_gemroc_id, TDCcons_SW_FEB_id;
  float Tac0_Tfine_min, Tac0_Tfine_max, Tac1_Tfine_min, Tac1_Tfine_max, Tac2_Tfine_min, Tac2_Tfine_max, Tac3_Tfine_min, Tac3_Tfine_max, Tac0_Efine_min, Tac0_Efine_max, Tac1_Efine_min, Tac1_Efine_max, Tac2_Efine_min, Tac2_Efine_max, Tac3_Efine_min, Tac3_Efine_max;;
  
  TDCconstree->SetBranchAddress("channel_id", &TDCcons_channel_id);
  TDCconstree->SetBranchAddress("gemroc_id", &TDCcons_gemroc_id);
  TDCconstree->SetBranchAddress("SW_FEB_id", &TDCcons_SW_FEB_id);
  TDCconstree->SetBranchAddress("Tac0_Tfine_min", &Tac0_Tfine_min);
  TDCconstree->SetBranchAddress("Tac0_Tfine_max", &Tac0_Tfine_max);
  TDCconstree->SetBranchAddress("Tac1_Tfine_min", &Tac1_Tfine_min);
  TDCconstree->SetBranchAddress("Tac1_Tfine_max", &Tac1_Tfine_max);
  TDCconstree->SetBranchAddress("Tac2_Tfine_min", &Tac2_Tfine_min);
  TDCconstree->SetBranchAddress("Tac2_Tfine_max", &Tac2_Tfine_max);
  TDCconstree->SetBranchAddress("Tac3_Tfine_min", &Tac3_Tfine_min);
  TDCconstree->SetBranchAddress("Tac3_Tfine_max", &Tac3_Tfine_max);
  TDCconstree->SetBranchAddress("Tac0_Efine_min", &Tac0_Efine_min);
  TDCconstree->SetBranchAddress("Tac0_Efine_max", &Tac0_Efine_max);
  TDCconstree->SetBranchAddress("Tac1_Efine_min", &Tac1_Efine_min);
  TDCconstree->SetBranchAddress("Tac1_Efine_max", &Tac1_Efine_max);
  TDCconstree->SetBranchAddress("Tac2_Efine_min", &Tac2_Efine_min);
  TDCconstree->SetBranchAddress("Tac2_Efine_max", &Tac2_Efine_max);
  TDCconstree->SetBranchAddress("Tac3_Efine_min", &Tac3_Efine_min);
  TDCconstree->SetBranchAddress("Tac3_Efine_max", &Tac3_Efine_max);
  
  float TDCcons_Tmin[11][8][64][4];
  float TDCcons_Tbin[11][8][64][4];
  float TDCcons_Emin[11][8][64][4];
  float TDCcons_Ebin[11][8][64][4];
  memset(TDCcons_Tmin, -1, sizeof(TDCcons_Tmin));
  memset(TDCcons_Tbin, -1, sizeof(TDCcons_Tbin));
  memset(TDCcons_Emin, -1, sizeof(TDCcons_Emin));
  memset(TDCcons_Ebin, -1, sizeof(TDCcons_Ebin));
  
  for (int i = 0; i < constree->GetEntries(); i++) {
    TDCconstree->GetEntry(i);
    if(TDCcons_gemroc_id<0||TDCcons_SW_FEB_id<-1) continue;
    TDCcons_Tmin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][0] = Tac0_Tfine_min;
    TDCcons_Tmin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][1] = Tac1_Tfine_min;
    TDCcons_Tmin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][2] = Tac2_Tfine_min;
    TDCcons_Tmin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][3] = Tac3_Tfine_min;
    TDCcons_Tbin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][0] = 6.25 / (Tac0_Tfine_max - Tac0_Tfine_min);
    TDCcons_Tbin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][1] = 6.25 / (Tac1_Tfine_max - Tac1_Tfine_min);
    TDCcons_Tbin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][2] = 6.25 / (Tac2_Tfine_max - Tac2_Tfine_min);
    TDCcons_Tbin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][3] = 6.25 / (Tac3_Tfine_max - Tac3_Tfine_min);
    TDCcons_Emin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][0] = Tac0_Efine_min;
    TDCcons_Emin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][1] = Tac1_Efine_min;
    TDCcons_Emin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][2] = Tac2_Efine_min;
    TDCcons_Emin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][3] = Tac3_Efine_min;
    TDCcons_Ebin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][0] = 6.25 / (Tac0_Efine_max - Tac0_Efine_min);
    TDCcons_Ebin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][1] = 6.25 / (Tac1_Efine_max - Tac1_Efine_min);
    TDCcons_Ebin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][2] = 6.25 / (Tac2_Efine_max - Tac2_Efine_min);
    TDCcons_Ebin[TDCcons_gemroc_id][TDCcons_SW_FEB_id][TDCcons_channel_id][3] = 6.25 / (Tac3_Efine_max - Tac3_Efine_min);
  }
  
  
  auto datafile = new TFile(iname.c_str());
  auto datatree = (TTree*)datafile->Get("tree");
  
  int dchannel, dgemroc, dFEB, dcount, dtimestamp, dl1ts_min_tcoarse, dlasttigerframenum, dtac, drunNo, dlayer;
  float dcharge_SH, ddelta_coarse, dtcoarse, decoarse, dtfine, define;
  
  datatree->SetBranchAddress("layer_id", &dlayer);
  datatree->SetBranchAddress("runNo", &drunNo);
  datatree->SetBranchAddress("channel_id", &dchannel);
  datatree->SetBranchAddress("gemroc_id", &dgemroc);
  datatree->SetBranchAddress("tiger_id", &dFEB);
  datatree->SetBranchAddress("charge_SH", &dcharge_SH);
  datatree->SetBranchAddress("delta_coarse", &ddelta_coarse);
  datatree->SetBranchAddress("count", &dcount);
  datatree->SetBranchAddress("timestamp", &dtimestamp);
  datatree->SetBranchAddress("tcoarse", &dtcoarse);
  datatree->SetBranchAddress("ecoarse", &decoarse);
  datatree->SetBranchAddress("tfine", &dtfine);
  datatree->SetBranchAddress("efine", &define);
  datatree->SetBranchAddress("l1ts_min_tcoarse", &dl1ts_min_tcoarse);
  datatree->SetBranchAddress("lasttigerframenum", &dlasttigerframenum);
  datatree->SetBranchAddress("tac_id", &dtac);
  
  
  auto ofile = new TFile(oname.c_str(),"RECREATE");
  auto otree = new TTree("tree","tree");
  
  int channel, gemroc, FEB, strip_x, strip_v, count, timestamp, l1ts_min_tcoarse, lasttigerframenum, chip, FEB_label, tac, runNo, layer, max_count, trigg_flag;
  float charge_SH, charge_SH_uncal, constant, slope, qmax, delta_coarse, pos_phi, tcoarse, ecoarse, tfine, efine, ttrigg = -99999, trigg_tcoarse = -99999, tfine_uncal, efine_uncal, time, radius; 
  
  
  otree->Branch("runNo",&runNo,"runNo/I");
  otree->Branch("layer",&layer,"layer/I");
  otree->Branch("channel",&channel,"channel/I");
  otree->Branch("gemroc",&gemroc,"gemroc/I");
  otree->Branch("FEB",&FEB,"FEB/I");
  otree->Branch("max_count",&max_count,"max_count/I");
  otree->Branch("strip_x",&strip_x,"strip_x/I");
  otree->Branch("strip_v",&strip_v,"strip_v/I");
  otree->Branch("radius",&radius,"radius/F");
  otree->Branch("charge_SH_uncal",&charge_SH_uncal,"charge_SH_uncal/F");
  otree->Branch("charge_SH",&charge_SH,"charge_SH/F");
  otree->Branch("QDCcali_constant",&constant,"constant/F");
  otree->Branch("QDCcali_slope",&slope,"slope/F");
  otree->Branch("QDCcali_qmax",&qmax,"qmax/F");
  otree->Branch("delta_coarse",&delta_coarse,"delta_coarse/F");
  otree->Branch("pos_phi",&pos_phi,"pos_phi/F");
  otree->Branch("count",&count,"count/I");
  otree->Branch("timestamp",&timestamp,"timestamp/I");
  otree->Branch("l1ts_min_tcoarse",&l1ts_min_tcoarse,"l1ts_min_tcoarse/I");
  otree->Branch("tcoarse",&tcoarse,"tcoarse/F");
  otree->Branch("ecoarse",&ecoarse,"ecoarse/F");
  otree->Branch("tfine_uncal",&tfine_uncal,"tfine_uncal/F");
  otree->Branch("efine_uncal",&efine_uncal,"efine_uncal/F");
  otree->Branch("ttrigg",&ttrigg,"ttrigg/F");
  otree->Branch("trigg_tcoarse",&trigg_tcoarse,"trigg_tcoarse/F");
  otree->Branch("lasttigerframenum",&lasttigerframenum,"lasttigerframenum/I");
  otree->Branch("chip",&chip,"chip/I");
  otree->Branch("FEB_label",&FEB_label,"FEB_label/I");
  otree->Branch("tfine",&tfine,"tfine/F");
  otree->Branch("efine",&efine,"efine/F");
  otree->Branch("tac",&tac,"tac/I");
  otree->Branch("trigg_flag",&trigg_flag,"trigg_flag/I");
  otree->Branch("time",&time,"time/F");
  
  
  max_count = 0;
  for (int i = 0; i < datatree->GetEntries(); i++) {
    datatree->GetEntry(i);
    if(dcount>max_count) max_count = dcount;
  }
  
  for (int i = 0; i < datatree->GetEntries(); i++) {
    datatree->GetEntry(i);
    runNo = drunNo;
    layer = dlayer;
    if(layer==1){
      radius = 90.223;
    }
    else if(layer==2){
      radius = 129.8;
    }
    
    channel = dchannel;
    gemroc = dgemroc;
    FEB = dFEB;
    count =  dcount;
    timestamp =  dtimestamp;
    charge_SH_uncal = dcharge_SH;
    delta_coarse = ddelta_coarse;
    tcoarse = dtcoarse;
    ecoarse = decoarse;
    tfine_uncal = dtfine;
    efine_uncal = define;
    tac = dtac;
    l1ts_min_tcoarse = dl1ts_min_tcoarse;
    lasttigerframenum = dlasttigerframenum;
    constant = cons[dgemroc][dFEB][dchannel];
    slope = mslope[dgemroc][dFEB][dchannel];
    qmax = mqmax[dgemroc][dFEB][dchannel];
    strip_x = mx[dgemroc][dFEB][dchannel];
    strip_v = mv[dgemroc][dFEB][dchannel];
    pos_phi = mphi[dgemroc][dFEB][dchannel];
    chip = mchip_id[dgemroc][dFEB][dchannel];
    FEB_label = mFEB_label_id[dgemroc][dFEB][dchannel];
    if(charge_SH_uncal >= 1008) charge_SH = ((-1*constant)-(1024-charge_SH_uncal))/slope; //-1*(constant/slope);
    else charge_SH = (-1 * constant + charge_SH_uncal) / slope;
    
    if(tfine_uncal<150||tfine_uncal>600) continue;
    //the units for tfine and efine with correction are in ns
    tfine = (tfine_uncal - TDCcons_Tmin[dgemroc][dFEB][dchannel][tac]) * TDCcons_Tbin[dgemroc][dFEB][dchannel][tac];
    efine = (efine_uncal - TDCcons_Emin[dgemroc][dFEB][dchannel][tac]) * TDCcons_Ebin[dgemroc][dFEB][dchannel][tac];
    
    //globle time in ns
    time = (timestamp * pow(2,15) + tcoarse)*6.25 - tfine;
    
    trigg_flag = 0;
    if(dchannel == trigg_channel){ // {&& dFEB == trigg_FEB ){//&& dgemroc == trigg_gemroc) {
      trigg_flag = 1;
      ttrigg = l1ts_min_tcoarse;
      trigg_tcoarse = tcoarse;
      //charge_SH = 78;
    }
    if(!save_TP && dchannel == trigg_channel) continue;
    otree->Fill();
  }
  
  ofile->Write();
  ofile->Close();
  
}

