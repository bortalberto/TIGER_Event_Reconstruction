#include "ana.h"

const int NRoc = 14;
const int NChannel = 64;
const int NFeb = 60;

void ana(int run, int subrun){
  std::string iname=ANADIR;
  iname=iname+std::to_string(run)+"/Sub_RUN_dec_"+std::to_string(subrun)+".root";
  std::string oname=ANADIR;
  oname=oname+std::to_string(run)+"/Sub_RUN_ana_"+std::to_string(subrun)+".root";

  std::ifstream inStream(iname, std::ios::binary);
  if (!inStream) {
    std::cerr << "File " << iname << " not found" << std::endl;
    return;
  }


  int trigg_channel=20;
  if(run>=118) trigg_channel=62;
  if(run==337) trigg_channel=-1;
  TString map_file = "mapping_IHEP.root";
  TString qdc_file = "QDCcalib.root";
  TString tdc_file = "TDCcalib.root";
  if(run>=250 && run<264) {
    map_file ="mapping_IHEP_planar_ROC3.root";
    qdc_file = "QDCcalib_planar_ROC3.root";
    tdc_file = "TDCcalib_planar_ROC3.root";
  }
  if(run>=264) {
    map_file="mapping_IHEP_L2_2planari.root";
    qdc_file = "QDCcalib_L2_2planari.root";
    tdc_file = "TDCcalib_L2_2planari.root";
  }
  if(run>=274) {
    map_file="mapping_IHEP_L2_2planari_bis.root";
    qdc_file = "QDCcalib_L2_2planari_bis.root";
    tdc_file = "TDCcalib_L2_2planari_bis.root";
  }
  if(run>=277) {
    map_file="mapping_IHEP_L2_2planari_tris.root";
    qdc_file = "QDCcalib_L2_2planari_tris.root";
    tdc_file = "TDCcalib_L2_2planari_tris.root";
  }
  if(run>=281){
    map_file="mapping_IHEP_L2_2planari_quad.root";
    qdc_file = "QDCcalib_L2_2planari_quad.root";
    tdc_file = "TDCcalib_L2_2planari_quad.root";
  }
  if(run==273 || run==267){ // TEMPORARY - 16/11/2019 - GM - TO BE REMOVED
    map_file="mapping_IHEP_L2_2planari_quad.root";
    qdc_file="QDCcalib_L2_2planari_quad.root";
    tdc_file="TDCcalib_l2_2planari_quad.root";
  }
  if(run>=286){
    map_file=  "mapping_IHEP_L2_2planari_penta.root";
    qdc_file =     "QDCcalib_L2_2planari_penta.root";
    tdc_file =     "TDCcalib_L2_2planari_penta.root";
  }
  if(run>=317){
    map_file= "mapping_IHEP_L1_L2_2planar.root";
    qdc_file=     "QDCcalib_L1_L2_2planar.root";
    tdc_file=     "TDCcalib_L1_L2_2planar.root";
  }

  TString tot_file1 = "delta_vth.root";
  TString tot_file2 = "ToT_calib.root";
  bool save_TP = true;
  //int trigg_FEB=4;
  //int trigg_gemroc=4;
  auto mapfile = new TFile(map_file);
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
  
  int mx[NRoc][8][NChannel];
  float mphi[NRoc][8][NChannel];
  int mv[NRoc][8][NChannel];
  int mchip_id[NRoc][8][NChannel];
  int mFEB_label_id[NRoc][8][NChannel];
  memset(mx, -1, sizeof(mx));
  memset(mv, -1, sizeof(mv));
  memset(mphi, -1, sizeof(mphi));
  memset(mchip_id, -1, sizeof(mchip_id));
  memset(mFEB_label_id, -1, sizeof(mFEB_label_id));

  cout<<maptree->GetEntries()<<endl;
  int ciao=0;
  for (int i = 0; i < maptree->GetEntries(); i++) {

    maptree->GetEntry(i);
    if(pos_x*pos_v<0 && (mx[gemroc_id][SW_FEB_id][channel_id])==-1 && mv[gemroc_id][SW_FEB_id][channel_id]==-1 && gemroc_id>=11){
      ciao++;
    }
    mx[gemroc_id][SW_FEB_id][channel_id] = pos_x;
    mv[gemroc_id][SW_FEB_id][channel_id] = pos_v;
    mphi[gemroc_id][SW_FEB_id][channel_id] = phi;
    mchip_id[gemroc_id][SW_FEB_id][channel_id] = chip_id;
    mFEB_label_id[gemroc_id][SW_FEB_id][channel_id] = FEB_label_id;
    //if(gemroc_id==0 && channel_id<5) cout<<gemroc_id<<" "<<SW_FEB_id<<" "<<channel_id<<" "<<mx[gemroc_id][SW_FEB_id][channel_id]<<" "<<mv[gemroc_id][SW_FEB_id][channel_id]<<endl;
  }
  cout<<ciao<<"------------"<<endl;
  if(ciao!=maptree->GetEntries()) cout<<"Check the mapping"<<endl;
  //for(int i=0;i<8;i++){
  //  for(int j=0;j<5;j++){
  //    cout<<i<<" "<<j<<" "<<mx[0][i][j]<<" "<<mv[0][i][j]<<endl;
  //  }
  //}
  
  auto consfile = new TFile(qdc_file);
  auto constree = (TTree*)consfile->Get("tree");
  
  int cons_channel_id, cons_gemroc_id, cons_SW_FEB_id;
  float cons_constant, cons_slope, cons_qmax;
  
  constree->SetBranchAddress("channel_id", &cons_channel_id);
  constree->SetBranchAddress("gemroc_id", &cons_gemroc_id);
  constree->SetBranchAddress("SW_FEB_id", &cons_SW_FEB_id);
  constree->SetBranchAddress("constant", &cons_constant);
  constree->SetBranchAddress("slope", &cons_slope);
  constree->SetBranchAddress("qmax", &cons_qmax);
  
  float cons[NRoc][8][NChannel];
  float mslope[NRoc][8][NChannel];
  float mqmax[NRoc][8][NChannel];
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
  
  
  auto TDCconsfile = new TFile(tdc_file);
  auto TDCconstree = (TTree*)TDCconsfile->Get("tree");
  
  int TDCcons_channel_id, TDCcons_gemroc_id, TDCcons_SW_FEB_id;
  float Tac0_Tfine_min, Tac0_Tfine_max, Tac1_Tfine_min, Tac1_Tfine_max, Tac2_Tfine_min, Tac2_Tfine_max, Tac3_Tfine_min, Tac3_Tfine_max, Tac0_Efine_min, Tac0_Efine_max, Tac1_Efine_min, Tac1_Efine_max, Tac2_Efine_min, Tac2_Efine_max, Tac3_Efine_min, Tac3_Efine_max;
  
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
  
  float TDCcons_Tmin[NRoc][8][NChannel][4];
  float TDCcons_Tbin[NRoc][8][NChannel][4];
  float TDCcons_Emin[NRoc][8][NChannel][4];
  float TDCcons_Ebin[NRoc][8][NChannel][4];
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
  
  auto TOTconsfile1 = new TFile(tot_file1);
  auto TOTconstree1 = (TTree*)TOTconsfile1->Get("tree");
  auto TOTconsfile2 = new TFile(tot_file2);
  auto TOTconstree2 = (TTree*)TOTconsfile2->Get("tree");

  int qlayer, qgemroc, qfeb, qchannel, qdelta_vth1_baseline;
  int delta_vth1_baseline[NRoc][8][NChannel];

  TOTconstree1->SetBranchAddress("layer_id", &qlayer);
  TOTconstree1->SetBranchAddress("gemroc_id", &qgemroc);
  TOTconstree1->SetBranchAddress("software_feb_id", &qfeb);
  TOTconstree1->SetBranchAddress("channel_id", &qchannel);
  TOTconstree1->SetBranchAddress("delta_vth1_baseline", &qdelta_vth1_baseline);
  
  for(int i=0; i<TOTconstree1->GetEntries(); i++){
    TOTconstree1->GetEntry(i);
    if(qgemroc<0 || qfeb<0 || qchannel<0) continue;
    delta_vth1_baseline[qgemroc][qfeb][qchannel] = qdelta_vth1_baseline;
  }

  int delta_vth;
  float par_A, par_B, par_C, par_D, par_E;
  float tot_parameters[100][5];

  TOTconstree2->SetBranchAddress("vthr",&delta_vth);
  TOTconstree2->SetBranchAddress("par_a",&par_A);
  TOTconstree2->SetBranchAddress("par_b",&par_B);
  TOTconstree2->SetBranchAddress("par_c",&par_C);
  TOTconstree2->SetBranchAddress("par_d",&par_D);
  TOTconstree2->SetBranchAddress("par_e",&par_E);

  for(int i=0; i<TOTconstree2->GetEntries(); i++){
    TOTconstree2->GetEntry(i);
    if(delta_vth<0) continue;
    tot_parameters[delta_vth][0]=par_A;
    tot_parameters[delta_vth][1]=par_B;
    tot_parameters[delta_vth][2]=par_C;
    tot_parameters[delta_vth][3]=par_D;
    tot_parameters[delta_vth][4]=par_E;
  }


  auto datafile = new TFile(iname.c_str());
  auto datatree = (TTree*)datafile->Get("tree");
  
  int dchannel, dgemroc, dtiger, dcount, dtimestamp, dl1ts_min_tcoarse, dlasttigerframenum, dtac, drunNo, dlayer, dsubRunNo;
  float dcharge_SH, ddelta_coarse, dtcoarse, decoarse, dtfine, define;
  
  datatree->SetBranchAddress("layer", &dlayer);
  datatree->SetBranchAddress("runNo", &drunNo);
  datatree->SetBranchAddress("subRunNo", &dsubRunNo);
  datatree->SetBranchAddress("channel", &dchannel);
  datatree->SetBranchAddress("gemroc", &dgemroc);
  datatree->SetBranchAddress("tiger", &dtiger);
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
  datatree->SetBranchAddress("tac", &dtac);
  
  
  auto ofile = new TFile(oname.c_str(),"RECREATE");
  auto otree = new TTree("tree","tree");
  
  int channel, gemroc, tiger, strip_x, strip_v, count, timestamp, l1ts_min_tcoarse, lasttigerframenum, chip, FEB_label, tac, runNo, layer, max_count, trigg_flag, tcoarse_min_ts, subRunNo;
  float charge_SH, charge_SH_uncal, charge_TOT, charge_TOT_uncal, constant, slope, qmax, delta_coarse, pos_phi, tcoarse, ecoarse, tfine, efine, ttrigg = -99999, trigg_tcoarse = -99999, tfine_uncal, efine_uncal, time, radius; 
  bool saturated = false;
  
  otree->Branch("runNo",&runNo,"runNo/I");
  otree->Branch("subRunNo",&subRunNo,"subRunNo/I");
  otree->Branch("layer",&layer,"layer/I");
  otree->Branch("channel",&channel,"channel/I");
  otree->Branch("gemroc",&gemroc,"gemroc/I");
  otree->Branch("tiger",&tiger,"tiger/I");
  otree->Branch("max_count",&max_count,"max_count/I");
  otree->Branch("strip_x",&strip_x,"strip_x/I");
  otree->Branch("strip_v",&strip_v,"strip_v/I");
  otree->Branch("radius",&radius,"radius/F");
  otree->Branch("charge_SH_uncal",&charge_SH_uncal,"charge_SH_uncal/F");
  otree->Branch("charge_SH",&charge_SH,"charge_SH/F");
  otree->Branch("charge_TOT_uncal",&charge_TOT_uncal,"charge_TOT_uncal/F");
  otree->Branch("charge_TOT",&charge_TOT,"charge_TOT/F");
  otree->Branch("QDCcali_constant",&constant,"constant/F");
  otree->Branch("QDCcali_slope",&slope,"slope/F");
  otree->Branch("QDCcali_qmax",&qmax,"qmax/F");
  otree->Branch("delta_coarse",&delta_coarse,"delta_coarse/F");
  otree->Branch("pos_phi",&pos_phi,"pos_phi/F");
  otree->Branch("count",&count,"count/I");
  otree->Branch("timestamp",&timestamp,"timestamp/I");
  otree->Branch("l1ts_min_tcoarse",&l1ts_min_tcoarse,"l1ts_min_tcoarse/I");
  otree->Branch("tcoarse_min_ts",&tcoarse_min_ts,"tcoarse_min_ts/I");
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
  otree->Branch("saturated",&saturated,"saturated/O");
  
  max_count = 0;
  for (int i = 0; i < datatree->GetEntries(); i++) {
    datatree->GetEntry(i);
    if(dcount>max_count) max_count = dcount;
  }
  
  for (int i = 0; i < datatree->GetEntries(); i++) {
    datatree->GetEntry(i);

    runNo = drunNo;
    subRunNo = dsubRunNo;
    layer = dlayer;
    if(layer==1){
      radius = 90.223;
    }
    else if(layer==2){
      radius = 129.8;
    }
    else if(layer==0){
      radius = 0;
    }
    
    channel = dchannel;
    gemroc = dgemroc;
    tiger = dtiger;
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
    tcoarse_min_ts = -dl1ts_min_tcoarse;
    lasttigerframenum = dlasttigerframenum;
    constant = cons[dgemroc][dtiger][dchannel];
    slope = mslope[dgemroc][dtiger][dchannel];
    qmax = mqmax[dgemroc][dtiger][dchannel];
    strip_x = mx[dgemroc][dtiger][dchannel];
    strip_v = mv[dgemroc][dtiger][dchannel];
    pos_phi = mphi[dgemroc][dtiger][dchannel];
    chip = mchip_id[dgemroc][dtiger][dchannel];
    FEB_label = mFEB_label_id[dgemroc][dtiger][dchannel];
    //SH Calibration
    if(charge_SH_uncal == 1008) saturated = true;
    else saturated = false;
    if(charge_SH_uncal >= 1008) charge_SH = ((-1*constant)-(1024-charge_SH_uncal))/slope; //-1*(constant/slope);
    else charge_SH = (-1 * constant + charge_SH_uncal) / slope;
    
    //if(tfine_uncal<150||tfine_uncal>600) continue;

    //the units for tfine and efine with correction are in ns
    tfine = (tfine_uncal - TDCcons_Tmin[dgemroc][dtiger][dchannel][tac]) * TDCcons_Tbin[dgemroc][dtiger][dchannel][tac];
    if(tfine>10 || tfine<-5) tfine=0;
    efine = (efine_uncal - TDCcons_Emin[dgemroc][dtiger][dchannel][tac]) * TDCcons_Ebin[dgemroc][dtiger][dchannel][tac];
    if(efine>10 || efine<-5) efine=0;

    //TOT Calibration
    TF1 *f_tot = new TF1("f","[4]*exp([0]*x+[1])+[2]+[3]*x",0,500);
    int delta = delta_vth1_baseline[dgemroc][dtiger][dchannel];
    f_tot->SetParameters(tot_parameters[delta][0],tot_parameters[delta][1],tot_parameters[delta][2],tot_parameters[delta][3],tot_parameters[delta][4]);
    charge_TOT_uncal = delta_coarse * 6.25 - efine + tfine;
    charge_TOT = f_tot->GetX(charge_TOT_uncal);
    f_tot->~TF1();

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

