import ROOT
import copy
import sys
from math import *
from array import array
from operator import itemgetter, attrgetter

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()

from DataFormats.FWLite import Events, Handle

ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.gROOT.ProcessLine("namespace edm {typedef edm::Wrapper<vector<float> > Wrapper<vector<float,allocator<float> > >; }");
ROOT.gROOT.ProcessLine("namespace edm {typedef edm::Wrapper<vector<double> > Wrapper<vector<double,allocator<double> > >; }");

############################################################
############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')

#parser.add_option('--makeCards', action='store_true', dest='makeCards', default=False,
#                  help='no X11 windows')
#
#parser.add_option('--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=50)
parser.add_option('--nBx',action="store",type="int",dest="nBx",default=50)
parser.add_option('--nEv',action="store",type="int",dest="nEv",default=500)
parser.add_option('--index',action="store",type="int",dest="index",default=1)



(options, args) = parser.parse_args()
############################################################
############################################################

def PlotCaloJet(histograms, calojets):
  for calojet in calojets:
    if calojet.pt() > 15:
      histograms["h_calojet_pt"].Fill( calojet.pt() );
      histograms["h_calojet_eta"].Fill( calojet.eta() );
  

def PlotPFJet(histograms, pfjets):
  for pfjet in pfjets:
    if pfjet.pt() > 30:
      histograms["h_pfjet_pt"].Fill( pfjet.pt() );
      histograms["h_pfjet_eta"].Fill( pfjet.eta() );

def PlotCHSJet(histograms, chsjets):
  for chsjet in chsjets:
    if chsjet.pt() > 30:
      histograms["h_chsjet_pt"].Fill( chsjet.pt() );
      histograms["h_chsjet_eta"].Fill( chsjet.eta() );

def PlotGenJet(histograms, genjets):
  for genjet in genjets:
    if genjet.pt() > 30:
      histograms["h_genjet_pt"].Fill( genjet.pt() );
      histograms["h_genjet_eta"].Fill( genjet.eta() );
  
def PlotHbHeRechit(histograms, hbherechits):
  histograms["h_n_rechits"].Fill( hbherechits.size() );
  for hbherechit in hbherechits:
    histograms["h_rechits_e"].Fill( hbherechit.energy() );
    histograms["h_rechits_ieta"].Fill( hbherechit.id().ieta() );
    histograms["h_rechits_time"].Fill( hbherechit.time() );
    histograms["h_rechits_timeVe"].Fill( hbherechit.time(), hbherechit.energy() );
    histograms["h_rechits_ietaVe"].Fill( hbherechit.id().ieta(), hbherechit.energy() );
    histograms["h_rechits_ietaVdepth"].Fill( hbherechit.id().ieta(), hbherechit.id().depth() );

def PlotEBRechit(histograms, EBrechits):
  histograms["h_n_EBrechits"].Fill( EBrechits.size() );
  for EBrechit in EBrechits:
    histograms["h_EBrechits_e"].Fill( EBrechit.energy() );
    #histograms["h_EBrechits_ieta"].Fill( EBrechit.id().eta() );
    histograms["h_EBrechits_time"].Fill( EBrechit.time() );

def PlotEERechit(histograms, EErechits):
  histograms["h_n_EErechits"].Fill( EErechits.size() );
  for EErechit in EErechits:
    histograms["h_EErechits_e"].Fill( EErechit.energy() );
    #histograms["h_EErechits_ieta"].Fill( EErechit.id().ieta() );
    histograms["h_EErechits_time"].Fill( EErechit.time() );

def PlotEKRechit(histograms, EKrechits):
  histograms["h_n_EKrechits"].Fill( EKrechits.size() );
  for EKrechit in EKrechits:
    histograms["h_EKrechits_e"].Fill( EKrechit.energy() );
    #histograms["h_EKrechits_ieta"].Fill( EKrechit.id().ieta() );
    histograms["h_EKrechits_time"].Fill( EKrechit.time() );

def PlotHFRechit(histograms, hfrechits):
  histograms["h_n_hfrechits"].Fill( hfrechits.size() );
  for hfrechit in hfrechits:
    histograms["h_hfrechits_e"].Fill( hfrechit.energy() );
    histograms["h_hfrechits_ieta"].Fill( hfrechit.id().ieta() );
    histograms["h_hfrechits_time"].Fill( hfrechit.time() );
    histograms["h_hfrechits_timeVe"].Fill( hfrechit.time(), hfrechit.energy() );
    histograms["h_hfrechits_ietaVe"].Fill( hfrechit.id().ieta(), hfrechit.energy() );
    histograms["h_hfrechits_ietaVdepth"].Fill( hfrechit.id().ieta(), hfrechit.id().depth() );



def PlotCaloTowers(histograms, calotowers):
  histograms["h_n_calotowers"].Fill( calotowers.size() );
  for calotower in calotowers:
    histograms["h_calotowers_pt"].Fill( calotower.pt() );
    histograms["h_calotowers_eta"].Fill( calotower.eta() );
    histograms["h_calotowers_etaVe"].Fill( calotower.eta(), calotower.energy() );
    histograms["h_calotowers_ietaViphi"].Fill( calotower.ieta(), calotower.iphi() );                                
    histograms["h_calotowers_ietaVhadEm"].Fill( calotower.ieta(), calotower.hadEnergy() ); 

    if calotower.eta() < 5 and calotower.eta() > -5:
      curBin = histograms["h_calotowers_etaVsume"].FindBin(calotower.eta());
      histograms["h_calotowers_etaVsume"].SetBinContent( curBin, histograms["h_calotowers_etaVsume"].GetBinContent(curBin) + calotower.energy() );
      histograms["h_calotowers_etaVsumEt"].SetBinContent( curBin, histograms["h_calotowers_etaVsumEt"].GetBinContent(curBin) + calotower.et() );
      histograms["h_calotowers_etaVsumEMEt"].SetBinContent( curBin, histograms["h_calotowers_etaVsumEMEt"].GetBinContent(curBin) + calotower.emEt() );
      histograms["h_calotowers_etaVsumHCEt"].SetBinContent( curBin, histograms["h_calotowers_etaVsumHCEt"].GetBinContent(curBin) + calotower.hadEt() );


def PlotPFCands(histograms, PFCands):
    histograms["h_n_PFcands"].Fill( PFCands.size() );
   #enum ParticleType {
      #X=0, // undefined
      #h, // charged hadron
      #e, // electron
      #mu, // muon
      #gamma, // photon
      #h0, // neutral hadron
      #h_HF, // HF tower identified as a hadron
      #egamma_HF // HF tower identified as an EM particle
    #};

    n_PFCand_h = 0;
    n_PFCand_e = 0;
    n_PFCand_m = 0;
    n_PFCand_g = 0;
    n_PFCand_n = 0;
    for PFCand in PFCands:
      #print PFCand.particleId()
      histograms["h_PFcands_pt"].Fill(PFCand.pt())
      histograms["h_PFcands_eta"].Fill(PFCand.eta())
      histograms["h_PFcands_phi"].Fill(PFCand.phi())
      if PFCand.particleId() == 1:
        n_PFCand_h += 1;
        histograms["h_PFcand_h_pt"].Fill(PFCand.pt())
        histograms["h_PFcand_h_eta"].Fill(PFCand.eta())
        histograms["h_PFcand_h_phi"].Fill(PFCand.phi())
      if PFCand.particleId() == 2:
        n_PFCand_e += 1;
        histograms["h_PFcand_e_pt"].Fill(PFCand.pt())
        histograms["h_PFcand_e_eta"].Fill(PFCand.eta())
        histograms["h_PFcand_e_phi"].Fill(PFCand.phi())
      if PFCand.particleId() == 3:
        n_PFCand_m += 1;
        histograms["h_PFcand_m_pt"].Fill(PFCand.pt())
        histograms["h_PFcand_m_eta"].Fill(PFCand.eta())
        histograms["h_PFcand_m_phi"].Fill(PFCand.phi())
      if PFCand.particleId() == 4:
        n_PFCand_g += 1;
        histograms["h_PFcand_g_pt"].Fill(PFCand.pt())
        histograms["h_PFcand_g_eta"].Fill(PFCand.eta())
        histograms["h_PFcand_g_phi"].Fill(PFCand.phi())
      if PFCand.particleId() == 5:
        n_PFCand_n += 1;
        histograms["h_PFcand_n_pt"].Fill(PFCand.pt())
        histograms["h_PFcand_n_eta"].Fill(PFCand.eta())
        histograms["h_PFcand_n_phi"].Fill(PFCand.phi())
  
    histograms["h_n_PFcand_h"].Fill(n_PFCand_h);
    histograms["h_n_PFcand_e"].Fill(n_PFCand_e);
    histograms["h_n_PFcand_m"].Fill(n_PFCand_m);
    histograms["h_n_PFcand_g"].Fill(n_PFCand_g);
    histograms["h_n_PFcand_n"].Fill(n_PFCand_n);


def PlotNeuCands(histograms, PFCands, Vertexes):
    histograms["h_n_offlinePrimaryVertices"].Fill(Vertexes.size());
    n_PFCand_n = 0;
    for PFCand in PFCands:
      if PFCand.particleId() == 5:
        n_PFCand_n += 1;
    histograms["h_n_PFcand_n_Vertex"].Fill(Vertexes.size(), n_PFCand_n);
    histograms["h_n_PFNeu_Vertex"].Fill(Vertexes.size(), n_PFCand_n);
    


def PlotTracks(histograms, Tracks):
  for Track in Tracks:
    histograms["h_Track_OuterEta"].Fill(Track.outerEta())

  



if __name__ == '__main__':

    nEvents = options.nEv;
    nPU = options.nPU; 
    nBx = options.nBx;
    index = options.index;
    ofilename = "hists/fileOut_PU"+str(nPU)+"bx"+str(nBx)+"_"+str(options.index)+".root";

    files = []
    dirEos = "file:/eos/uscms/store/user/ntran/upgradeJet/files/";
    if nPU == 19 and nBx == 0: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/reco_calotow_phase1_0PU.root")
    elif nPU == 19 and nBx == 3: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/reco_3d_phase1_0PU.root")
    elif nPU == 23 and nBx == 0: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/reco_calotow_phase2shcal_0PU.root")
    elif nPU == 23 and nBx == 3: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/reco_3d_phase2shcal_0PU.root")
    elif nPU == 19 and nBx == 10: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/testcalotow_phase1_PU140bx25.root")
    elif nPU == 19 and nBx == 13: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Mark_SLHC16/PhaseI_3D.root")
        #files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/test3dcluster_phase1_PU140bx25.root")
    elif nPU == 23 and nBx == 10: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/testcalotow_phase2_PU140bx25.root")
    elif nPU == 23 and nBx == 13: 
        files.append("/data/nbay04/c/benwu/JetMET_TP/Mark_SLHC16/PhaseII_3D.root")
        #files.append("/data/nbay04/c/benwu/JetMET_TP/Gaelle_SLHC15/test3dcluster_phase2_PU140bx25.root")
    else:
        print "invalid file settings...";
        sys.exit();
        
    events = Events( files );

    # loop over events
    count = 0
    ntotal = events.size()
    #print "Nevents = "+str(ntotal)
    print "Start looping"

    photonHandle = Handle( "vector<reco::Photon>" );
    photonLabel = "photons";

    calojetHandle = Handle( "vector<reco::CaloJet>" );
    calojetLabel = "ak5CaloJets";

    pfjetHandle = Handle( "vector<reco::PFJet>" );
    pfjetLabel = "ak5PFJets";

    chsjetHandle = Handle( "vector<reco::PFJet>" );
    chsjetLabel = "ak5PFJetsCHS";

    genjetHandle = Handle( "vector<reco::GenJet>" );
    genjetLabel = "ak5GenJets";
    
    rechitHandle_EB= Handle( "edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >" );
    rechitHandle_EE= Handle( "edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >" );
    rechitHandle_EK= Handle( "edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >" );
    ecalrechitLabel = "ecalRecHit";
    rechitLabel_EB = "EcalRecHitsEB";
    rechitLabel_EE = "EcalRecHitsEE";
    rechitLabel_EK = "EcalRecHitsEK";

    simhitHandle_EE = Handle( "vector<PCaloHit>" );
    simhitHandle_EK = Handle( "vector<PCaloHit>" );
    simhitLabel = "g4SimHits";
    simhitLabel_EE = "EcalHitsEE";
    simhitLabel_EK = "EcalHitsEK";

    rechitHandle_hbhe = Handle( "edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >" );
    rechitHandle_hf = Handle( "edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >" );    
    hcalrechitLabel = "reducedHcalRecHits";
    rechitLabel2_hbhe = "hbheUpgradeReco";
    rechitLabel2_hf = "hfUpgradeReco";    


    calotowerHandle = Handle( "edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower> >" );
    calotowerLabel = "towerMaker";

    ## some pileup information
    pileupHandle = Handle( "vector<PileupSummaryInfo>" );
    pileupLabel = "addPileupInfo";

    VertexHandle = Handle( "vector<reco::Vertex>" );
    VertexLabel = "offlinePrimaryVertices";

    PFcandHandle = Handle("vector<reco::PFCandidate>");
    PFcandLabel = "particleFlow";
    
    TrackHandle = Handle("vector<reco::Track>");
    TrackLabel = "generalTracks";
    
    ##output
    fileOut = ROOT.TFile(ofilename,"RECREATE");
    histograms = {};
    histograms["h_calojet_pt"]               = ROOT.TH1F("h_calojet_pt"               , "; pT; N"                       , 40  , 0   , 100) ;                   # 0th
    histograms["h_calojet_eta"]              = ROOT.TH1F("h_calojet_eta"              , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_pfjet_pt"]                 = ROOT.TH1F("h_pfjet_pt"                 , "; pT; N"                       , 40  , 0   , 100) ;
    histograms["h_pfjet_eta"]                = ROOT.TH1F("h_pfjet_eta"                , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_chsjet_pt"]                = ROOT.TH1F("h_chsjet_pt"                , "; pT; N"                       , 40  , 0   , 100) ;
    histograms["h_chsjet_eta"]               = ROOT.TH1F("h_chsjet_eta"               , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_genjet_pt"]                = ROOT.TH1F("h_genjet_pt"                , "; pT; N"                       , 40  , 0   , 100) ;
    histograms["h_genjet_eta"]               = ROOT.TH1F("h_genjet_eta"               , "; eta; N"                      , 100 , -5  , 5) ;                  # 5th

    histograms["h_n_rechits"]                = ROOT.TH1F("h_n_rechits"                , "; N hbhe rechits; N"           , 100 , 0   , 1000) ;
    histograms["h_rechits_e"]                = ROOT.TH1F("h_rechits_e"                , "; energy; N"                   , 40  , 0   , 20);
    histograms["h_rechits_ieta"]             = ROOT.TH1F("h_rechits_ieta"             , "; ieta; N"                     , 100 , -50 , 50) ;
    histograms["h_rechits_time"]             = ROOT.TH1F("h_rechits_time"             , "; time; N"                     , 50  , -50 , 50) ;

    histograms["h_n_EBrechits"]              = ROOT.TH1F("h_n_EBrechits"              , "; N EB rechits; N"             , 100 , 0   , 1000) ;
    histograms["h_EBrechits_e"]              = ROOT.TH1F("h_EBrechits_e"              , "; EB energy; N"                , 40  , 0   , 20);
    histograms["h_EBrechits_ieta"]           = ROOT.TH1F("h_EBrechits_ieta"           , "; EB ieta; N"                  , 100 , -50 , 50) ;
    histograms["h_EBrechits_time"]           = ROOT.TH1F("h_EBrechits_time"           , "; EBtime; N"                   , 50  , -50 , 50) ;

    histograms["h_n_EErechits"]              = ROOT.TH1F("h_n_EErechits"              , "; N EE rechits; N"             , 100 , 0   , 1000) ;
    histograms["h_EErechits_e"]              = ROOT.TH1F("h_EErechits_e"              , "; EE energy; N"                , 40  , 0   , 20);
    histograms["h_EErechits_ieta"]           = ROOT.TH1F("h_EErechits_ieta"           , "; EE ieta; N"                  , 100 , -50 , 50) ;
    histograms["h_EErechits_time"]           = ROOT.TH1F("h_EErechits_time"           , "; EEtime; N"                   , 50  , -50 , 50) ;

    histograms["h_n_EKrechits"]              = ROOT.TH1F("h_n_EKrechits"              , "; N EK rechits; N"             , 100 , 0   , 1000) ;
    histograms["h_EKrechits_e"]              = ROOT.TH1F("h_EKrechits_e"              , "; EK energy; N"                , 40  , 0   , 20);
    histograms["h_EKrechits_ieta"]           = ROOT.TH1F("h_EKrechits_ieta"           , "; EK ieta; N"                  , 100 , -50 , 50) ;
    histograms["h_EKrechits_time"]           = ROOT.TH1F("h_EKrechits_time"           , "; EKtime; N"                   , 50  , -50 , 50) ;

    histograms["h_n_hfrechits"]              = ROOT.TH1F("h_n_hfrechits"              , "; N hf rechits; N"             , 100 , 0   , 1000) ;      #10th
    histograms["h_hfrechits_e"]              = ROOT.TH1F("h_hfrechits_e"              , "; energy; N"                   , 50  , 0   , 100) ;
    histograms["h_hfrechits_ieta"]           = ROOT.TH1F("h_hfrechits_ieta"           , "; ieta; N"                     , 100 , -50 , 50) ;
    histograms["h_hfrechits_time"]           = ROOT.TH1F("h_hfrechits_time"           , "; time; N"                     , 50  , -50 , 50) ;

    histograms["h_n_calotowers"]             = ROOT.TH1F("h_n_calotowers" , "; N calotowers; N" , 300 , 3000   , 6000) ;
    histograms["h_calotowers_pt"]            = ROOT.TH1F("h_calotowers_pt"            , "; pT; N"                       , 50  , 0   , 100) ;                 #15th
    histograms["h_calotowers_eta"]           = ROOT.TH1F("h_calotowers_eta" , "; CaloTower #eta; N"                      , 100 , -5  , 5) ;

    histograms["h_calotowers_etaVe"]         = ROOT.TH2F("h_calotowers_etaVe"         , "; eta; e"                      , 100 , -5  , 5                              , 50 , 0 , 200) ;
    histograms["h_rechits_timeVe"]           = ROOT.TH2F("h_rechits_timeVe"           , "; time; e"                     , 50  , -50 , 50                             , 50 , 0 , 100) ;
    histograms["h_rechits_ietaVe"]           = ROOT.TH2F("h_rechits_ietaVe"           , "; ieta; e"                     , 100 , -50 , 50                             , 50 , 0 , 100) ;
    histograms["h_hfrechits_timeVe"]         = ROOT.TH2F("h_hfrechits_timeVe"         , "; time; e"                     , 50  , -50 , 50                             , 50 , 0 , 100) ;    #20th
    histograms["h_hfrechits_ietaVe"]         = ROOT.TH2F("h_hfrechits_ietaVe"         , "; ieta; e"                     , 100 , -50 , 50                             , 50 , 0 , 100) ;
    histograms["h_calotowers_ietaViphi"]     = ROOT.TH2F("h_calotowers_ietaViphi"     , "; ieta; iphi"                  , 100 , -50 , 50                             , 80 , 0 , 80) ;
    histograms["h_calotowers_ietaVhadEm"]    = ROOT.TH2F("h_calotowers_ietaVhadEm"    , "; ieta; hadronic energy"       , 100 , -50 , 50                             , 50 , 0 , 200) ;
    histograms["h_rechits_ietaVdepth"]       = ROOT.TH2F("h_rechits_ietaVdepth"       , "; ieta; depth"                 , 100 , -50 , 50                             , 5  , 0 , 5) ;
    histograms["h_hfrechits_ietaVdepth"]     = ROOT.TH2F("h_hfrechits_ietaVdepth"     , "; ieta; depth"                 , 100 , -50 , 50                             , 5  , 0 , 5) ; #25th
    histograms["h_calotowers_etaVsume"]      = ROOT.TH1F("h_calotowers_etaVsume"      , "; CaloTower #eta; sum e"                  , 100 , -5  , 5) ;
    histograms["h_calotowers_etaVsumEt"]     = ROOT.TH1F("h_calotowers_etaVsumEt"     , "; CaloTower #eta; sum et"                 , 100 , -5  , 5) ;
    histograms["h_calotowers_etaVsumEMEt"]   = ROOT.TH1F("h_calotowers_etaVsumEMEt"   , "; CaloTower #eta; sum EM et"              , 100 , -5  , 5) ;
    histograms["h_calotowers_etaVsumHCEt"]   = ROOT.TH1F("h_calotowers_etaVsumHCEt"   , "; CaloTower #eta; sum HC et"              , 100 , -5  , 5) ;
    
#============================================================================#
#------------------------------     PFCand     ------------------------------#
#============================================================================#
    histograms["h_n_PFcands"]                = ROOT.TH1F("h_n_PFcands" , "; N PFcands; N" , 500 , 3000  , 8000) ;
    histograms["h_PFcands_pt"]               = ROOT.TH1F("h_PFcands_pt"               , "; pT; N"                       , 50  , 0   , 100) ;                 #15th
    histograms["h_PFcands_eta"]              = ROOT.TH1F("h_PFcands_eta"              , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_PFcands_phi"]              = ROOT.TH1F("h_PFcands_phi"              , "; phi; N"                      , 100 , -5  , 5) ;

    histograms["h_n_PFcand_h"]               = ROOT.TH1F("h_n_PFcand_h"               , "; N PFcand_h; N"               , 100 , 0   , 5000) ;
    histograms["h_PFcand_h_pt"]              = ROOT.TH1F("h_PFcand_h_pt"              , "; pT; N"                       , 50  , 0   , 100) ;                 #15th
    histograms["h_PFcand_h_eta"]             = ROOT.TH1F("h_PFcand_h_eta"             , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_PFcand_h_phi"]             = ROOT.TH1F("h_PFcand_h_phi"             , "; phi; N"                      , 100 , -5  , 5) ;

    histograms["h_n_PFcand_e"]               = ROOT.TH1F("h_n_PFcand_e"               , "; N PFcand_e; N"               , 100 , 0   , 5000) ;
    histograms["h_PFcand_e_pt"]              = ROOT.TH1F("h_PFcand_e_pt"              , "; pT; N"                       , 50  , 0   , 100) ;                 #15th
    histograms["h_PFcand_e_eta"]             = ROOT.TH1F("h_PFcand_e_eta"             , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_PFcand_e_phi"]             = ROOT.TH1F("h_PFcand_e_phi"             , "; phi; N"                      , 100 , -5  , 5) ;

    histograms["h_n_PFcand_m"]               = ROOT.TH1F("h_n_PFcand_m"               , "; N PFcand_m; N"               , 100 , 0   , 5000) ;
    histograms["h_PFcand_m_pt"]              = ROOT.TH1F("h_PFcand_m_pt"              , "; pT; N"                       , 50  , 0   , 100) ;                 #15th
    histograms["h_PFcand_m_eta"]             = ROOT.TH1F("h_PFcand_m_eta"             , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_PFcand_m_phi"]             = ROOT.TH1F("h_PFcand_m_phi"             , "; phi; N"                      , 100 , -5  , 5) ;
    histograms["h_n_PFcand_g"]               = ROOT.TH1F("h_n_PFcand_g" , "; N PFcand_g; N" , 200 , 0   , 2000) ;
    histograms["h_PFcand_g_pt"]              = ROOT.TH1F("h_PFcand_g_pt"              , "; pT; N"                       , 50  , 0   , 100) ;                 #15th
    histograms["h_PFcand_g_eta"]             = ROOT.TH1F("h_PFcand_g_eta"             , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_PFcand_g_phi"]             = ROOT.TH1F("h_PFcand_g_phi"             , "; phi; N"                      , 100 , -5  , 5) ;

    histograms["h_n_PFcand_n"]               = ROOT.TH1F("h_n_PFcand_n"               , "; N PFcand_n; N"               , 100 , 0   , 1000) ;
    histograms["h_PFcand_n_pt"]              = ROOT.TH1F("h_PFcand_n_pt"              , "; pT; N"                       , 50  , 0   , 100) ;                 #15th
    histograms["h_PFcand_n_eta"]             = ROOT.TH1F("h_PFcand_n_eta"             , "; eta; N"                      , 100 , -5  , 5) ;
    histograms["h_PFcand_n_phi"]             = ROOT.TH1F("h_PFcand_n_phi"             , "; phi; N"                      , 100 , -5  , 5) ;

    histograms["h_n_offlinePrimaryVertices"] = ROOT.TH1F("h_n_offlinePrimaryVertices" , "; N offlinePrimaryVertices; N" , 100 , 60  , 160) ;
    
#============================================================================#
#-------------------------------     Track     ------------------------------#
#============================================================================#
    histograms["h_Track_OuterEta"]             = ROOT.TH1F("h_Track_OuterEta" , "; Outer eta; N"                      , 100 , -5  , 5) ;

    if nBx >= 10:
      print "plotting"
      histograms["h_n_PFcand_n_Vertex"]          = ROOT.TProfile("h_n_PFcand_n_Vertex"           , "; N vertex; N PFcand_n; N"  , 60 , 60 , 120, 's' ) ;
    else:
      histograms["h_n_PFcand_n_Vertex"]          = ROOT.TProfile("h_n_PFcand_n_Vertex"           , "; N vertex; N PFcand_n; N"  , 20 , 0   , 20, 's' ) ;
    histograms["h_n_PFNeu_Vertex"]          = ROOT.TProfile("h_n_PFNeu_Vertex" , "; N vertex; N PFcand_neutral; N" , 140 , 0   , 140, 's' ) ;
    for i in range(1,101):
        histograms["h_calotowers_etaVsume"].SetBinContent(i,0.);
        histograms["h_calotowers_etaVsumEt"].SetBinContent(i,0.);
        histograms["h_calotowers_etaVsumEMEt"].SetBinContent(i,0.);
        histograms["h_calotowers_etaVsumHCEt"].SetBinContent(i,0.);
    
    #for event in events:
    ctr = 0;
    for event in events:

        if ctr == nEvents: break;
        #if ctr % 500 == 0: 
          #print "ctr = ", ctr;
        print "ctr = ", ctr;
        ctr = ctr + 1;

        if ctr < 10:
            event.getByLabel( pileupLabel, pileupHandle );
            pileupInfos = pileupHandle.product();
            print "pileupInfos.size() = ",pileupInfos.size();
            #for pv in pileupInfos:
                #print "PU INFO: bx = ",pv.getBunchCrossing(),", nPU = ",pv.getPU_NumInteractions();      
        
        # photons
        #print "photons..."
        event.getByLabel( photonLabel, photonHandle );
        photons = photonHandle.product();

        # calo jets
        #print "calo jets..."    
        event.getByLabel( calojetLabel, calojetHandle );
        calojets = calojetHandle.product();

        # pf jets
        #print "pf jets..."    
        event.getByLabel( pfjetLabel, pfjetHandle );
        pfjets = pfjetHandle.product();

        event.getByLabel( chsjetLabel, chsjetHandle );
        chsjets = chsjetHandle.product();

        # gen jets
        event.getByLabel( genjetLabel, genjetHandle );
        genjets = genjetHandle.product();

        # calo towers
        #print "calo towers..."             
        event.getByLabel( calotowerLabel, calotowerHandle );
        calotowers = calotowerHandle.product();
        
        ## ecal
        event.getByLabel(ecalrechitLabel, rechitLabel_EB, rechitHandle_EB);
        EBrechits = rechitHandle_EB.product();        

        event.getByLabel(ecalrechitLabel, rechitLabel_EE, rechitHandle_EE);
        EErechits = rechitHandle_EE.product();        

        event.getByLabel(ecalrechitLabel, rechitLabel_EK, rechitHandle_EK);
        EKrechits = rechitHandle_EK.product();        

        event.getByLabel(simhitLabel, simhitLabel_EE, simhitHandle_EE);
        EEsimhits = simhitHandle_EE.product();        

        event.getByLabel(simhitLabel, simhitLabel_EK, simhitHandle_EK);
        EKsimhits = simhitHandle_EK.product();        

        # hbhe rechits
        #print "rechits..."    
        event.getByLabel(rechitLabel2_hbhe, "", rechitHandle_hbhe );
        hbherechits = rechitHandle_hbhe.product();        

        event.getByLabel(rechitLabel2_hf, "", rechitHandle_hf );        
        hfrechits = rechitHandle_hf.product();   
             
        event.getByLabel(PFcandLabel, PFcandHandle );        
        PFCands = PFcandHandle.product();   

        event.getByLabel(TrackLabel, TrackHandle );        
        Tracks = TrackHandle.product();   

        event.getByLabel(VertexLabel, VertexHandle );        
        Vertexes = VertexHandle.product();   
#        print "-----------------------------------"
#        print "photon size = ",photons.size();
#        print "calojets size = ",calojets.size(); 
#        print "pfjets size = ",pfjets.size();           
#        print "hbherechits size = ",hbherechits.size();           
#        print "calo tower size = ",calotowers.size();           

                    
        PlotCaloJet(histograms, calojets)
        PlotPFJet(histograms, pfjets)
        PlotCHSJet(histograms, chsjets)
        PlotGenJet(histograms, genjets)

#        print "n rechits = ", hbherechits.size();
        PlotEBRechit(histograms, EBrechits)
        PlotEERechit(histograms, EErechits)
        PlotEKRechit(histograms, EKrechits)
        PlotHbHeRechit(histograms, hbherechits)
        PlotHFRechit(histograms, hfrechits)

        PlotCaloTowers(histograms, calotowers)
        PlotPFCands(histograms, PFCands)
        PlotNeuCands(histograms, PFCands, Vertexes)
                
        PlotTracks(histograms, Tracks)

#        print "n calotowers = ", calotowers.size();    
                        
    fileOut.cd();
    for i in histograms.keys():
        #histograms[i].Scale( 1./histograms[i].Integral() );
        histograms[i].Write();
    fileOut.Close();




