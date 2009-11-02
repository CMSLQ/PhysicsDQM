#ifndef LQeDQM_H
#define LQeDQM_H


/** \class LQeDQM
 *
 *  DQM offline for first generation leptoquarks
 *
 *  $Date: 2009/10/22 $
 *  $Revision: 1.00 $
 *  \author Ellie Twedt, University of Maryland 
 */


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"

class DQMStore;
class MonitorElement;

class LQeDQM : public edm::EDAnalyzer {
 public:

  /// Constructor
  LQeDQM(const edm::ParameterSet&);
  
  /// Destructor
  virtual ~LQeDQM();
  
  /// Inizialize parameters for histo binning
  void beginJob(edm::EventSetup const& iSetup);

  /// Get the analysis
  void analyze(const edm::Event&, const edm::EventSetup&);

  /// Save the histos
  void endJob(void);

  double calcDeltaPhi(double phi1, double phi2);

 private:

  // ----------member data ---------------------------
  
  DQMStore* theDbe;
  // Switch for verbosity
  std::string logTraceName;

  // Variables from config file
  std::string   theElecTriggerPathToPass;
  std::string   theMuonTriggerPathToPass;
  double    jetPtCut;
  double    elePtCut;
  double    metCut;
  edm::InputTag theTriggerResultsCollection;
  edm::InputTag theElectronCollectionLabel;
  edm::InputTag theCaloJetCollectionLabel;
  edm::InputTag theCaloMETCollectionLabel;

  // Histograms
  MonitorElement* h_e1j_invMass;
  MonitorElement* h_e2j_invMass;
  MonitorElement* h_ee_invMass;
  MonitorElement* h_ST;
  MonitorElement* h_STmet;
  MonitorElement* h_metj_transMass;
  MonitorElement* h_jet1_et;
  MonitorElement* h_jet2_et;
  MonitorElement* h_jet1_eta;
  MonitorElement* h_jet2_eta;
  MonitorElement* h_jet_count;
  MonitorElement* h_e1_et;
  MonitorElement* h_e2_et;
  MonitorElement* h_e1_eta;
  MonitorElement* h_e2_eta;
  MonitorElement* h_e1_phi;
  MonitorElement* h_ele_count;
  MonitorElement* h_met;
};
#endif
