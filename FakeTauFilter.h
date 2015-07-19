// --------------------------------------------------
// 
// File:  GeneratorFilters/FakeTauFilter.h
// Description: Filters based on the presence of a jet 
//              that is likely to fake an hadronic tau
//
// Authors:
//         Q Buat:  July 2015

#ifndef GENERATORFILTERSFAKETAUFILTER_H
#define GENERATORFILTERSFAKETAUFILTER_H

#include "GeneratorModules/GenFilter.h"
#include "TruthUtils/FastJet.h"

class IdentifiedPseudoJet : public fastjet::PseudoJet {

 public:
  IdentifiedPseudoJet(const HepMC::GenParticle * p); 
  virtual ~IdentifiedPseudoJet() { }

  int pdg_id() const {return m_pdgid;}
  int charge() const {return m_charge;}

 private:
  int m_pdgid;
  int m_charge;
};

/// Filter events based on presence of charged leptons
class FakeTauFilter : public GenFilter {

public:

  /// Constructor
  FakeTauFilter(const std::string& name, ISvcLocator* pSvcLocator);

  /// Destructor
  virtual ~FakeTauFilter() { }

  /// Initialize
  virtual StatusCode filterInitialize() {
    return StatusCode::SUCCESS;
  }

  /// Finalize
  virtual StatusCode filterFinalize() {
    return StatusCode::SUCCESS;
  }

  /// Do the filtering
  virtual StatusCode filterEvent();


private:

  double m_fastjet_cone_size;
  double m_fastjet_pt_min;
  double m_fastjet_eta_max;
  double m_true_track_pt;
  int m_max_trk_core;
  int m_max_trk_iso;
  double m_core_dr;
  double m_iso_dr;
  int m_n_truthfakes;

  bool is_good_jet(const fastjet::PseudoJet & jet);
  bool is_good_track(const IdentifiedPseudoJet & track);
  bool is_core_track(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet);
  bool is_iso_track(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet);

};


#endif
