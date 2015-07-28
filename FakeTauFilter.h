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
#include "TruthUtils/PIDUtils.h"
#include "TLorentzVector.h"
/// simple class to decorate a fastjet::PseudoJet
/// with charge and pdg_id info
class IdentifiedPseudoJet : public fastjet::PseudoJet {

 public:
  IdentifiedPseudoJet(const HepMC::GenParticle * p)
    : fastjet::PseudoJet(p->momentum().px(), 
			 p->momentum().py(), 
			 p->momentum().pz(), 
			 p->momentum().e())
    {
      m_pdgid = p->pdg_id();
      m_charge = MC::PID::charge(m_pdgid);
    }

  virtual ~IdentifiedPseudoJet() { }

  int pdg_id() const {return m_pdgid;}
  int charge() const {return m_charge;}

 private:
  int m_pdgid;
  int m_charge;
};

// class to hold simple EDM of the built tau 
class TruthFakeTau : public fastjet::PseudoJet 
{
 public:
  TruthFakeTau(const TLorentzVector & vec)
    : fastjet::PseudoJet(vec.Px(), vec.Py(), vec.Pz(), vec.E())
    {}

  int nTracks() const {return m_nTracks;}
  int nWideTracks() const {return m_nWideTracks;}
  
  void set_ntracks(const int & n) {m_nTracks=n;}
  void set_nwidetracks(const int & n) {m_nWideTracks=n;}

 private:
  int m_nTracks;
  int m_nWideTracks;
};

typedef std::vector<TruthFakeTau> TruthFakeTaus;

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
  double m_pt_core_min;
  double m_true_track_pt;
  int m_min_trk_core;
  int m_max_trk_core;
  int m_min_trk_iso;
  int m_max_trk_iso;
  double m_core_dr;
  double m_iso_dr;

  bool m_use_dr;
  double m_dr_taus;

  unsigned int m_n_truthfakes;
  TruthFakeTaus m_TruthFakeTaus;

  bool pass_deltaR(const TruthFakeTaus & taus, const double & drcut);

  bool is_good_jet(const fastjet::PseudoJet & jet);
  bool is_good_track(const IdentifiedPseudoJet & track);
  bool is_core(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet);
  bool is_iso(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet);

};


#endif
