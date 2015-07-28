// GeneratorFilters/FakeTauFilter

#include "GeneratorFilters/FakeTauFilter.h"
#include "TruthUtils/PIDCodes.h"
#include "TruthUtils/HepMCHelpers.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <cmath>

//  ////////////
// IdentifiedPseudoJet::IdentifiedPseudoJet(const HepMC::GenParticle * p) 
//   : fastjet::PseudoJet(p->momentum().px(), 
// 		       p->momentum().py(), 
// 		       p->momentum().pz(), 
// 		       p->momentum().e())
// {
//   m_pdgid = p->pdg_id();
//   m_charge = MC::PID::charge(m_pdgid);
// }

FakeTauFilter::FakeTauFilter(const std::string& name, ISvcLocator* pSvcLocator)
  : GenFilter(name,pSvcLocator) 
{
  declareProperty("FastJetConeSize", m_fastjet_cone_size=0.4);
  declareProperty("FastJetPtmin", m_fastjet_pt_min=25000.);
  declareProperty("FastJetEtamax", m_fastjet_eta_max=2.8);
  declareProperty("TrueTrackPt", m_true_track_pt=1000.);
  declareProperty("MinPtCore", m_pt_core_min=25000.);
  declareProperty("MinTracksCore", m_min_trk_core=1);
  declareProperty("MaxTracksCore", m_max_trk_core=4);
  declareProperty("MinTracksIso", m_min_trk_iso=0);
  declareProperty("MaxTracksIso", m_max_trk_iso=1);
  declareProperty("CoreDr", m_core_dr=0.2);
  declareProperty("IsoDr", m_iso_dr=0.4);
  declareProperty("NumberOfFakeTaus", m_n_truthfakes=2);
  declareProperty("UseDr", m_use_dr=true);
  declareProperty("DrTaus", m_dr_taus=2.8);

}


StatusCode FakeTauFilter::filterEvent() {

  // clear the list of built taus
  m_TruthFakeTaus.clear();

  // Loop over all particles in the event and build up the truth jets
  const HepMC::GenEvent* genEvt = event();

  std::vector<IdentifiedPseudoJet> identified_pseudo_jets;
  for(HepMC::GenEvent::particle_const_iterator part = genEvt->particles_begin(); part != genEvt->particles_end(); ++part) {

    // remove unstables
    if (not MC::isGenStable(*part))
      continue;

    // No muons
    if ((*part)->pdg_id() == MC::PID::MUON or (*part)->pdg_id() == MC::PID::ANTIMUON)
      continue;

    // No electron neutrinos 
    if ((*part)->pdg_id() == MC::PID::NU_E or (*part)->pdg_id() == MC::PID::NU_EBAR)
      continue;

    // No muon neutrinos 
    if ((*part)->pdg_id() == MC::PID::NU_MU or (*part)->pdg_id() == MC::PID::NU_MUBAR)
      continue;

    // No tau neutrinos 
    if ((*part)->pdg_id() == MC::PID::NU_TAU or (*part)->pdg_id() == MC::PID::NU_TAUBAR)
      continue;

    identified_pseudo_jets.push_back(IdentifiedPseudoJet((*part)));
  }  // end of the loop over the truth particles

  // fastjet clustering
  const fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, m_fastjet_cone_size); 
  fastjet::ClusterSequence cseq(identified_pseudo_jets, jet_def);
  auto jets = sorted_by_pt(cseq.inclusive_jets(m_fastjet_pt_min));
  auto indices = cseq.particle_jet_indices(jets);

  int n_good_jets = 0;
  // loop over all the candidates built with fastjet
  for (unsigned int ijet = 0; ijet < jets.size(); ijet++) {
    auto jet = jets[ijet];

    if (not is_good_jet(jet))
      continue;

    ATH_MSG_INFO("Consider jet "<< ijet 
		 << " out of "<< jets.size()
		 << " with pT = "<< jet.pt()
		 << " eta = " << jet.eta()
		 << " pseudorap = " << jet.pseudorapidity()
		 << " phi = " << jet.phi());

    int n_tracks_core = 0;
    int n_tracks_iso = 0;
    TLorentzVector jet_core;
    jet_core.SetPxPyPzE(0., 0., 0., 0.);
    // loop over all selected input particles
    for (unsigned int ip=0; ip < identified_pseudo_jets.size(); ip++) {
      auto part = identified_pseudo_jets[ip];
      // check that the true particle is in the considered jet.
      if (indices[ip] == (int)ijet) {
	// compute the core pT
	if (is_core(part, jet)) {
	  TLorentzVector jet_temp;
	  jet_temp.SetPtEtaPhiM(part.pt(), 
				part.pseudorapidity(),
				part.phi(),
				part.m());
	  jet_core += jet_temp;
	}
	// check that the selected particle 
	// makes a good track candidate 
	if (is_good_track(part)) {
	  ATH_MSG_INFO("\t particle "<< ip
		       << " with pT = "<< part.pt()
		       << " eta = " << part.eta()
		       << " phi = " << part.phi()
		       << " pdg id = " << part.pdg_id()
		       << " good track = "<< is_good_track(part));
	  ATH_MSG_INFO("\t\t particle "<< ip
		       << " is core = "<< is_core(part, jet)
		       << " is iso = " << is_iso(part, jet));
	  n_tracks_core += (int)is_core(part, jet);
	  n_tracks_iso += (int)is_iso(part, jet);
	}
      }
    } // end of the loop over all the input particles

    ATH_MSG_INFO("\t\t\t jet "<< ijet 
		 << " out of "<< jets.size()
		 << " with pT = "<< jet.pt()
		 << " and core pt "<<jet_core.Pt());

    TruthFakeTau tau(jet_core);
    tau.set_ntracks(n_tracks_core);
    tau.set_nwidetracks(n_tracks_iso);

    if (tau.pt() < m_pt_core_min)
      continue;

    // check that the jet passes the track counting requirements
    if (tau.nTracks() >= m_min_trk_core and tau.nTracks() <= m_max_trk_core) {
      if (tau.nWideTracks() >= m_min_trk_iso and tau.nWideTracks() <= m_max_trk_iso) {
	n_good_jets++;
	m_TruthFakeTaus.push_back(tau);
      }
    }

  }// end of the loop over the fastjet jets

  ATH_MSG_INFO("simple couting = " << n_good_jets
	       << " and vector size = "<< m_TruthFakeTaus.size());
  

  if (m_TruthFakeTaus.size() < m_n_truthfakes) {
    setFilterPassed(false);
  } else {
    if (m_use_dr) {
      if (pass_deltaR(m_TruthFakeTaus, m_dr_taus))
	setFilterPassed(true);
      else
	setFilterPassed(false);
    } else {
	setFilterPassed(true);
    }
  }

  return StatusCode::SUCCESS;
}

bool FakeTauFilter::pass_deltaR(const TruthFakeTaus & taus, const double & drcut)
{

  for (unsigned int i1 = 0; i1 < taus.size(); i1++) {
    auto p1 = taus[i1];
    TVector3 v1;
    v1.SetPtEtaPhi(p1.pt(), p1.pseudorapidity(), p1.phi());
    for (unsigned int i2 = i1 + 1; i2 < taus.size(); i2++) {
      auto p2 = taus[i2];
      TVector3 v2;
      v2.SetPtEtaPhi(p2.pt(), p2.pseudorapidity(), p2.phi());
      if (v1.DeltaR(v2) < drcut)
	return true;
    }
  }
  return false;


}



bool FakeTauFilter::is_good_jet(const fastjet::PseudoJet & jet)
{
  // pt cut
  if (jet.pt() < m_fastjet_pt_min)
    return false;

  // eta cut
  if (fabs(jet.pseudorapidity()) > m_fastjet_eta_max)
    return false;

  return true;
}


bool FakeTauFilter::is_good_track(const IdentifiedPseudoJet & track)
{
  // charge cut
  if (track.charge() == 0)
    return false;

  // pt cut
  if (track.pt() < m_true_track_pt)
    return false;

  return true;
}


bool FakeTauFilter::is_core(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet)
{
  TVector3 jet_vec;
  TVector3 track_vec;
  jet_vec.SetPtEtaPhi(jet.pt(), jet.pseudorapidity(), jet.phi());
  track_vec.SetPtEtaPhi(track.pt(), track.pseudorapidity(), track.phi());
  ATH_MSG_INFO("\t\t\t Delta R = " << jet_vec.DeltaR(track_vec));
  if (jet_vec.DeltaR(track_vec) < m_core_dr) 
    return true;
  else
    return false;
}


bool FakeTauFilter::is_iso(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet)
{
  TVector3 jet_vec;
  TVector3 track_vec;
  jet_vec.SetPtEtaPhi(jet.pt(), jet.pseudorapidity(), jet.phi());
  track_vec.SetPtEtaPhi(track.pt(), track.pseudorapidity(), track.phi());
  if (jet_vec.DeltaR(track_vec) < m_iso_dr and not is_core(track, jet))
    return true;
  else
    return false;
}
