// GeneratorFilters/FakeTauFilter

#include "GeneratorFilters/FakeTauFilter.h"
#include "TruthUtils/PIDCodes.h"
#include "TruthUtils/HepMCHelpers.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

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
  declareProperty("FastJetPtmin", m_fastjet_pt_min=20000.);
  declareProperty("FastJetEtamax", m_fastjet_eta_max=2.8);
  declareProperty("TrueTrackPt", m_true_track_pt=5000.);
  declareProperty("MinTracksCore", m_min_trk_core=1);
  declareProperty("MaxTracksCore", m_max_trk_core=4);
  declareProperty("MinTracksIso", m_min_trk_iso=0);
  declareProperty("MaxTracksIso", m_max_trk_iso=2);
  declareProperty("CoreDr", m_core_dr=0.2);
  declareProperty("IsoDr", m_iso_dr=0.4);
  declareProperty("NumberOfFakeTaus", m_n_truthfakes=2);
}


StatusCode FakeTauFilter::filterEvent() {

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
  for (unsigned int ijet = 0; ijet < jets.size(); ijet++) {
    auto jet = jets[ijet];

    if (not is_good_jet(jet))
      continue;

    ATH_MSG_INFO("Consider jet "<< ijet 
		 << " out of "<< jets.size()
		 << " with pT = "<< jet.pt()
		 << " eta = " << jet.eta()
		 << " phi = " << jet.phi());

    int n_tracks_core = 0;
    int n_tracks_iso = 0;
    // loop over all selected input particles
    for (unsigned int ip=0; ip < identified_pseudo_jets.size(); ip++) {
      auto part = identified_pseudo_jets[ip];
      // check that the true particle is in the considered jet.
      if (indices[ip] == (int)ijet) {
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
		       << " is core = "<< is_core_track(part, jet)
		       << " is iso = " << is_iso_track(part, jet));
	  n_tracks_core += (int)is_core_track(part, jet);
	  n_tracks_iso += (int)is_iso_track(part, jet);
	}
      }
    }

    // check that the jet passes the track counting requirements
    if (n_tracks_core >= m_min_trk_core and n_tracks_core <= m_max_trk_core)
      if (n_tracks_iso >= m_min_trk_iso and n_tracks_iso <= m_max_trk_iso) 
      n_good_jets += 1;
  }

  if (n_good_jets >= m_n_truthfakes)
    setFilterPassed(true);
  else
    setFilterPassed(false);

  return StatusCode::SUCCESS;
}


bool FakeTauFilter::is_good_jet(const fastjet::PseudoJet & jet)
{
  // pt cut
  if (jet.pt() < m_fastjet_pt_min)
    return false;

  // eta cut
  if (fabs(jet.eta()) > m_fastjet_eta_max)
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


bool FakeTauFilter::is_core_track(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet)
{
  if (jet.delta_R(track) < m_core_dr) 
    return true;
  else
    return false;
}


bool FakeTauFilter::is_iso_track(const fastjet::PseudoJet & track, const fastjet::PseudoJet & jet)
{
  if (jet.delta_R(track) < m_iso_dr and not is_core_track(track, jet))
    return true;
  else
    return false;
}
