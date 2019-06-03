#pragma once

#include <cmath>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "INTEGER.hh"

namespace {

  /// calculate Rapidity
  double rap(INTEGER::vec4 const & mom){
    return 1./2.*log( (mom[3]+mom[2]) / (mom[3]-mom[2]) );
  }

  /// calculate transverse momentum
  double pt(INTEGER::vec4 const & mom){
    return sqrt( mom[0]*mom[0] + mom[1]*mom[1] );
  }

  /// calculate center of mass energy
  double Ecms(INTEGER::incoming_vec const & incoming){
    return sqrt(2*(incoming[0][3]*incoming[0][3]+incoming[1][3]*incoming[1][3]));
  }

  fastjet::PseudoJet to_PseudoJet(INTEGER::vec4 const & p){
    return {p[0], p[1], p[2], p[3]};
  }

  template<class iterator_cbegin, class iterator_cend>
  std::vector<fastjet::PseudoJet> to_PseudoJets(
      iterator_cbegin const & begin, iterator_cend const & end
  ){
    std::vector<fastjet::PseudoJet> jets;
    for(auto p=begin; p<end; ++p){
      jets.emplace_back(to_PseudoJet(*p));
    }
    return jets;
  }
}

/// cuts on a single jet
bool cuts_jets(INTEGER::vec4 const & mom){
  constexpr double minpt = 30.;
  constexpr double max_rap = 5.;
  return pt(mom)>minpt && rap(mom) < max_rap;
}

/// cuts on the Higgs boson
bool cuts_higgs(INTEGER::vec4 const & mom){
  constexpr double max_rap = 5.;
  return rap(mom) < max_rap;
}

/// cuts on final configuration
bool cuts_global(INTEGER::incoming_vec const & incoming,
                 INTEGER::particle_vec const & outgoing
){
  constexpr double max_Ecms = 10000.;
  if(Ecms(incoming) > max_Ecms) return false;
  auto out = to_PseudoJets(outgoing.cbegin()+1, outgoing.cend());
  constexpr double jet_R = .4;
  fastjet::ClusterSequence cs{out,
    fastjet::JetDefinition(fastjet::JetAlgorithm::antikt_algorithm, jet_R)};
  return cs.inclusive_jets().size() == outgoing.size()-1;
}
