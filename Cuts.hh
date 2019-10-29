#pragma once

#include <cmath>
#include <random>

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

  /// Prints a vector/array
  template<class T>
  void print_array(T const & a){
    for(auto p=begin(a); p<end(a)-1; ++p)
      std::cout <<std::setw(4)<< *p << ", ";
    std::cout <<std::setw(4)<< *(end(a)-1) << "}";
  }

}

/// cuts on a single jet
bool cuts_jets(INTEGER::vec4 const & mom){
  // constexpr double minpt = 30.;
  constexpr double minpt = 0.;
  constexpr double max_rap = 5.;
  return pt(mom)>minpt && rap(mom) < max_rap;
}

/// cuts on the Higgs boson
bool cuts_higgs(INTEGER::vec4 const & mom){
  static std::mt19937_64 ran;
  // if(ran()/double(ran.max()) > 1) return false;
  constexpr double max_rap = 5.;
  return rap(mom) < max_rap;
}

/// cuts on final configuration
bool cuts_global(INTEGER::incoming_vec const & incoming,
                 INTEGER::particle_vec const & outgoing
){
  // static std::mt19937_64 ran;
  // if(ran()/double(ran.max()) > 1) return false;
  constexpr double max_Ecms = 10000.;
  if(Ecms(incoming) > max_Ecms) return false;

  // INTEGER::vec4 boson = outgoing[0];
  // if(rap(boson) != rap(outgoing[5])) return false;

  // auto out = to_PseudoJets(outgoing.cbegin()+1, outgoing.cend());
  auto out = sorted_by_rapidity(to_PseudoJets(outgoing.cbegin(), outgoing.cend()));
  constexpr double jet_R = .4;
  constexpr double minpt = 30.;
  static const std::vector<int> search_idx{0, 1, 2, 3, 4, 5, 5};
  fastjet::ClusterSequence cs{out,
    fastjet::JetDefinition(fastjet::JetAlgorithm::antikt_algorithm, jet_R)};
  auto const jets = sorted_by_rapidity(cs.inclusive_jets(minpt));
  // if(jets.size() < 4) return false;
  std::vector<int> const jet_idx{ cs.particle_jet_indices(jets) };
  for(size_t i=0; i<search_idx.size(); ++i){
    if(search_idx[i]<-1) continue;
    if(search_idx[i] != jet_idx[i]) return false;
  }
  std::cout << "#jets " << jets.size() << " jet idx: ";
  for(auto const & idx: jet_idx)
    std::cout << idx << " ";
  std::cout << std::endl;
  for(auto const & p: out){
    std::cout << "          ev.outgoing.push_back({gluon, {";
    print_array(p.four_mom());
    std::cout << ", {}});" << std::endl;
  }
  for(size_t i=0; i<incoming.size(); ++i){
    std::cout << "          ev.incoming[" << i << "] = {gluon, {";
    print_array(incoming[i]);
    std::cout << ", {}};" << std::endl;
  }
  return true;
}
