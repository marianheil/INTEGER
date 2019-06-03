#pragma once

#include <cmath>

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

  /// calculate polar angle phi
  double phi(INTEGER::vec4 const & p) {
    double phi = atan2(p[1],p[2]);
    if (phi < 0.0) phi += 2.*M_PI;
    if (phi >= 2.*M_PI) phi -= 2.*M_PI;
    return phi;
  }

  /// calculate jet radius
  double dR(INTEGER::vec4 const & a, INTEGER::vec4 const & b){
    return std::abs( std::abs(rap(a)-rap(b)) - std::abs(phi(a)-phi(b)) );
  }

  /// calculate center of mass energy
  double Ecms(INTEGER::incoming_vec const & incoming){
    return sqrt(2*(incoming[0][3]*incoming[0][3]+incoming[1][3]*incoming[1][3]));
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
  constexpr double max_rap = 2.;
  return rap(mom) < max_rap;
}

/// cuts on final configuration
bool cuts_global(INTEGER::incoming_vec const & incoming, INTEGER::particle_vec const & outgoing){
  constexpr double min_dR = .4;
  constexpr double max_Ecms = 10000.;
  for(size_t i=0; i<outgoing.size(); ++i){
    for(size_t j=0; j<i; ++j){
      if(dR(outgoing[i], outgoing[j]) < min_dR)
        return false;
    }
  }
  return Ecms(incoming) < max_Ecms;
}
