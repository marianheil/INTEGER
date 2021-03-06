#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <unordered_map>
#include <vector>

namespace {
  using vec4 = std::array<double,4>;
  using t_result = std::vector<vec4>;
  using mom_zE = std::vector<std::array<double,2>>; // (multiple) z & E
  using mom_yzE = std::unordered_map<double, mom_zE>; // (multiple) y
  using mom_map = std::unordered_map<double, mom_yzE>; // (multiple) x

  bool is_int(double d){
    double intpart;
    return std::modf(d, &intpart) == 0;
  }

  template<class T>
  bool is_int(T const & d){
    for(auto p: d)
      if(!is_int(p)) return false;
    return true;
  }

  double rap(double const pz, double const E){
    return 1./2.*log(1.*(E+pz)/(E-pz));
  }

  double rap(vec4 const & p) {
    return rap(p[2], p[3]);
  }

  double phi(vec4 const & p) {
    double phi = atan2(p[1],p[2]);
    if (phi < 0.0) phi += 2.*M_PI;
    if (phi >= 2.*M_PI) phi -= 2.*M_PI;
    return phi;
  }

  double dR(vec4 const & a, vec4 const & b){
    return std::abs(std::abs(rap(a)-rap(b))-std::abs(phi(a)-phi(b)));
  }

  double get_E2(double const p0, double const p1, double const p2, double const m2 = 0.){
    return sqrt(+p0*p0+p1*p1+p2*p2+m2);
  }

  mom_map possible_mom_map(
      double const mass2,double const maxp, double const max_rap, double const minpt){
    mom_map results;
    for(double i0=-maxp; i0<maxp;++i0){
      for(double i1=-maxp; i1<maxp;++i1){
        if(i0*i0+i1*i1>minpt*minpt){
          for(double i2=-maxp; i2<maxp;++i2){
            const double E(get_E2(i0, i1, i2, mass2));
            if(is_int(E) && rap(i2, E) < max_rap)
              results[i0][i1].emplace_back(std::array<double,2>{i2,E});
          }
        }
      }
    }
    return std::move(results);
  }

  std::array<vec4, 2> construct_incoming(std::vector<vec4>::const_iterator const begin, std::vector<vec4>::const_iterator const end){
    double xa = 0.;
    double xb = 0.;
    for(auto out = begin; out<end; ++out){
      xa += (*out)[3] - (*out)[2];
      xb += (*out)[3] + (*out)[2];
    }
    return {{{{0,0,-xa/2.,xa/2.}}, {{0,0,xb/2.,xb/2.}}}};
  }

  template<class T>
  void print_array(T const & a){
    std::cout << "{";
    for(auto p=a.begin(); p<a.end()-1; ++p)
      std::cout << *p << ", ";
    std::cout << *(a.end()-1) << "}";
  }

  size_t mom_map_size(mom_map const & map){
    size_t ret=0;
    for(const auto & x: map){
      for(const auto & y: x.second){
        ret+=y.second.size();
      }
    }
    return ret;
  }

  std::vector<t_result> find_possible_map(double const maxp, bool const check_in){
    constexpr double max_rap = 5.; // maximal rapidity
    constexpr double min_Rjj = 0.4; // minimal R distance between 2 jets
    constexpr double max_Ecm = 10000.; // maximal center of mass Energy
    constexpr double higgs_m2 = 88.*88.*2.; // mass square of the Higgs boson
    const auto p_parton{possible_mom_map(0, maxp, max_rap, 30)};
    const auto p_higgs{possible_mom_map(higgs_m2, maxp, max_rap, 0)};
    std::cout << p_parton.size() << " " << p_higgs.size() << "\nreal: "
      << mom_map_size(p_parton) << " " << mom_map_size(p_higgs) << std::endl;
    std::vector<t_result> results;
    for(const auto & ax: p_higgs)
      for(const auto & bx: p_parton){
        // enforce different px for higgs and the first jet (make it look more random)
        if(bx.first == ax.first) continue;
        for(const auto & cx: p_parton){
          double const tx = -(ax.first+bx.first+cx.first);
          auto const & dx = p_parton.find(tx);
          if(dx == p_parton.end())
            continue;
          for(const auto & ay: ax.second)
            for(const auto & by: bx.second)
              for(const auto & cy: cx.second){
                double const ty = -(ay.first+by.first+cy.first);
                auto const & dy = dx->second.find(ty);
                if(dy == dx->second.end())
                  continue;

                t_result particles;
                particles.resize(6);
                for(const auto & azE: ay.second){
                  particles[0]=vec4{ax.first,  ay.first,  azE[0], azE[1]};
                  for(const auto & bzE: by.second){
                    particles[1]=vec4{bx.first,  by.first,  bzE[0], bzE[1]};
                    for(const auto & czE: cy.second){
                      particles[2]=vec4{cx.first,  cy.first,  czE[0], czE[1]};
                      if( dR( particles[2], particles[1] ) < min_Rjj) continue;
                      for(const auto & dzE: dy->second){
                        particles[3]=vec4{dx->first, dy->first, dzE[0], dzE[1]};
                        if(( dR(particles[3], particles[2]) < min_Rjj) || ( dR(particles[3], particles[1]) < min_Rjj)) continue;
                        auto const incoming(construct_incoming(particles.begin(), particles.end()-2));
                        if(!check_in ||
                            ( is_int(incoming)
                              && (2*(incoming[0][3]*incoming[0][3]+incoming[1][3]*incoming[1][3])<max_Ecm*max_Ecm) )
                        ){
                          particles[particles.size()-2]=std::move(incoming[0]);
                          particles[particles.size()-1]=std::move(incoming[1]);
                          std::cout << "POSSIBLE FOUND" << std::endl;
                          for(auto const & p: particles){
                            print_array(p);
                            std::cout << std::endl;
                          }
                          results.push_back(particles);
                        }
                      }
                    }
                  }
                }
              }
        }
      }
    return results;
  }
}

int main(){
  std::vector<t_result> results{find_possible_map(100, true)};
  std::cout << "Map Found " << results.size() << " possible momenta" << std::endl;
}
