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

  double get_E2(double const p0, double const p1, double const p2, double const m2 = 0.){
    return sqrt(+p0*p0+p1*p1+p2*p2+m2);
  }

  std::vector<vec4> possible_mom(
      double const mass2, double const maxp,double const minpt){
    std::vector<vec4> results;
    for(double i0=-maxp; i0<maxp;++i0){
      for(double i1=-maxp; i1<maxp;++i1){
        if(i0*i0+i1*i1>minpt*minpt){
          for(double i2=-maxp; i2<maxp;++i2){
            const double E(get_E2(i0, i1, i2, mass2));
            if(is_int(E))
              results.emplace_back(vec4{i0,i1,i2,E});
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

  std::array<vec4, 2> construct_incoming(std::vector<vec4> const & outgoing){
    return construct_incoming(outgoing.begin(), outgoing.end());
  }

  template<class T>
  void print_array(T const & a){
    std::cout << "{";
    for(auto p=a.begin(); p<a.end()-1; ++p)
      std::cout << *p << ", ";
    std::cout << *(a.end()-1) << "}";
  }

  std::vector<t_result> find_possible(double const maxp, bool const check_in){
    const auto p_parton{possible_mom(0,maxp,30)};
    const auto p_higgs{possible_mom(126*126,maxp,0)};
    std::cout << "Possible integer phase space points for the partons: " << p_parton.size() << " & higgs: " << p_higgs.size() << std::endl;
    std::vector<t_result> results;
    for(const auto & a: p_higgs){
      for(const auto & b: p_parton){
        for(const auto & c: p_parton){
          // for(const auto & d: p_parton){
            double const t = a[0]+b[0]+c[0]/*+d[0]*/;
            // first component are sorted => t always increases
            if(t>0)
              break;
            if((t == 0) && (a[1]+b[1]+c[1]/*+d[1]*/ == 0)){
              t_result particles{a, b, c/*, d*/};
              auto const incoming(construct_incoming(particles));
              if(!check_in || is_int(incoming)){
                particles.insert( particles.end(), incoming.begin(), incoming.end() );
                std::cout << "POSSIBLE FOUND" << std::endl;
                for(auto const & p: particles){
                  print_array(p);
                  std::cout << std::endl;
                }
                results.emplace_back(std::move(particles));
              }
            }
          // }
        }
      }
    }
    return std::move(results);
  }
}

int main(){
  std::vector<t_result> results_normal{find_possible(200, true)};
  std::cout << "Found " << results_normal.size() << " possible momenta" << std::endl;
}
