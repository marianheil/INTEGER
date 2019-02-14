#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <typeinfo>

namespace {
  using vec4 = std::array<double,4> ;
  using t_result = std::vector<vec4> ;
  using mom_map = std::unordered_map<double, // x
      std::unordered_map<double,  // (multiple) y
        std::vector<std::array<double,2>> // (multiple) (z, E)
      >
    >;

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

  // bool is_int(vec4 const & d){
  //   for(auto p: d)
  //     if(!is_int(p)) return false;
  //   return true;
  // }

  // bool is_int(t_result const & d){
  //   for(auto p: d)
  //     if(!is_int(p)) return false;
  //   return true;
  // }

  double m2(vec4 const & p){
    return -p[0]*p[0]-p[1]*p[1]-p[2]*p[2]+p[3]*p[3];
  }

  double m(vec4 const & p) {
    return sqrt(m2(p));
  }

  double get_E(double const p0, double const p1, double const p2, double const m = 0.){
    return sqrt(+p0*p0+p1*p1+p2*p2+m*m);
  }

  std::vector<vec4> possible_mom(
      double const mass,double const maxp,double const minpt){
    std::vector<vec4> results;
    for(double i0=-maxp; i0<maxp;++i0){
      for(double i1=-maxp; i1<maxp;++i1){
        if(i0*i0+i1*i1>minpt*minpt){
          for(double i2=-maxp; i2<maxp;++i2){
            const double E(get_E(i0, i1, i2, mass));
            if(is_int(E))
              results.emplace_back(vec4{i0,i1,i2,E});
          }
        }
      }
    }
    return std::move(results);
  }

  mom_map possible_mom_map(
      double const mass,double const maxp,double const minpt){
    mom_map results;
    for(double i0=-maxp; i0<maxp;++i0){
      for(double i1=-maxp; i1<maxp;++i1){
        if(i0*i0+i1*i1>minpt*minpt){
          for(double i2=-maxp; i2<maxp;++i2){
            const double E(get_E(i0, i1, i2, mass));
            if(is_int(E))
              results[i0][i1].emplace_back(std::array<double,2>{i2,E});
          }
        }
      }
    }
    return std::move(results);
  }

  std::array<vec4, 2> construct_incoming(std::vector<vec4> const & outgoing){
    double xa = 0.;
    double xb = 0.;
    for(auto const & out: outgoing){
      xa += out[3] - out[2];
      xb += out[3] + out[2];
    }
    return {{{{0,0,-xa/2.,xa/2.}}, {{0,0,xb/2.,xb/2.}}}};
  }

  template<class T>
  void print_array(T const & a){
    std::cout << "{";
    for(const auto p: a)
      std::cout << p << ", ";
    std::cout << "}";
  }

  std::vector<t_result> find_possible(double const maxp, bool const check_in){
    const auto p_parton{possible_mom(0,maxp,30)};
    const auto p_higgs{possible_mom(125,maxp,0)};
    std::cout << p_parton.size() << " " << p_higgs.size() << std::endl;
    std::vector<t_result> results;
    for(const auto & a: p_higgs){
      for(const auto & b: p_parton){
        // for(const auto & c: p_parton){
          // for(const auto & d: p_parton){
            double const t = a[0]+b[0]/*+c[0]*//*+d[0]*/;
            // first component are sorted => t always increases
            if(t>0)
              break;
            if((t == 0) && (a[1]+b[1]/*+c[1]*//*+d[1]*/ == 0)){
              t_result particles{a,b/*,c*//*,d*/};
              auto const incoming(construct_incoming(particles));
              if(!check_in || is_int(incoming)){
                particles.insert( particles.end(), incoming.begin(), incoming.end() );
                // std::cout << "POSSIBLE FOUND" << std::endl;
                // for(auto const & p: particles){
                //   print_array(p);
                //   std::cout << std::endl;
                // }
                results.emplace_back(std::move(particles));
              }
            }
          // }
        // }
      }
    }
    return std::move(results);
  }

  // template<class T>
  // std::vector<T> possible_xy(
  //     std::vector<T> const & outs,
  //     double const val /*momentum previous particle*/,
  //     size_t const pos /*current position*/
  // ){
  //   if(outs.size()==pos+1){ // if last particle
  //     auto const & res = outs[pos].first.find(val);
  //     // check if value exists
  //     if(res == outs[pos].first.end()) return {};
  //     // and return it
  //     return {res};
  //   }
  //   std::vector<T> returns;
  //   // else find all possible results of one stage up
  //   for(auto const & it: outs[pos]){
  //     auto const pref = possible_xy(outs, val+it.first, pos+1);
  //     if(pref.size()!=0)
  //       returns.emplace_back({it, pref});
  //   }
  //   return std::move(returns);
  // }

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
    const auto p_parton{possible_mom_map(0,maxp,30)};
    const auto p_higgs{possible_mom_map(125,maxp,0)};
    std::cout << p_parton.size() << " " << p_higgs.size() << "\nreal: "
      << mom_map_size(p_parton) << " " << mom_map_size(p_higgs) << std::endl;
    std::vector<t_result> results;
    for(const auto & ax: p_higgs)
      for(const auto & bx: p_parton){
        // TODO make this calls nicer
        double const tx = -(ax.first+bx.first/*+cx.first*/);
        auto const & dx = p_parton.find(tx);
        if(dx == p_parton.end())
          continue;
        for(const auto & ay: ax.second)
          for(const auto & by: bx.second){
            double const ty = -(ay.first+by.first/*+cy.first*/);
            auto const & dy = dx->second.find(ty);
            if(dy == dx->second.end())
              continue;
            for(const auto & azE: ay.second)
              for(const auto & bzE: by.second)
                for(const auto & dzE: dy->second){
                  t_result particles;
                  particles.emplace_back(vec4{ax.first,  ay.first,  azE[0], azE[1]});
                  particles.emplace_back(vec4{bx.first,  by.first,  bzE[0], bzE[1]});
                  particles.emplace_back(vec4{dx->first, dy->first, dzE[0], dzE[1]});
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
          }
      }
    return std::move(results);
  }
}

/**
 * optimisation:
 * 1. for "d" use (unordered) map( x->map( y->vector(array(z,E)) ) )
 * 2. if(d.find(a[0]+b[0]+c[0]) != d.end()) ...
 * 3. if(d.find(a[1]+b[1]+c[1]) != d.end()) ...
 *
 */

int main(){
  std::vector<t_result> results{find_possible_map(100, true)};
  std::cout << "Map Found " << results.size() << " possible momenta" << std::endl;
  // std::vector<t_result> results_normal{find_possible(200, false)};
  // std::cout << "Found " << results_normal.size() << " possible momenta" << std::endl;
}
