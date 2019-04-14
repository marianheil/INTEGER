#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <typeinfo>

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

  size_t mom_map_size(mom_map const & map){
    size_t ret=0;
    for(const auto & x: map){
      for(const auto & y: x.second){
        ret+=y.second.size();
      }
    }
    return ret;
  }

  std::vector<t_result> construct_incoming(
    std::vector<double> const & vec_x, std::vector<double> const & vec_y,
    std::vector<double> const & vec_z, std::vector<double> const & vec_E
  ){
    t_result particles(vec_x.size());
    for(size_t i=0; i<vec_x.size(); ++i){
      particles[i] = {vec_x[i], vec_y[i], vec_z[i], vec_E[i]};
    }
    auto const incoming(construct_incoming(particles));
    if(is_int(incoming)){
      particles.insert( particles.end(), incoming.begin(), incoming.end() );
      std::cout << "POSSIBLE FOUND" << std::endl;
      for(auto const & p: particles){
        print_array(p);
        std::cout << std::endl;
      }
      return {std::move(particles)};
    }
    return std::vector<t_result>();
  }

  std::vector<t_result> iterate_zE(
    std::vector<mom_zE const *> const & all_zE, size_t const index,
    std::vector<double> const & prev_x, std::vector<double> const & prev_y,
    std::vector<double> & prev_z, std::vector<double> & prev_E
  ){
    std::vector<t_result> results;
    for(const auto & azE: *(all_zE[index])){
      prev_z.push_back(azE[0]);
      prev_E.push_back(azE[1]);
      std::vector<t_result> new_result;
      if(index >= all_zE.size()-1){
        // stop iteration over z & E and calculate incoming
        new_result = construct_incoming(prev_x, prev_y, prev_z, prev_E);
      } else{
        new_result = iterate_zE(all_zE, index+1, prev_x, prev_y, prev_z, prev_E);
      }
      if(new_result.size() > 0)
        results.insert(results.end(), new_result.begin(), new_result.end());
      prev_z.pop_back();
      prev_E.pop_back();
    }
    return results;
  }

  std::vector<t_result> iterate_y(
    std::vector<mom_yzE const *> const & all_yzE, size_t const index,
    std::vector<double> const & prev_x, std::vector<double> & prev_y,
    std::vector<mom_zE const *> & prev_zE
  ){
    if(index >= all_yzE.size()-1){
      // stop iteration over y and iterate over z & E
      double ty = 0.;
      for( auto const & y: prev_y)
        ty+=y;
      auto const & ay = all_yzE[index]->find(ty);
      if(ay == all_yzE[index]->end())
        return std::vector<t_result>();
      prev_y.push_back(ay->first);
      prev_zE.push_back(&(ay->second));
      std::vector<double> new_z, new_E;
      auto const & results = iterate_zE(prev_zE, 0, prev_x, prev_y, new_z, new_E);
      prev_zE.pop_back();
      prev_y.pop_back();
      return results;
    }
    std::vector<t_result> results;
    for(const auto & ay: *(all_yzE[index])){
      prev_y.push_back(ay.first);
      prev_zE.push_back(&(ay.second));
      auto const & new_result = iterate_y(all_yzE, index+1, prev_x, prev_y, prev_zE);
      if(new_result.size() > 0)
        results.insert(results.end(), new_result.begin(), new_result.end());
      prev_zE.pop_back();
      prev_y.pop_back();
    }
    return results;

  }

  std::vector<t_result> iterate_x(
    std::vector<mom_map const *> const & all_mom, size_t const index,
    std::vector<double> & prev_x,
    std::vector<mom_yzE const *> & prev_y
  ){
    if(index >= all_mom.size()-1){
      // stop iteration of x and begin iteration on y
      double tx = 0.;
      for( auto const & x: prev_x)
        tx-=x;
      auto const & ax = all_mom[index]->find(tx);
      if(ax == all_mom[index]->end())
        return std::vector<t_result>();
      prev_x.push_back(ax->first);
      prev_y.push_back(&(ax->second));
      std::vector<double> new_y;
      std::vector<mom_zE const * > new_zE;
      auto const & new_result = iterate_y(prev_y, 0, prev_x, new_y, new_zE);
      prev_y.pop_back();
      prev_x.pop_back();
      return new_result;
    }
    // call recusively
    std::vector<t_result> results;
    for(const auto & ax: *(all_mom[index])){
      prev_x.push_back(ax.first);
      prev_y.push_back(&(ax.second));
      auto const & new_result = iterate_x(all_mom, index+1, prev_x, prev_y);
      results.insert(results.end(), new_result.begin(), new_result.end());
      prev_y.pop_back();
      prev_x.pop_back();
    }
    return results;
  }

  std::vector<t_result> find_possible_map(double const maxp, bool const check_in){
    const auto p_parton{possible_mom_map(0,maxp,30)};
    const auto p_higgs{possible_mom_map(88.*sqrt(2.),maxp,0)};
    std::cout << p_parton.size() << " " << p_higgs.size() << "\nreal: "
      << mom_map_size(p_parton) << " " << mom_map_size(p_higgs) << std::endl;
    std::vector<t_result> results;
    for(const auto & ax: p_higgs)
      for(const auto & bx: p_parton)
      for(const auto & cx: p_parton){
        // TODO make this calls nicer
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
            for(const auto & azE: ay.second)
              for(const auto & bzE: by.second)
              for(const auto & czE: cy.second)
                for(const auto & dzE: dy->second){
                  t_result particles;
                  particles.emplace_back(vec4{ax.first,  ay.first,  azE[0], azE[1]});
                  particles.emplace_back(vec4{bx.first,  by.first,  bzE[0], bzE[1]});
                  particles.emplace_back(vec4{cx.first,  cy.first,  czE[0], czE[1]});
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
