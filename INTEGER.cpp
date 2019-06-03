#include <algorithm>
#include <array>
#include <iostream>
#include <unordered_map>

#include "INTEGER.hh"

#include "Cuts.hh"

namespace INTEGER{
  namespace {
    using mom_zE = std::vector<std::array<double,2>>; //!< (multiple) z & E
    using mom_yzE = std::unordered_map<double, mom_zE>; //!< (multiple) y
    using mom_map = std::unordered_map<double, mom_yzE>; //!< (multiple) x

    /// return true if double only has an integer part
    bool is_int(double d){
      double intpart;
      return std::modf(d, &intpart) == 0;
    }

    /// returns true if all component of a vector are pure integers
    template<class T>
    bool is_int(T const & d){
      for(auto p: d)
        if(!is_int(p)) return false;
      return true;
    }

    /// returns E from m^2
    double get_E2(double const p0, double const p1, double const p2, double const m2 = 0.){
      return sqrt(+p0*p0+p1*p1+p2*p2+m2);
    }

    /// generates mom_map of all allowed integer PSP up to @param maxp
    mom_map possible_mom_map(
        double const mass2,double const maxp, cuts_momenta const & pass_cuts){
      mom_map results;
      for(double i0=-maxp; i0<maxp;++i0){
        for(double i1=-maxp; i1<maxp;++i1){
          for(double i2=-maxp; i2<maxp;++i2){
            const double E(get_E2(i0, i1, i2, mass2));
            if(is_int(E) && pass_cuts({i0, i1, i2, E}))
              results[i0][i1].emplace_back(std::array<double,2>{i2,E});
          }
        }
      }
      return std::move(results);
    }

    /// Prints a vector/array
    template<class T>
    void print_array(T const & a){
      std::cout << "{";
      for(auto p=a.begin(); p<a.end()-1; ++p)
        std::cout << *p << ", ";
      std::cout << *(a.end()-1) << "}";
    }

    /// Gives the sive of a unfolded mom_map
    size_t mom_map_size(mom_map const & map){
      size_t ret=0;
      for(const auto & x: map){
        for(const auto & y: x.second){
          ret+=y.second.size();
        }
      }
      return ret;
    }

    /// generate incoming states
    incoming_vec construct_incoming(
      std::vector<vec4>::const_iterator const begin,
      std::vector<vec4>::const_iterator const end
    ){
      double xa = 0.;
      double xb = 0.;
      for(auto out = begin; out<end; ++out){
        xa += (*out)[3] - (*out)[2];
        xb += (*out)[3] + (*out)[2];
      }
      return {{{{0,0,-xa/2.,xa/2.}}, {{0,0,xb/2.,xb/2.}}}};
    }

    incoming_vec construct_incoming(std::vector<vec4> const & outgoing){
      return construct_incoming(outgoing.begin(), outgoing.end());
    }

    std::vector<particle_vec> construct_incoming(
      std::vector<double> const & vec_x, std::vector<double> const & vec_y,
      std::vector<double> const & vec_z, std::vector<double> const & vec_E,
      cuts_final const & cuts
    ){
      particle_vec particles(vec_x.size());
      for(size_t i=0; i<vec_x.size(); ++i){
        particles[i] = {vec_x[i], vec_y[i], vec_z[i], vec_E[i]};
      }
      auto const incoming(construct_incoming(particles));
      if(is_int(incoming) && cuts(incoming, particles)){
        particles.insert( particles.end(), incoming.begin(), incoming.end() );
        std::cout << "POSSIBLE FOUND" << std::endl;
        for(auto const & p: particles){
          print_array(p);
          std::cout << std::endl;
        }
        return {std::move(particles)};
      }
      return std::vector<particle_vec>();
    }

    /// recursive iteration over z & E
    std::vector<particle_vec> iterate_zE(
      std::vector<mom_zE const *> const & all_zE, size_t const index,
      std::vector<double> const & prev_x, std::vector<double> const & prev_y,
      std::vector<double> & prev_z, std::vector<double> & prev_E,
      cuts_final const & cuts
    ){
      std::vector<particle_vec> results;
      for(const auto & azE: *(all_zE[index])){
        prev_z.push_back(azE[0]);
        prev_E.push_back(azE[1]);
        std::vector<particle_vec> new_result;
        if(index >= all_zE.size()-1){
          // stop iteration over z & E and calculate incoming
          new_result = construct_incoming(prev_x, prev_y, prev_z, prev_E, cuts);
        } else{
          // recursion
          new_result = iterate_zE(all_zE, index+1, prev_x, prev_y, prev_z, prev_E, cuts);
        }
        if(new_result.size() > 0)
          results.insert(results.end(), new_result.begin(), new_result.end());
        prev_z.pop_back();
        prev_E.pop_back();
      }
      return results;
    }

    /// recursive iteration over y
    std::vector<particle_vec> iterate_y(
      std::vector<mom_yzE const *> const & all_yzE, size_t const index,
      std::vector<double> const & prev_x, std::vector<double> & prev_y,
      std::vector<mom_zE const *> & prev_zE,
      cuts_final const & cuts
    ){
      if(index >= all_yzE.size()-1){
        // stop iteration over y and iterate over z & E
        double ty = 0.;
        for( auto const & y: prev_y)
          ty-=y;
        auto const & ay = all_yzE[index]->find(ty);
        if(ay == all_yzE[index]->end()
            || std::find(prev_y.begin(), prev_y.end(), ay->first) != prev_y.end())
          return std::vector<particle_vec>();
        prev_y.push_back(ay->first);
        prev_zE.push_back(&(ay->second));
        std::vector<double> new_z, new_E;
        auto const & results = iterate_zE(prev_zE, 0, prev_x, prev_y, new_z, new_E, cuts);
        prev_zE.pop_back();
        prev_y.pop_back();
        return results;
      }
      // recursion
      std::vector<particle_vec> results;
      for(const auto & ay: *(all_yzE[index])){
        if(std::find(prev_y.begin(), prev_y.end(), ay.first) != prev_y.end()){
          continue;
        }
        prev_y.push_back(ay.first);
        prev_zE.push_back(&(ay.second));
        auto const & new_result = iterate_y(all_yzE, index+1, prev_x, prev_y, prev_zE, cuts);
        if(new_result.size() > 0)
          results.insert(results.end(), new_result.begin(), new_result.end());
        prev_zE.pop_back();
        prev_y.pop_back();
      }
      return results;

    }

    /// recursive iteration over x
    std::vector<particle_vec> iterate_x(
      std::vector<mom_map const *> const & all_mom, size_t const index,
      std::vector<double> & prev_x,
      std::vector<mom_yzE const *> & prev_y,
      cuts_final const & cuts
    ){
      if(index >= all_mom.size()-1){
        // stop iteration of x and begin iteration on y
        double tx = 0.;
        for( auto const & x: prev_x)
          tx-=x;
        auto const & ax = all_mom[index]->find(tx);
        if(ax == all_mom[index]->end()
            || std::find(prev_x.begin(), prev_x.end(), ax->first) != prev_x.end())
          return std::vector<particle_vec>();
        prev_x.push_back(ax->first);
        prev_y.push_back(&(ax->second));
        std::vector<double> new_y;
        std::vector<mom_zE const * > new_zE;
        auto const & new_result = iterate_y(prev_y, 0, prev_x, new_y, new_zE, cuts);
        prev_y.pop_back();
        prev_x.pop_back();
        return new_result;
      }
      // recursion
      std::vector<particle_vec> results;
      for(const auto & ax: *(all_mom[index])){
        if(std::find(prev_x.begin(), prev_x.end(), ax.first) != prev_x.end()){
          continue;
        }
        prev_x.push_back(ax.first);
        prev_y.push_back(&(ax.second));
        auto const & new_result = iterate_x(all_mom, index+1, prev_x, prev_y, cuts);
        results.insert(results.end(), new_result.begin(), new_result.end());
        prev_y.pop_back();
        prev_x.pop_back();
      }
      return results;
    }
  }

  std::vector<particle_vec> find_possible_recusive(
    double const maxp, size_t const njets,
    cuts_momenta const & jet_cuts, cuts_momenta const & higgs_cuts, cuts_final const & global_cuts
  ){
    constexpr double higgs_m2 = 88.*88.*2.; // mass square of the Higgs boson
    const auto p_parton{possible_mom_map(0, maxp, jet_cuts)};
    const auto p_higgs{possible_mom_map(higgs_m2, maxp, higgs_cuts)};
    std::cout << p_parton.size() << " " << p_higgs.size() << "\nreal: "
      << mom_map_size(p_parton) << " " << mom_map_size(p_higgs) << std::endl;
    std::vector<mom_map const *> all_maps;
    all_maps.push_back(&p_higgs);
    for(size_t i=0; i<njets; ++i){
      all_maps.push_back(&p_parton);
    }
    std::vector<double> new_x;
    std::vector<mom_yzE const * > new_y;
    return iterate_x(all_maps, 0, new_x, new_y, global_cuts);
  }
}

int main(){
  using namespace INTEGER;
  std::vector<particle_vec> results{find_possible_recusive(200, 4, cuts_jets, cuts_higgs, cuts_global)};
  std::cout << "Map Found " << results.size() << " possible momenta" << std::endl;
}
