#pragma once

#include <array>
#include <functional>
#include <vector>

namespace INTEGER{
  using vec4 = std::array<double,4>; //!< (x,y,z,E)
  using particle_vec = std::vector<vec4>;
  using incoming_vec = std::array<vec4, 2>;
  using cuts_momenta = std::function<bool(vec4)>;
  using cuts_final = std::function<bool(incoming_vec, particle_vec)>;

  /// wrapper for recursive call
  std::vector<particle_vec> find_possible_recusive(
    double const maxp, size_t const njets,
    cuts_momenta const & jet_cuts, cuts_momenta const & boson_cuts, cuts_final const & global_cuts
  );
}
