#pragma once

#include <string>
#include <vector>

#include <boost/serialization/access.hpp>

#include <graphtyper/constants.hpp>


namespace gyper
{

class RefNode;
class VarNode;
class Alt;


class Label
{
  friend class boost::serialization::access;
  friend class gyper::RefNode;
  friend class gyper::VarNode;

public:
  uint32_t order{0};
  std::vector<char> dna{};
  uint32_t variant_num{0};

  Label() noexcept;
  Label(Label const & l) noexcept;
  Label(Label && l) noexcept;
  Label(uint32_t const order, std::vector<char> && dna, uint32_t const variant_num) noexcept;
  Label(uint32_t const order, Alt && alt, uint32_t const variant_num);

  /**
   * CLASS INFORMATION
   */
  uint32_t reach() const;

private:
  template <typename Archive>
  void serialize(Archive & ar, const unsigned int);
};

} // namespace gyper
