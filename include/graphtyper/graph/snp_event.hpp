#pragma once

#include <cstdint>
#include <sstream>
#include <string>

#include <cereal/access.hpp>

namespace gyper
{
class SnpEvent
{
  friend class cereal::access;

public:
  uint32_t pos{0};
  char base{'N'};

  SnpEvent() = default;
  SnpEvent(uint32_t _pos, char _base);

  std::string to_string() const;

  /*********************
   * OPERATOR OVERLOAD *
   *********************/
  bool operator==(SnpEvent const & b) const;
  bool operator!=(SnpEvent const & b) const;
  bool operator<(SnpEvent const & b) const;

private:
  template <class Archive>
  void serialize(Archive & ar, const unsigned int);
};

} // namespace gyper
