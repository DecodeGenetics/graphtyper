#include <iostream>

#include <graphtyper/utilities/options.hpp>

namespace gyper
{
Options * Options::instance()
{
  return _instance;
}

const Options * Options::const_instance()
{
  return _instance;
}

Options::Options()
{
}

Options * Options::_instance = new Options;

} // namespace gyper
