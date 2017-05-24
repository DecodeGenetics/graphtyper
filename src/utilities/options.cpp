#include <iostream>

#include <graphtyper/utilities/options.hpp>

namespace gyper
{

Options *
Options::instance()
{
  // if (!_instance)
  //   _instance = new Options;

  return _instance;
}


Options::Options(){}


void
Options::print()
{
  std::cout << "[options] INFO: Verbosity is " << verbosity << std::endl;
}


Options * Options::_instance = new Options;

} // namespace gyper
