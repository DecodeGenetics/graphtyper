#pragma once

#include <iostream>

#include <graphtyper/utilities/logging.hpp>

inline void setup_logging()
{
  if (gyper::log_singleton == nullptr)
    gyper::log_singleton =
      std::unique_ptr<gyper::log_singleton_t>{new gyper::log_singleton_t{gyper::log_severity::debug, std::clog}};
}

// this is a global variable during whose initialisation the log_singleton is also initialised
inline int IGNOREME = []()
{
  setup_logging();
  return 0;
}();
