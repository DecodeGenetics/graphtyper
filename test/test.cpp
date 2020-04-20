#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>

#include <graphtyper/constants.hpp>

#include "test.hpp"

namespace gyper
{

Logger::Logger()
{
  namespace blog = boost::log;

  blog::core::get()->set_filter
  (
    blog::trivial::severity >= blog::trivial::debug
  );

  blog::add_common_attributes();
  blog::register_simple_formatter_factory<blog::trivial::severity_level, char>("Severity");
  std::string log_format = "[%TimeStamp%] <%Severity%> %Message%";

  boost::log::add_console_log(std::clog,
                              boost::log::keywords::auto_flush = true,
                              boost::log::keywords::format = log_format);

  BOOST_LOG_TRIVIAL(debug) << "[" << __HERE__ << "] Boost.Log has been set up";
}


Logger * Logger::_instance = new Logger;

} // namespace gyper
