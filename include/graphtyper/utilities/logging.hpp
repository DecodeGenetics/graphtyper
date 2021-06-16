#pragma once

#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <mutex>
#include <ostream>

namespace gyper
{
// LOG SEVERITY
enum class log_severity
{
  debug,
  info,
  warning,
  error
};

inline std::string_view severity2string(log_severity const severity)
{
  switch (severity)
  {
  case log_severity::debug:
    return "debug";
  case log_severity::info:
    return "info";
  case log_severity::warning:
    return "warning";
  case log_severity::error:
    return "error";
  }
  return "debug";
}

// TYPE FOR STORING LOGGING STATE
class log_singleton_t
{
private:
  // Pointer with customisable delete behaviour
  using stream_ptr_t = std::unique_ptr<std::ostream, std::function<void(std::ostream *)>>;
  // Stream deleter that does nothing (no ownership assumed).
  static void stream_deleter_noop(std::ostream *)
  {
  }
  // Stream deleter with default behaviour (ownership assumed).
  static void stream_deleter_default(std::ostream * ptr)
  {
    delete ptr;
  }

public:
  log_singleton_t() = delete;
  log_singleton_t(log_singleton_t const &) = delete;
  log_singleton_t(log_singleton_t &&) = delete;
  log_singleton_t & operator=(log_singleton_t const &) = delete;
  log_singleton_t & operator=(log_singleton_t &&) = delete;

  // construct from existing stream
  log_singleton_t(log_severity _severity, std::ostream & _sink) : severity{_severity}, sink{&_sink, stream_deleter_noop}
  {
  }

  // construct from filename
  log_singleton_t(log_severity _severity, std::string const & filename) :
    severity{_severity}, sink{new std::ofstream{filename, std::ios::binary}, stream_deleter_default}
  {
  }

  // The object's state
  log_severity severity;
  stream_ptr_t sink;
  std::mutex mutex;
};

// GLOBAL INSTANCE OF LOGGING STATE
inline std::unique_ptr<log_singleton_t> log_singleton = nullptr;

// THE PRINT FUNCTION
template <typename... args_t>
void print_log(log_severity const severity, args_t &&... args)
{
  assert(log_singleton != nullptr);
  assert(log_singleton->sink != nullptr);

  if (severity < log_singleton->severity)
    return;

  std::lock_guard guard{log_singleton->mutex};

  // time
  auto now = std::chrono::system_clock::now();
  auto in_time_t = std::chrono::system_clock::to_time_t(now);
  *log_singleton->sink << std::put_time(std::localtime(&in_time_t), "[%Y-%m-%d %H:%M:%S");

  auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
  *log_singleton->sink << '.' << std::setfill('0') << std::setw(3) << milliseconds.count() % 1000 << ']';

  // severity
  *log_singleton->sink << " <" << severity2string(severity) << "> ";

  // args
  ((*log_singleton->sink << args), ...);

  *log_singleton->sink << '\n';
}

} // namespace gyper
