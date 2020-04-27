#include <algorithm> // std::generate_n
#include <cassert> // assert
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string> // std::string
#include <sstream>
#include <utility>
#include <vector>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include <graphtyper/utilities/system.hpp>

#include <boost/log/trivial.hpp>

namespace
{

std::string
get_env_var(std::string const & key, std::string const & _default = "")
{
  char * val = getenv(key.c_str());
  return val ? val : _default;
}


std::string
get_random_string(long length)
{
  // from https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
  char const charset[] = "0123456789"
                         "ABCDEFGHIJ"
                         "KLMNOPQRST"
                         "UVWXYZ"
                         "abcdefghijklmnopqrstuvwxyz";

  std::vector<char> char_vec;
  char_vec.reserve(sizeof(charset) - 1); // skip the last NULL
  std::copy(charset, charset + sizeof(charset) - 1, std::back_inserter(char_vec));
  assert(char_vec.size() > 0);

  std::default_random_engine rng(std::random_device{} ());
  std::uniform_int_distribution<> dist(0, char_vec.size() - 1);

  std::string str(length, 0);
  std::generate_n(str.begin(), length, [&](){
      return char_vec[dist(rng)];
    });

  assert(std::find(str.begin(), str.end(), '\0') == str.end()); // No NULLs in final string
  return str;
}


} // namespace anon


namespace gyper
{

void
create_dir(std::string const & dir, unsigned mode)
{
  mkdir(dir.c_str(), mode);
}


std::string
current_sec()
{
  time_t now = time(0);
  struct tm time_structure;
  char buf[32];
  time_structure = *localtime(&now);
  strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &time_structure);
  std::string const sec(buf);
  return sec;
}


std::string
create_temp_dir(GenomicRegion const & region)
{
  std::ostringstream ss;
  ss << get_env_var("TMPDIR", "/tmp") << "/graphtyper_" << current_sec() << "_"
     << region.chr << "_" << std::setw(9) << std::setfill('0') << (region.begin + 1)
     << "." << get_random_string(6);

  std::string tmp = ss.str();
  create_dir(tmp, 0700);
  return tmp;
}


void
remove_file_tree(std::string const & path)
{
  std::string full_path;
  DIR * dir;
  struct stat stat_path, stat_entry;
  struct dirent * entry;

  // stat for the path
  stat(path.c_str(), &stat_path);

  // if path does not exists or is not dir - exit with status -1
  if (S_ISDIR(stat_path.st_mode) == 0)
  {
    BOOST_LOG_TRIVIAL(error) << "'" << path << "' is not directory.";
    std::exit(1);
  }

  // if not possible to read the directory for this user
  if ((dir = opendir(path.c_str())) == NULL)
  {
    BOOST_LOG_TRIVIAL(error) << "Can't open directory '" << path << "'";
    std::exit(1);
  }

  // iteration through entries in the directory
  while ((entry = readdir(dir)) != NULL)
  {
    // skip entries "." and ".."
    if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, ".."))
      continue;

    // determinate a full path of an entry
    std::string full_path = path + "/" + entry->d_name;

    // stat for the entry
    stat(full_path.c_str(), &stat_entry);

    // recursively remove a nested directory
    if (S_ISDIR(stat_entry.st_mode) != 0)
    {
      remove_file_tree(full_path);
      continue;
    }

    // remove a file object
    if (unlink(full_path.c_str()) == 0)
    {
      //printf("Removed a file: %s\n", full_path);
    }
    else
    {
      BOOST_LOG_TRIVIAL(warning) << "Can't remove a file '" << full_path << "'";
    }
  }

  // remove the devastated directory and close the object of it
  if (rmdir(path.c_str()) == 0)
  {
    //printf("Removed a directory: %s\n", path);
  }
  else
  {
    BOOST_LOG_TRIVIAL(warning) << "Can't remove a directory '" << path << "'";
  }

  closedir(dir);
}


bool
is_file(std::string const & filename)
{
  struct stat sb;
  return stat(filename.c_str(), &sb) == 0 && (S_ISREG(sb.st_mode) || S_ISLNK(sb.st_mode));
}


void
check_file_exists(std::string const & filename)
{
  if (filename.size() == 0 || !is_file(filename))
  {
    BOOST_LOG_TRIVIAL(error) << "No file '" << filename << "' exists.";
    std::exit(1);
  }
}


void
check_file_exists_or_empty(std::string const & filename)
{
  if (filename.size() > 0 && !is_file(filename))
  {
    BOOST_LOG_TRIVIAL(error) << "No file '" << filename << "' exists.";
    std::exit(1);
  }
}


bool
is_directory(std::string const & filename)
{
  struct stat sb;
  return stat(filename.c_str(), &sb) == 0;
}


bool
is_defined_in_env(std::string const & var)
{
  char * ref_cache;
  ref_cache = getenv(var.c_str());

  if (ref_cache != NULL)
    return true;

  return false;
}


} // namespace gyper
