#pragma once

// Source: https://askubuntu.com/questions/1256440/how-to-get-libstdc-with-c17-filesystem-headers-on-ubuntu-18-bionic
#if __has_include(<filesystem>)
#  include <filesystem>
namespace gyper
{
namespace filesystem = std::filesystem;
}
#elif __has_include(<experimental/filesystem>)
#  include <experimental/filesystem>
namespace gyper
{
namespace filesystem = std::experimental::filesystem;
}
#else
#  error "Missing the <filesystem> header."
#endif
