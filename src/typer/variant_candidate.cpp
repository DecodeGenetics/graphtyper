#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include <graphtyper/utilities/logging.hpp>

#include <graphtyper/typer/variant_candidate.hpp>


namespace gyper
{

bool
VariantCandidate::add_base_in_front(bool const add_N)
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = std::move(seqs);
  bool ret = new_var.add_base_in_front(add_N);
  abs_pos = new_var.abs_pos;
  seqs = std::move(new_var.seqs);
  return ret;
}


bool
VariantCandidate::add_base_in_back(bool const add_N)
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = std::move(seqs);
  bool ret = new_var.add_base_in_back(add_N);
  abs_pos = new_var.abs_pos;
  seqs = std::move(new_var.seqs);
  return ret;
}


void
VariantCandidate::expanded_normalized()
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = std::move(seqs);
  new_var.expanded_normalized();
  abs_pos = new_var.abs_pos;
  seqs = std::move(new_var.seqs);
}


void
VariantCandidate::normalize()
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = std::move(seqs);

  //std::cerr << new_var.abs_pos << " " << new_var.seqs.size()
  //          << " " << std::string(new_var.seqs[0].begin(), new_var.seqs[0].end())
  //          << " " << std::string(new_var.seqs[1].begin(), new_var.seqs[1].end()) << "\n";

  new_var.normalize();
  //std::cerr << "ok\n";
  abs_pos = new_var.abs_pos;
  seqs = std::move(new_var.seqs);
}


bool
VariantCandidate::is_normalized() const
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = seqs;
  return new_var.is_normalized();
}


bool
VariantCandidate::is_snp_or_snps() const
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = seqs;
  return new_var.is_snp_or_snps();
}


int
VariantCandidate::is_transition_or_transversion() const
{
  if (seqs.size() == 2 && seqs[0].size() == 1 && seqs[1].size() == 1)
  {
    if ((seqs[0][0] == 'A' && seqs[1][0] == 'G') ||
        (seqs[0][0] == 'G' && seqs[1][0] == 'A') ||
        (seqs[0][0] == 'C' && seqs[1][0] == 'T') ||
        (seqs[0][0] == 'T' && seqs[1][0] == 'C'))
    {
      return 1; // transition
    }
    else
    {
      return 2; // transversion
    }
  }

  return false;
}


std::string
VariantCandidate::print() const
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = seqs;
  return new_var.to_string();
}


bool
VariantCandidate::operator==(VariantCandidate const & v) const
{
  return abs_pos == v.abs_pos && seqs == v.seqs;
}


bool
VariantCandidate::operator!=(VariantCandidate const & b) const
{
  return !(*this == b);
}


bool
VariantCandidate::operator<(VariantCandidate const & b) const
{
  return abs_pos < b.abs_pos || (abs_pos == b.abs_pos && seqs < b.seqs);
}


std::size_t
VariantCandidateHash::operator()(VariantCandidate const & v) const
{
  assert(v.seqs.size() == 2);
  std::size_t h1 = std::hash<uint32_t>()(v.abs_pos);
  std::size_t h2 = std::hash<std::string_view>{}(std::string_view{v.seqs[0].data(), v.seqs[0].size()});
  std::size_t h3 = 42 + std::hash<std::string_view>{}(std::string_view{v.seqs[1].data(), v.seqs[1].size()});
  return h1 ^ (h2 << 1) ^ (h3 + 0x9e3779b9);
}


} // namespace gyper
