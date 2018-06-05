#include <boost/log/trivial.hpp>
#include <boost/functional/hash.hpp> // boost::hash_range

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
VariantCandidate::normalize()
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = std::move(seqs);
  new_var.normalize();
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


std::string
VariantCandidate::print() const
{
  Variant new_var;
  new_var.abs_pos = abs_pos;
  new_var.seqs = seqs;
  return new_var.print();
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
  std::size_t h2 = boost::hash_range(v.seqs[0].begin(), v.seqs[0].end());
  std::size_t h3 = boost::hash_range(v.seqs[1].begin(), v.seqs[1].end());
  return h1 ^ (h2 << 1) ^ (h3 + 0x9e3779b9);
}

} // namespace gyper
