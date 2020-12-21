#pragma once
#include <iostream>
#include <unordered_set>
#include <vector>

#include <boost/serialization/access.hpp>

#include <graphtyper/graph/label.hpp>


namespace gyper
{

class Ref;
class Alt;

class RefNode
{
  friend class boost::serialization::access;

public:
  Label label{};
  std::vector<TNodeIndex> out_var_ids{};

  RefNode() = default;
  RefNode(Label && l, std::vector<TNodeIndex> && c) noexcept;
  RefNode(RefNode const & rn) noexcept;
  RefNode(RefNode && o) noexcept;
  ~RefNode() = default;

  std::size_t out_degree() const;
  Label const & get_label() const;
  TNodeIndex get_var_index(unsigned const & index) const;
  std::vector<TNodeIndex> const & get_vars() const;

private:
  template <class Archive>
  void serialize(Archive & ar, const unsigned int);
};


class VarNode
{
  friend class boost::serialization::access;

public:
  Label label{};
  TNodeIndex out_ref_id{};
  std::unordered_set<long> events;
  std::unordered_set<long> anti_events;

  VarNode() = default;
  //VarNode(Label && l, TNodeIndex && ori, std::unordered_set<long> && events) noexcept;

  VarNode(Label && l, TNodeIndex && ori,
          std::unordered_set<long> && events,
          std::unordered_set<long> && anti_events) noexcept;

  VarNode(VarNode const & vn) noexcept;
  VarNode(VarNode && o) noexcept;
  ~VarNode() = default;

  std::size_t out_degree() const;
  Label const & get_label() const;
  TNodeIndex get_out_ref_index() const;

private:
  template <class Archive>
  void serialize(Archive & ar, const unsigned int);
};


} // namespace gyper
