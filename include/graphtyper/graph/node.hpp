#pragma once
#include <iostream>
#include <vector>

#include <boost/serialization/access.hpp>

#include <graphtyper/graph/label.hpp>


namespace gyper
{

class RefNode
{
  friend class boost::serialization::access;

public:
  Label label{};
  std::vector<TNodeIndex> out_var_ids{};

  RefNode();
  RefNode(Label && l, std::vector<TNodeIndex> && c) noexcept;
  RefNode(RefNode const & rn) noexcept;

  void change_label_order(uint32_t change);

  std::size_t out_degree() const;
  Label const & get_label() const;
  TNodeIndex get_var_index(unsigned const & index) const;
  std::vector<TNodeIndex> get_vars() const;

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

  VarNode();
  VarNode(Label && l, TNodeIndex && ori) noexcept;
  VarNode(VarNode const & vn) noexcept;

  void change_label_order(uint32_t change);

  std::size_t out_degree() const;
  Label const & get_label() const;
  TNodeIndex get_out_ref_index() const;

private:
  template <class Archive>
  void serialize(Archive & ar, const unsigned int);
};


} // namespace gyper
