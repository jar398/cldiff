# Given two source checklists and an alignment between them,
# compute the merge of the two checklists in the following sense:

# 1. Every node of the merge comes either from a mutual match between the 
#    sources, or from an unmatched node in either source.
# 2. The parent of a merge node is a merge node.
# 3. The second checklist has priority in the sense that nodes in the first
#    that conflict with the second are not included in the merge.

# The merged tree is represented by the parent mapping.  The inverse
# (children) is computed on demand.

import checklist as cl
import relation as rel

def merge_checklists(A, B, al, xmrcas):
  parents = {}
  roots = []
  def half_compute_parents(check, inject, al):
    def process(node):
      merged = inject(node, al)
      if not merged in parents:
        p = merged_parent(merged, inject, al, xmrcas)
        if p:
          parents[merged] = p     # Otherwise it's a root
        else:
          if not merged in roots:
            roots.append(merged)
      for child in cl.get_children(node):
        process(child)
    for root in cl.get_roots(check):
      process(root)
  half_compute_parents(A, inject_A, al)
  half_compute_parents(B, inject_B, al)
  return (parents, roots)

def merged_parent(merged, inject, al, xmrcas):
  (x, y) = merged    # False if node is inconsistent
  if x:
    u = consistent_parent(x, al)
  else:
    u = None
  if y:
    v = cl.get_parent(y)
  else:
    v = None
  if not u:
    if not v:
      # node is a root in both checklists.  no merged parent
      return None
    # node is only a root in high priority checklist
    return inject_B(v, al)
  elif not v:
    # node is only a root in low priority checklist
    return inject_A(u, al)
  else:
    over = xmrcas.get(u)
    if over and cl.mrca(over, v) == v:
      # parent is in low priority checklist - "insertion"
      return inject_A(u, al)
  return inject_B(v, al)

# Skip ancestors that are inconsistent with high priority checklist

def consistent_parent(x, al):
  u = cl.get_parent(x)
  u_art = al.get(u)
  if not u_art:
    return u
  elif not rel.is_variant(u_art.relation, rel.conflict):
    return u
  else:
    # u conflicts with checklist B.  Continue up lineage
    return consistent_parent(u, al)

# Inject node in low priority checklist into merged checklist

def inject_A(node, al):
  m1 = al.get(node)
  if m1:
    if rel.is_variant(m1.relation, rel.eq):
      m2 = al.get(m1.cod)
      if m2 and rel.is_variant(m2.relation, rel.eq):
        if m2.cod == node:
          return (node, m1.cod)
  return (node, None)

# Inject node in high priority checklist into merged checklist

def inject_B(node, al):
  m1 = al.get(node)
  if m1:
    if rel.is_variant(m1.relation, rel.eq):
      m2 = al.get(m1.cod)
      if m2 and rel.is_variant(m2.relation, rel.eq):
        if m2.cod == node:
          return (m1.cod, node)
  return (None, node)
