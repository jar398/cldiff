# Given two source checklists A and B and an alignment between them,
# compute the merge of the two checklists in the following sense:

# 1. The set of 'merge nodes' form a forest, i.e. every merge node has a
#    merge node as a parent, if it's not a root.
# 2. Every node of the merge comes either from a mutual match between the 
#    sources, or from an unmatched node in either source.
# 3. The second or 'B' checklist has priority in the sense that nodes in the first
#    that conflict with the second are not included in the merge.

# The merge tree is represented by the node-to-parent mapping.  The inverse
# of this mapping (children) is computed on demand.

import checklist as cl
import relation as rel

def merge_checklists(A, B, al, xmrcas):
  parents = {}
  roots = []
  def half_compute_parents(check, inject, al):
    def process(node):
      for child in cl.get_children(node):
        process(child)
      merged = inject(node, al)
      if not merged in parents:
        p = merged_parent(merged, al)
        if p:
          parents[merged] = p     # Otherwise it's a root
        else:
          if not merged in roots:
            roots.append(merged)
    for root in cl.get_roots(check):
      process(root)
  half_compute_parents(A, inject_A, al)
  half_compute_parents(B, inject_B, al)
  return (parents, roots)

# A is low priority

def merged_parent(merged, al):
  (x, y) = merged    # False if node is inconsistent

  if x and y:                   # x ≈ y
    assert al[x].cod == y
    p = cl.get_parent(x)
    q = cl.get_parent(y)
    ar = get_consistent_articulation(p, al)
    if ar:
      if cl.how_related(ar.cod, q) == rel.lt:
        return inject_A(ar.dom, al)
      else:
        return inject_B(q, al)
    else:
      return inject_A(p, al)

  elif x:
    ar = get_consistent_articulation(x, al)
    if ar:
      q = ar.cod
      return inject_B(q, al)
    else:
      return inject_A(cl.get_parent(x), al)

  elif y:
    return inject_B(cl.get_parent(y), al)

  else:
    assert False

# Scan up through lineage to find non-conflict articulation

def get_consistent_articulation(node, al):
  ar = al.get(node)
  while ar and rel.is_variant(ar.relation, rel.conflict):
    node = cl.get_parent(node)
    ar = al.get(node)
  return ar

# Inject node in low priority checklist into merged checklist

def inject_A(node, al):
  if node == cl.forest_tnu: return None
  m1 = al.get(node)
  if m1:
    if is_eq(m1.relation):
      m2 = al.get(m1.cod)
      if m2 and is_eq(m2.relation):
        if m2.cod == node:
          return (node, m1.cod)
  return (node, None)

# Inject node in high priority checklist into merged checklist

def inject_B(node, al):
  if node == cl.forest_tnu: return None
  m1 = al.get(node)
  if m1:
    if is_eq(m1.relation):
      m2 = al.get(m1.cod)
      if m2 and is_eq(m2.relation):
        if m2.cod == node:
          return (m1.cod, node)
  return (None, node)

def is_eq(re):
  return rel.is_variant(re, rel.eq)  # == rel.intensional
