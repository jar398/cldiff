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
import dribble

def merge_checklists(A, B, al):
  parents = {}
  roots = []
  def half_compute_parents(check, inject, al):
    def process(node):
      merged = inject(node, al)
      if not merged in parents:
        p = merged_parent(merged, al)
        if p:
          if dribble.watch(node):
            (x, y) = p
            dribble.log("# Merged parent(%s) = (%s, %s)" %
                        (cl.get_unique(node), cl.get_unique(x), cl.get_unique(y)))
          parents[merged] = p     # Otherwise it's a root
        else:
          if dribble.watch(node):
            dribble.log("# No merge(%s)" % cl.get_unique(node))
          if not merged in roots:
            roots.append(merged)
      for child in cl.get_children(node):
        process(child)
    for root in cl.get_roots(check):
      process(root)
  half_compute_parents(B, inject_B, al)
  half_compute_parents(A, inject_A, al)    # these will not override
  return (parents, roots)

# A is low priority

def merged_parent(merged, al):
  (x, y) = merged    # False if node is inconsistent

  if y:
    q = cl.get_parent(y)
    if x:
      p = cl.get_parent(x)
    else:
      p = partner(y, al)    # cannot be =, must be <

    # If p takes us back to q, then parent should be p, else q
    scan = p
    while scan and scan != cl.forest_tnu:
      z = partner(scan, al)
      if z:
        if z == q:
          return inject_A(p, al)
        else:
          return inject_B(q, al)
      scan = cl.get_parent(scan)
    return inject_B(q, al)

  else:
    assert x
    q = partner(x, al)
    if q:
      return inject_B(q, al)
    else:
      return inject_A(cl.get_parent(x), al)

# Inject node in A (low priority checklist) into merged checklist

def inject_A(node, al):
  if node == cl.forest_tnu: return None
  return (node, eq_partner(node, al))

# Inject node in B (high priority checklist) into merged checklist

def inject_B(node, al):
  if node == cl.forest_tnu: return None
  return (eq_partner(node, al), node)

# Aligned node, if relation is < or =

def partner(x, al):
  ar = al.get(x)
  if ar and (ar.relation == rel.eq or ar.relation == rel.refines):
    return ar.cod
  else:
    return None

# Aligned node, if relation is =

def eq_partner(node, al):
  ar = al.get(node)
  if ar and ar.relation == rel.eq:
    return ar.cod
  else:
    return None
