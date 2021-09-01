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

# 'al' is a proposal (alignment)

def merge_checklists(A, B, al):
  # Retractions.
  retractions = find_incompatibilities(A, B, al)
  parents = {}
  roots = []
  def half_compute_parents(check, inject, al, retractions):
    def process(node):
      merged = inject(node, al)
      if merged in parents:
        # how can this happen?
        pass
      else:
        p = merged_parent(merged, al, retractions)
        if p:
          if dribble.watch(node):
            (x, y) = p
            dribble.log("# Merged parent(%s) = (%s, %s)" %
                        (cl.get_unique(node), cl.get_unique(x), cl.get_unique(y)))
          if merged in p:
            def w(m):
              (x,y)=m
              return ("%s/%s" % (cl.get_unique(x), cl.get_unique(y)))
            dribble.log("** Setting merged parent twice?\n  %s -> %s then %s" %
                        (w(merged), w(p), w(parents[merged])))
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
  half_compute_parents(B, inject_B, al, {})
  half_compute_parents(A, inject_A, al, retractions)    # these will not override
  return (parents, roots)

# A is low priority

def merged_parent(merged, al, retract):

  # Like get_parent except skip retracted nodes
  def compatible_ancestor(x):
    p = cl.get_parent(x)
    if p in retract:
      return compatible_ancestor(p)
    return x

  (x, y) = merged    # False if node is inconsistent

  # Suppose x = y.  There are four candidates for merged parent:
  #   p (in A) = x's parent
  #   n (in B) = x's partner
  #   q (in B) = y's parent
  #   m (in A) = y's partner
  # The goal is to pick the smallest one of the four.
  # x < p <= n, y < q <= m, either p = m or q = n (??)

  if y:
    if x:
      p = compatible_ancestor(x)
      q = compatible_ancestor(y)
      if partner(y, al) == p: # m
        return inject_B(q, al)
      # If p's partner is x's partner, then ... don't use x's partner
      if partner(x, al) == q: # n
        return inject_A(p, al)
      else:
        assert False
    else:
      return inject_A(partner(y, al), al)
  elif x:
    return inject_B(partner(x, al), al)
  else:
    assert False

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
  if ar and (ar.relation == rel.eq or ar.relation == rel.lt):
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

# "Incompatibility" check - is node x (in A) compatible with checklist B?
# N.b. proposal only contains = and < articulations

def find_incompatibilities(A, B, proposal):
  retractions = {}
  def process(x):
    ar = proposal.get(x)
    if ar and ar.relation == rel.lt:
      proof = test_compatibility(ar.dom, ar.cod, proposal)
      if proof:
        retractions[x] = proof
    for child in cl.get_children(x):
      process(child)
  for root in cl.get_roots(A):
    process(root)
  return retractions

# (yk, c, d, e) = test_compatibility(x, y, xmrcas)
# y = xmrcas(xmrcas(x))

def test_compatibility(x, y, proposal):
  # Look for an intersection between any y-child and x
  # x is in A checklist, y is in B checklist
  # Let y1, y2, ... yn be y's children.
  # We seek a child yk with descendants (d,e) with d<yk and and e!yk,
  # and a second child yj with c<yj in x (i.e. c and x intersect).
  yk = c = d = e = None
  for yi in cl.get_children(y):
    if not yk:
      (di, ei) = cross_compare(x, yi, proposal)
      if di and ei:
        yk = yi; d = di; e = ei
      elif not di and not ei:
        c = yi
    else:
      c = yi
  if yk:
    return (yk, c, d, e)
  else:
    return None

def express_proof(proof):
  (c, d, e) = proof
  # Assume resolution (x < y) until conflict is proven
      # assume potential child until proven otherwise
  if c and d and e:
    proof_expression = (">< %s [%s, %s, %s]" % 
                        (cl.get_unique(yk), cl.get_unique(c),
                         cl.get_unique(d), cl.get_unique(e)))
    dribble.log("** %s doesn't refine %s because\n   %s\n   yk [in x, in both, in yk]" %
                (cl.get_unique(x), cl.get_unique(y), proof_expression))
  # Should squirrel away the proof somewhere!
  return proof_expression


# ---------- Determine disjointness across checklists
# Doesn't seem to belong in this file.

# Compare node (in A) to conode (in B) according to particle sets.

# Returns disjoint descendants (x, y) of conode where
#    x <= conode ∩ node, or x = None if conode ! node
#    y <= conode - node, or y = None if conode <= node

# We're looking for two things
#  1. whether they intersect (if so then parent conflict can be detected) -
#      i.e. a conode descendant that *is* a node descendant.
#  2. why they are not equal (goes to conflict diagnosis) - example of
#      conode descendant that is *not* a node descendant.

# On each recursive call, we can obtain one or the other of these, or both.
# Or maybe none.

# The recursion is over conode, searching for sibling descendants for
# which one is in node and the other is not.
def cross_compare(node, conode, proposal):
  ar = proposal.get(conode) # y -> x'
  if ar == None:
    return (None, None)    
  how = cl.how_related(node, ar.cod)
  if how == rel.disjoint:
    return (None, conode)    # conode ! node
  elif how == rel.eq or how == rel.gt:
    return (conode, None)    # conode <= node
  assert how == rel.lt    # node < back is inconclusive
  d_seen = None           # < both node and conode
  e_seen = None           # < conode only
  for child in cl.get_children(conode):
    saw = cross_compare(node, child, proposal)
    if saw:
      (d, e) = saw
      if d and e:
        return saw
      elif d:
        if e_seen:
          return (d, e_seen)
        else:
          d_seen = d
      elif e:
        if d_seen:
          return (d_seen, e)
        else:
          e_seen = e
  return (d_seen, e_seen)

