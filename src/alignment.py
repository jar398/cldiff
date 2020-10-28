import sys

import checklist as cl
import relation as rel
import articulation as art
import intension
import dribble
from intension import choose_best_match

# For each B-record, we choose an articulation with the closest
# A-record that it matches (preferably but not necessarily an '='
# articulation).

def align(B, A):

  # Precompute all best matches
  best = intension.best_intensional_match_map(B, A)

  # Turn tipward best matches into = or < articulations as appropriate
  tipwards = intension.intensional_alignment(tipward(best, A, B))

  # Extensional analysis
  cross_mrcas = analyze_cross_mrcas(B, A, tipwards)
  dribble.log("# Number of cross-mrcas: %s" % len(cross_mrcas))
  ext_map = extensional_match_map(A, B, tipwards, cross_mrcas)
  dribble.log("# Number of extensional relationships: %s" % len(ext_map))

  # Add extensional matches to a draft that already has intensional matches
  the_alignment = assemble_alignment(tipwards, best, ext_map)
  return (the_alignment, cross_mrcas)

# Side-affects draft

def assemble_alignment(draft, best, ext_map):
  for node in ext_map:
    alignment_step(node, best, ext_map, draft)
  return draft

def alignment_step(node, best, ext_map, draft):
  def luup(x, y):
    if x == cl.forest_tnu or y == cl.forest_tnu: return
    assert cl.get_checklist(x) != cl.get_checklist(y)
    if in_chain(x, y0) and in_chain(y, x0):
      bar = best.get(x)
      if bar and find_in_chain(bar.cod, y, x0):
        if bar.cod == y:
          art.proclaim(draft, art.change_relation(bar, rel.eq, "extensional"))
          luup(cl.get_parent(x), cl.get_parent(y))
        else:
          art.proclaim(draft, art.extensional(y, x, rel.lt, "presumed 0"))
          luup(x, cl.get_parent(y))
      else:
        bar = best.get(y)
        if bar and find_in_chain(bar.cod, x, y0):
          if bar.cod == y:
            art.proclaim(draft, bar)
            luup(cl.get_parent(x), cl.get_parent(y))
          else:
            art.proclaim(draft, art.extensional(x, y, rel.lt, "presumed 1"))
            luup(cl.get_parent(x), y)
        else:
          # neither x nor y matches by name
          art.proclaim(draft, art.extensional(x, y, rel.eq, "presumed mutual"))
          luup(cl.get_parent(x), cl.get_parent(y))

  # See is b is in the chain (matching nodes in lineage)
  def find_in_chain(b, y, x0):
    if y == cl.forest_tnu:
      return False
    elif in_chain(y, x0):
      if b == y:
        return True
      else:
        return find_in_chain(b, cl.get_parent(y), x0)
    else:
      return False

  def in_chain(y, x0):
    if y in draft: return
    ar = ext_map.get(y)
    return ar and ar.cod == x0 and ar.relation == rel.matches

  if not draft.get(node):
    ar = ext_map.get(node)
    if ar:
      x0 = ar.dom
      y0 = ar.cod
      if in_chain(x0, y0) and in_chain(y0, x0):
        if x0 < y0:
          luup(x0, y0)
        else:
          pass    # pick it up later (or earlier)
      else:
        art.proclaim(draft, ar)

# ---------- Extensionality by particle set

# Suppress < transitivity

def extensional_match_map(A, B, draft, xmrcas):
  ext = {}
  def process(node, less):
    if not draft.get(node):
      e = extensional_match(node, xmrcas)
      if e:
        if e.relation == rel.lt:
          if less and e.cod == less.cod:
            print("# Suppressing %s because redundant with\n  %s" %
                  (art.express(e), art.express(less)))
            pass
          else:
            less = e
            ext[node] = e
        else:
          less = None
          ext[node] = e
      else:
        less = None
    for child in cl.get_children(node):
      process(child, less)
  for root in cl.get_roots(A): process(root, None)
  for root in cl.get_roots(B): process(root, None)
  return ext

# Guaranteed invertible, except for monotypic node chains
# This code is derived from 
#   reference-taxonomy/org/opentreeoflife/conflict/ConflictAnalysis.java

def extensional_match(node, xmrcas):
  partner = xmrcas.get(node)      # node in other checklist; 'conode'
  if not partner:
    # Descendant of a particle
    if dribble.watch(node):
      dribble.log("# EM: %s is not tipward." % cl.get_unique(node))
    return None
  back = xmrcas.get(partner)    # 'bounce'
  if not back:
    # Not sure how this can happen but it does (NCBI vs. GBIF)
    dribble.log("%s <= %s <= nowhere" % (cl.get_unique(node),
                                         cl.get_unique(partner)))
    if dribble.watch(node):
      dribble.log("# EM: %s killed because aborted round trip." % cl.get_unique(node))
    return None
  # node <= partner <= back
  how = cl.how_related(node, back)    # Alway rcc5
  if how == rel.eq:
    # Should end up being eq iff name match or unique match
    # Can test for unique match by looking at xmrca of parent

    # Could be part of a 'monotypic' chain; fix later
    how = rel.matches
    reason = "mutual-cross-mrca"
  elif how == rel.gt:
    how = rel.matches
    reason = "monotypic-inversion"
  elif how == rel.disjoint:
    reason = "particle-set-exclusion"
  else:               # must be rel.lt
    # Assume resolution (node < partner) until conflict is proven
    reason = "refinement"
    # Look for an intersection between any partner-child and node
    # x is in A checklist, y is in B checklist
    for pchild in cl.get_children(partner):
      pchild_back = xmrcas.get(pchild)
      if pchild_back == None:
        # pchild ! node
        pass
      else:
        (d, e) = cross_compare(node, pchild, xmrcas)
        # d < node while e ! node
        if d and e:
          how = rel.conflict
          reason = ("%s is in; its sibling %s is not" %
                    (cl.get_unique(d), cl.get_unique(e)))
          dribble.log("** %s (x) conflicts with %s (y) because\n"
                      "   %s in x isn't in y (its sibling %s is)" %
                      (cl.get_unique(node),
                       cl.get_unique(partner),
                       cl.get_unique(e),
                       cl.get_unique(d)))
          break
        elif e:
          reason = ("%s is not in it" % cl.get_unique(e))

  ar = art.extensional(node, partner, how, reason)
  if dribble.watch(node):
    dribble.log("# Extensional articulation %s" % art.express(ar))
  return ar

# ---------- Determine disjointness across checklists

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

def cross_compare(node, conode, xmrcas):
  back = xmrcas.get(conode)
  if back == None:
    return (None, None)
  how = cl.how_related(node, back)
  if how == rel.eq:
    return (conode, None)    # conode <= node
  elif how == rel.gt:
    return (conode, None)    # conode <= node
  elif how == rel.disjoint:
    return (None, conode)    # conode ! node
  assert how == rel.lt    # node < back is inconclusive
  x_seen = None
  y_seen = None
  for child in cl.get_children(conode):
    saw = cross_compare(node, child, xmrcas)
    if saw:
      (x, y) = saw
      if x and y:
        return saw
      elif x:
        if y_seen:
          return (x, y_seen)
        else:
          x_seen = x
      elif y:
        if x_seen:
          return (x_seen, y)
        else:
          y_seen = y
  return (x_seen, y_seen)

# ---------- Cross-MRCAs

def analyze_cross_mrcas(A, B, tipwards):
  cross_mrcas = {}
  def half_analyze_cross_mrcas(checklist, other):
    def subanalyze_cross_mrcas(node, other):
      result = None
      probe = tipwards.get(node)
      if probe:
        # Could be: = < or >
        result = probe.cod
      else:
        children = cl.get_children(node)
        if children:
          m = None      # None is the identity for mrca
          for child in children:
            m2 = subanalyze_cross_mrcas(child, other)
            if m2 != None:
              m = cl.mrca(m, m2) if m != None else m2
          if m != None:
            result = m
      if result:
        assert cl.get_checklist(result) != cl.get_checklist(node)
        if dribble.watch(node):
          dribble.log("# Cross-mrca(%s) = %s" %
                      (cl.get_unique(node), cl.get_unique(result)))
        cross_mrcas[node] = result
      return result             # in B
    for root in cl.get_roots(checklist):
      subanalyze_cross_mrcas(root, other)
  half_analyze_cross_mrcas(A, B)
  half_analyze_cross_mrcas(B, A)

  # Sanity check
  for node in cross_mrcas:
    cross = cross_mrcas[node]
    probe = cross_mrcas.get(cross)
    if probe:
      assert cl.get_checklist(probe) == cl.get_checklist(node)
    else:
      dribble.log("# No return cross-MRCA for %s -> %s -> ..." %\
                  (cl.get_unique(node), cl.get_unique(cross)))
  return cross_mrcas

# ---------- Tipwards

# A particle is a mutual =-articulation of tipward accepted nodes
# deriving from sameness of 'intrinsic' node properties: name, id,
# rank, parent, etc.

# This function returns a partial map from nodes to articulations.

# Filter out internal nodes (those having a matched descendant)

def tipward(amap, A, B):
  tw = {}
  def filter(node):
    debug = dribble.watch(node)
    found_match = None
    for child in cl.get_children(node):
      ar = filter(child)
      if ar:
        found_match = ar
    if found_match:    # Some descendant is a particle
      if debug: dribble.log("# %s: descendant matches, not keeping: %s" %
                            (cl.get_unique(node), art.express(found_match)))
      return found_match
    elif node in amap:
      ar = amap[node]
      tw[ar.dom] = ar
      if debug: dribble.log("# %s is a tipward match, keeping: %s" %
                            (cl.get_unique(node), art.express(ar)))
      return ar
    else:
      if debug: dribble.log("# %s is unmatched" % cl.get_unique(node))
      return None
  for root in cl.get_roots(A): filter(root)
  for root in cl.get_roots(B): filter(root)
  return tw

