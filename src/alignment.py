debug = False

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

def align(B, A, captured = {}):
  intension.clear_cache()
  best = intension.best_intensional_match_map(B, A, captured)
  tipward_best = tipward(best, B, A)
  particles = mutual(tipward_best)
  check_mutuality(particles)
  dribble.log("# Number of particles: %s" % (len(particles)>>1))

  cross_mrcas = analyze_cross_mrcas(B, A, particles)
  dribble.log("# Number of cross-mrcas: %s" % len(cross_mrcas))

  # Find membership-based relationships
  extensions = analyze_extensions(cross_mrcas)

  draft = fix_alignment(particles)    # change ≈ to =
      # should be derived from best tips, not particles

  the_alignment = finish_alignment(B, A, draft, cross_mrcas)
  dribble.log("# Number of articulations in alignment: %s" % len(the_alignment))

  return (the_alignment, cross_mrcas)

# This is a godawful mess

def analyze_extensions(cross_mrcas):
  return {}

def finish_alignment(B, A, particles, xmrcas):
  global cross_mrcas
  cross_mrcas = xmrcas
  draft = {}
  def half_finish(B, A):
    def process(node):
      m = articulate(node, A, particles, xmrcas)
      if m:
        assert cl.is_accepted(m.cod)
        draft[node] = m
      else:
        # Otherwise, what?
        pass
      for child in cl.get_children(node):
        process(child)
    for root in cl.get_roots(B):
      process(root)
  half_finish(B, A)
  half_finish(A, B)
  return draft

# Deal with ad hoc splits and merges.
# This currently doesn't do much of anything.

def fix_alignment(draft):
  alignment = {}

  # Suppose `node` x comes from the A checklist, and there is a split
  # such that x matches multiple nodes y1, y2 in the B checklist.

  # For each A-node x, find all B-nodes y, ... that are contending to
  # be equivalent to x.

  incoming = {}
  for node in draft:    # y node
    ar = draft[node]     # ar : y -> x
    if ar.cod in incoming:
      incoming[ar.cod].append(ar) # x -> (y->x, y2->x)
    else:
      incoming[ar.cod] = [ar]     # x -> (y->x)

  # Modify the relation for all approximate-match nodes.

  for y in incoming:          # y node
    # e.g. inc = {a ≈ y, b ≈ y}
    inc = incoming[y]    # cod = y for all articulations
    # mut: y ≈ x
    mut = draft.get(y)   # mut : y -> mutual
    sibs = []
    if mut:
      x0 = mut.cod
      rent = cl.get_parent(x0)    # parent of x
      for ar in inc:    # ar: x ≈ y
        if (rel.is_variant(ar.relation, rel.eq) and
            # parent of b is same as parent of a?
            cl.get_parent(ar.dom) == rent):
          sibs.append(ar)
    for ar in inc:
      x = ar.dom
      if ar in sibs:
        if len(sibs) == 1:
          alignment[x] = art.set_relation(ar, rel.eq) # x = y
        else:
          alignment[x] = art.change_relation(ar, rel.merge, "merge", "split")
      else:
        alignment[x] = ar       # x < y
    if len(sibs) > 1:
      # Report!
      dribble.log("# Split/join %s -> %s -> %s" %
                  (" ∨ ".join(map(lambda e:cl.get_unique(e.dom), sibs)),
                   cl.get_unique(y),
                   cl.get_unique(rent)))
      
  return alignment

# ---------- One-sided best match

def articulate(node, other, particles, xmrcas):     # B-node to A
  matches = filtered_matches(node, other, particles, xmrcas)

  # This is awkward
  if len(matches) == 0:
    return None
  elif len(matches) == 1:
    match = matches[0]
    rematches = filtered_matches(match.cod, cl.get_checklist(node), particles, xmrcas)
    if len(rematches) == 0:
      dribble.log("** Shouldn't happen: %s" % cl.get_unique(node))
      return None
    elif len(rematches) == 1:
      if rel.is_variant(match.relation, rel.eq):
        return art.set_relation(match, rel.eq)
      else:
        return match
    else:
      dribble.log("** Rehelp: %s -> %s" %
            (cl.get_unique(match.cod),
             (", ".join(map(lambda ar:cl.get_unique(ar.cod), rematches)))))
      return match
  else:
    dribble.log("** Help: %s -> %s" %
          (cl.get_unique(node),
             (", ".join(map(lambda ar:cl.get_unique(ar.cod), matches)))))
    return matches[0]

def filtered_matches(node, other, particles, xmrcas):
  assert node > 0
  assert cl.is_accepted(node)
  assert other
  assert cl.get_checklist(node) != other

  exties = extensional_matches(node, other, particles, xmrcas)    #monotypic chain (lineage)

  if len(exties) <= 1:
    return exties
  else:
    namies = intension.intensional_matches(node, other)
    # Nodes we match to intensionally
    targets = [namie.cod for namie in namies]
    bothies = [match for match in exties if match.cod in targets]
    if len(bothies) > 0:
      return [choose_best_match(bothies)]
    else:
      # print("# no extensional + intensional match for %s" % cl.get_unique(node))
      return exties

# ---------- EXTENSIONALITY / by split

# Starting with one match, extend to a set of matches by adding
# 'monotypic' ancestors

def extensional_matches(node, other, particles, xmrcas):       # B-node to A
  assert cl.get_checklist(node) != other
  match = extensional_match(node, other, particles, xmrcas)    # Single B/A extensional match
  if not match:
    return []
  matches = [match]
  if rel.is_variant(match.relation, rel.eq):

    anchor = xmrcas.get(match.cod)

    # Scan upwards from match.code looking for nodes that return back to anchor
    scan = match.cod    # node -> partner
    while True:
      scan = cl.get_parent(scan)    # in other
      if scan == cl.forest_tnu: break
      back = xmrcas.get(scan)
      if back != anchor:
        break
      matches.append(art.monotypic(node, scan, rel.eq))
  return matches

# Guaranteed invertible, except for monotypic node chains
# This code is derived from 
#   reference-taxonomy/org/opentreeoflife/conflict/ConflictAnalysis.java

def extensional_match(node, other, particles, xmrcas):
  part = particles.get(node)
  if part:
    return part
  partner = xmrcas.get(node)      # node in other checklist; 'conode'
  if not partner:
    # Descendant of a particle
    return None
  back = xmrcas.get(partner)    # 'bounce'
  if not back:
    # Not sure how this can happen but it does (NCBI vs. GBIF)
    dribble.log("%s <= %s <= nowhere" % (cl.get_unique(node),
                                         cl.get_unique(partner)))
    return None
  how = cl.how_related(back, node)    # Alway rcc5
  # assert how != rel.disjoint - oddly, not always true
  reason = "cross-mrca " + how.name
  if how == rel.eq:
    # Should end up being eq iff name match or unique match
    # Can test for unique match by looking at xmrca of parent

    # Could be part of a 'monotypic' chain; fix later
    re = rel.extensional
  elif how != rel.gt:
    re = how
  else:               # must be rel.lt
    # Assume resolution (node < partner) until conflict is proven
    re = rel.lt

    # Look for an intersection between any partner-child and node
    # x is in A checklist, y is in B checklist
    for pchild in cl.get_children(partner):
      pchild_back = xmrcas.get(pchild)
      if pchild_back == None:
        # pchild ! node
        pass
      else:
        (x, y) = cross_compare(node, pchild, xmrcas)
        if x and y:
          re = rel.conflict
          reason = ("%s is in the A taxon; B-sibling %s is not" %
                    (cl.get_unique(x), cl.get_unique(y)))
          break
        elif y:
          reason = ("%s is not in the A taxon" % cl.get_unique(y))

  return art.extensional(node, partner, re, reason)

# ---------- Determine disjointness across checklists

# Compare node (in A) to conode (in B) according to particle sets.
# They can bear any relationship: < = > >< !

# Returns subnodes (x, y) of conode where
#    x <= node ∩ conode, or is None if node ∩ conode = ∅
#    y <= conode - node, or is None if conode - node = ∅

# We're looking for two things
#  1. whether they intersect (if so then parent conflict can be detected) -
#      i.e. a conode descendant that *is* a node descendant.
#  2. why they are not equal (goes to conflict diagnosis) - example of
#      conode descendant that is *not* a node descendant.

# On each recursive call, we can provide one or the other of these, or both.
# Or maybe none.

def cross_compare(node, conode, xmrcas):
  back = xmrcas.get(conode)
  if back == None:
    return (None, None)
  how = cl.how_related(back, node)
  assert how != rel.conflict
  if how == rel.eq:
    return (conode, None)
  elif how == rel.lt:
    return (conode, None)
  elif how == rel.disjoint:
    return (None, conode)
  assert how == rel.gt
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

def analyze_cross_mrcas(A, B, particles):
  cross_mrcas = {}
  def half_analyze_cross_mrcas(checklist, other):
    def subanalyze_cross_mrcas(node, other):
      result = None
      probe = particles.get(node)
      if probe:
        assert probe.dom == node
        assert cl.is_accepted(probe.cod)
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
      dribble.log("# ** No return cross-MRCA for %s -> %s -> ..." %\
                  (cl.get_unique(node), cl.get_unique(cross)))
  return cross_mrcas

# ---------- Particles

# A particle is a mutual =-articulation of tipward accepted nodes
# deriving from sameness of 'intrinsic' node properties: name, id,
# rank, parent, etc.

# This function returns a partial map from nodes to articulations.

def check_mutuality(particles):
  for node in particles:
    ar = particles.get(node)
    if ar:
      rar = particles.get(ar.cod)
      if rar:
        if rar.cod == node:
          pass
        else:
          dribble.log(
              "** Round trip fail:\n  %s\n  %s\n" %
              (art.express(ar), art.express(rar)))
      else:
        dribble.log("** No return match: %s -> none" %
                    (art.express(ar)))

# Filter out nodes where the relationship is not mutual

def mutual(am):
  mut = {}
  for node in am:
    ar = am[node]               # articulation
    back = am.get(ar.cod)
    if back and back.cod == node:
      mut[node] = ar
  return mut

# Filter out internal nodes (those having a mapped descendant)

def tipward(amap, A, B):
  tw = {}
  def filter_checklist(c):
    for root in cl.get_roots(c):
      filter(root)
  def filter(node):
    found_match = False
    for child in cl.get_children(node):
      if filter(child):
        found_match = True
    if found_match:    # Some descendant is a particle
      return True
    elif node in amap:
      tw[node] = amap[node]
      return True
    else:
      return False
  filter_checklist(A)
  filter_checklist(B)
  return tw

