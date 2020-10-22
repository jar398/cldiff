# This is a godawful mess

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

  # Start to fill in the alignment
  tipward_best = tipward(best, B, A)
  draft = detect_lumps_splits(tipward_best)    # change ≈ to =

  # Extensional analysis
  # particles = mutual(tipward_best)
  particles = tipward(mutual(best), B, A)
  check_mutuality(particles)
  dribble.log("# Number of particles: %s" % (len(particles)>>1))
  cross_mrcas = analyze_cross_mrcas(B, A, particles)
  dribble.log("# Number of cross-mrcas: %s" % len(cross_mrcas))
  extensionals = extensional_match_map(A, B, particles, cross_mrcas)
  dribble.log("# Number of extensional relationships: %s" % len(extensionals))

  # Add extensional matches to a draft that already has intensional matches
  the_alignment = assemble_alignment(draft, best, extensionals)

  dribble.log("# Number of articulations in alignment: %s" % len(the_alignment))

  return (the_alignment, cross_mrcas)

def assemble_alignment(draft, best, ext_map):
  for node in ext_map:
    ar1 = draft.get(node)
    ar2 = articulate(node, best, ext_map)
    if ar1 and ar2 and ar1.cod != ar2.cod:
      dribble.log("** Intensional and extensional matches differ:\n  %s\n  %s" %
                  (art.express(ar1), art.express(ar2)))
    if ar2:
      if (dribble.watch(node)):
        dribble.log("# Matching extensionally: %s" %
                    art.express(ar2))
      draft[node] = ar2
    elif ar1:
      if (dribble.watch(node)):
        dribble.log("# Matching intensionally: %s" %
                    art.express(ar1))
      draft[node] = ar1
    else:
      if (dribble.watch(node)):
        dribble.log("# No match for %s" % cl.get_unique(node))
  return draft

# ---------- Final alignment

# Determine the appropriate relationship to put in the final
# alignment.  Note that best has all intensional best matches, not
# just the particles.

def articulate(node, best, ext_map):     # B-node to A
  matches = filtered_matches(node, best, ext_map)

  # This is awkward
  if len(matches) == 0:
    ar = ext_map.get(node)
    if dribble.watch(node):
      dribble.log("# Non-equal articulation: %s" % art.express(ar))
    return ar
  elif len(matches) == 1:
    match = matches[0]
    rematches = filtered_matches(match.cod, best, ext_map)
    if len(rematches) == 0:
      dribble.log("** Shouldn't happen: %s" % cl.get_unique(node))
      return None
    elif len(rematches) == 1:
      if rel.is_variant(match.relation, rel.eq):
        return art.set_relation(match, rel.eq)
      else:
        return match
    else:
      dribble.log("** Multiple returns: %s -> %s" %
            (cl.get_unique(match.cod),
             (", ".join(map(lambda ar:cl.get_unique(ar.cod), rematches)))))
      return match
  else:
    dribble.log("** Multiple matches: %s -> %s" %
          (cl.get_unique(node),
             (", ".join(map(lambda ar:cl.get_unique(ar.cod), matches)))))
    return matches[0]

# Find matches that are both extensional and intensional.
# Returns a list of =-like articulations, perhaps empty.

def filtered_matches(node, best, ext_map):
  assert node > 0
  assert cl.is_accepted(node)

  exties = extensional_chain(node, ext_map)    #monotypic chain (lineage)
  if dribble.watch(node):
    dribble.log("# Chain %s = [%s]" %
                (cl.get_unique(node),
                 " <- ".join(map(lambda ar:cl.get_unique(ar.cod), exties))))

  if len(exties) < 1:
    if dribble.watch(node):
      dribble.log("# No extensionally equal matches to %s" % cl.get_unique(node))
    return exties
  else:
    bestie = best.get(node)
    if bestie and bestie.cod in [e.cod for e in exties]:
      if dribble.watch(node):
        dribble.log("# Extensional + intensional match: %s" % art.express(bestie))
      return [art.change_relation(bestie, rel.eq, "particleset=", "particleset=")]
    else:
      if dribble.watch(node):
        dribble.log("# No applicable intensional match to %s" % cl.get_unique(node))
      return exties

# ---------- Find preimage of extensional map (a chain)

# Starting with one match, extend to a set of extensional matches by
# adding 'monotypic' ancestors.
# Returns a list of =-like articulations, perhaps empty.

def extensional_chain(x, ext_map):       # B-node to A
  match = ext_map.get(x)    # Single B/A extensional match
  if dribble.watch(x):
    dribble.log("# Chain(%s) - articulation is %s" %
                (cl.get_unique(x), art.express(match)))
  if not match:
    return []
  if not rel.is_variant(match.relation, rel.eq):
    return []

  # match.cod is the start of the chain
  # anchor is a descendant of x
  y = match.cod

  # Start at anchor's match and go rootward finding all matches to anchor
  matches = []
  x_anchor = None
  while True:
    back = ext_map.get(y)
    if (back == None or
        not rel.is_variant(back.relation, rel.eq)):
      if dribble.watch(x):
        dribble.log("# Chain %s - no back-articulation from %s: %s" %
                    (cl.get_unique(x), cl.get_unique(y), art.express(back)))
      break
    if x_anchor == None:
      x_anchor = back.cod
      if dribble.watch(x):
        dribble.log("# Chain - setting anchor: %s" % cl.get_unique(x_anchor))
    elif back.cod != x_anchor:
      if dribble.watch(x):
        dribble.log("# Chain - end: %s" % cl.get_unique(back.cod))
      break
    if x_anchor == x:
      if y == match.cod:
        m = match               # x -> y
      else:
        m = art.reverse(back)   # x -> y
    else:
      m = art.monotypic(x, y, rel.extensional)
    matches.append(m)
    y = cl.get_parent(y)
  #if len(matches) > 1:
  #  print ("# nontrivial chain for %s" % cl.get_unique(x))
  return matches

# ---------- Extensionality by particle set

# Extensional relationships only, no particle matches

def extensional_match_map(A, B, particles, xmrcas):
  ext = {}
  def process(node):
    for child in cl.get_children(node):
      process(child)
    ar = extensional_match(node, particles, xmrcas)
    if ar:
      ext[node] = ar
  for root in cl.get_roots(A): process(root)
  for root in cl.get_roots(B): process(root)
  return ext

# Guaranteed invertible, except for monotypic node chains
# This code is derived from 
#   reference-taxonomy/org/opentreeoflife/conflict/ConflictAnalysis.java

def extensional_match(node, particles, xmrcas):
  part = particles.get(node)
  if part:
    if dribble.watch(node):
      dribble.log("# EM: particle. %s" % art.express(part))
    return part
  partner = xmrcas.get(node)      # node in other checklist; 'conode'
  if not partner:
    # Descendant of a particle
    if dribble.watch(node):
      dribble.log("# EM: %s is descendant of particle." % cl.get_unique(node))
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
    how = rel.extensional
    reason = "mutual-cross-mrca"
  elif how == rel.gt:
    how = rel.extensional
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
        if d and e:
          how = rel.conflict
          reason = ("%s is in; its sibling %s is not" %
                    (cl.get_unique(d), cl.get_unique(e)))
          dribble.log("** %s (x) conflicts with %s (y) because:\n"
                      "   %s (c) is a child of y,\n"
                      "   %s (d) and %s (e) are descendants of c,\n"
                      "   and d < x while e ! x" %
                      (cl.get_unique(node),
                       cl.get_unique(partner),
                       cl.get_unique(pchild),
                       cl.get_unique(d),
                       cl.get_unique(e)))
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


# ---- Deal with ad hoc splits and merges among tipward best matches.

def detect_lumps_splits(best):
  result = {}

  # Suppose `node` x comes from the A checklist, and there is a split
  # such that x matches multiple nodes y1, y2 in the B checklist.

  # For each A-node x, find all B-nodes y, ... that are contending to
  # be equivalent to x.

  incoming = {}
  for node in best:    # y node
    ar = best[node]     # ar : y -> x
    if ar.cod in incoming:
      incoming[ar.cod].append(ar) # x -> (y->x, y2->x)
    else:
      incoming[ar.cod] = [ar]     # x -> (y->x)

  # Modify the relation for all approximate-match nodes.
  flush = []
  for y in incoming:          # y node
    # e.g. inc = {a ≈ y, b ≈ y}
    inc = incoming[y]    # cod = y for all articulations
    # mut: y ≈ x
    mut = best.get(y)   # mut : y -> mutual
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
          result[x] = art.set_relation(ar, rel.eq) # x = y
        else:
          result[x] = art.change_relation(ar, rel.merge, "merge", "split")
      else:
        result[x] = ar       # x < y
    if len(sibs) > 1:
      flush.append(y)
      # Report!
      dribble.log("# Split/lump %s -> %s -> %s" %
                  (" ∨ ".join(map(lambda e:cl.get_unique(e.dom), sibs)),
                   cl.get_unique(y),
                   cl.get_unique(rent)))

  for y in flush:
    del result[y]
      
  return result

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
    debug = dribble.watch(node)
    ar = am[node]               # articulation
    back = am.get(ar.cod)
    if back:
      if back.cod == node:
        mut[node] = ar
        if debug: dribble.log("# Mutually matched: %s" % 
                              (art.express(ar)))
      else:
        if debug: dribble.log("# Not mutual: %s, %s" %
                              (art.express(ar), art.express(back)))
    else:
      if debug: dribble.log("# No return: %s" % art.express(ar))
  return mut

# Filter out internal nodes (those having a mapped descendant)

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
      tw[node] = ar
      if debug: dribble.log("# %s is a tipward match, keeping: %s" %
                            (cl.get_unique(node), art.express(ar)))
      return ar
    else:
      if debug: dribble.log("# %s is unmatched" % cl.get_unique(node))
      return None
  for root in cl.get_roots(A): filter(root)
  for root in cl.get_roots(B): filter(root)
  return tw

