debug = False

# Temporary hack for experimenting with poorly formed EOL checklists
EOL = False

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
  particles = find_particles(B, A)
  particles = fix_alignment(particles)

  cross_mrcas = analyze_cross_mrcas(B, A, particles)
  print ("# Number of cross-mrcas:", len(cross_mrcas))

  the_alignment = finish_alignment(B, A, particles, cross_mrcas)
  print ("# Number of articulations in alignment:", len(the_alignment))

  return (the_alignment, cross_mrcas)

def finish_alignment(B, A, particles, xmrcas):
  global cross_mrcas
  cross_mrcas = xmrcas
  draft = {}
  def half_finish(B, A):
    def process(node):
      m = articulate(node, A, particles, xmrcas)
      if m:
        if not EOL:
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

# Deal with splits and merges.

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

  for node in incoming:          # y node
    # e.g. inc = {a ≈ node, b ≈ node}
    inc = incoming[node]    # articulations whose cod is node
    # mut: node ≈ a
    mut = draft.get(node)   # mut : node -> mutual
    sibs = []
    if mut:
      rent = cl.get_parent(mut.cod)    # parent of a
      for ar in inc:    # ar: a ≈ node
        if (rel.is_variant(ar.relation, rel.eq) and
            # parent of b is same as parent of a?
            cl.get_parent(ar.dom) == rent):
          sibs.append(ar)
    for ar in inc:
      target = ar.dom
      if ar in sibs:
        if len(sibs) == 1:
          alignment[target] = art.set_relation(ar, rel.eq)
        else:
          alignment[target] = art.change_relation(ar, rel.merge, "merge", "split")
    if len(sibs) > 1:
      # alignment[node] = art.change_relation(ar, rel.split, "split", "merge")
      # Report!
      print("# Split/join %s -> %s -> %s" %
            (" ∨ ".join(map(lambda e:cl.get_unique(e.dom), sibs)),
             cl.get_unique(node),
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
      print("** Shouldn't happen: %s" % cl.get_unique(node))
      return None
    elif len(rematches) == 1:
      if rel.is_variant(match.relation, rel.eq):
        return art.set_relation(match, rel.eq)
      else:
        return match
    else:
      print("** Rehelp: %s -> %s" %
            (cl.get_unique(match.cod),
             (", ".join(map(lambda ar:cl.get_unique(ar.cod), rematches)))))
      return match
  else:
    print("** Help: %s -> %s" %
          (cl.get_unique(node),
             (", ".join(map(lambda ar:cl.get_unique(ar.cod), matches)))))
    return matches[0]

def filtered_matches(node, other, particles, xmrcas):
  assert node > 0
  assert cl.is_accepted(node)
  assert other
  assert cl.get_checklist(node) != other

  exties = extensional_matches(node, particles, other)    #monotypic chain (lineage)

  if len(exties) <= 1:
    return exties
  else:
    namies = intension.intensional_matches(node, other)
    # Nodes we match to intensionally
    targets = [namie.cod for namie in namies]
    bothies = [match for match in exties if match.cod in targets]
    if len(bothies) > 0:
      return intension.skim_best_matches(bothies)
    else:
      # print("# no extensional + intensional match for %s" % cl.get_unique(node))
      return exties

# ---------- EXTENSIONALITY / by split

# Starting with one match, extend to a set of matches by adding
# 'monotypic' ancestors

def extensional_matches(node, particles, other):       # B-node to A
  assert cl.get_checklist(node) != other
  match = extensional_match(node, particles, other)    # Single B/A extensional match
  if not match:
    return []
  matches = [match]
  if rel.is_variant(match.relation, rel.eq):

    here = cl.get_checklist(node)
    anchor = cross_mrca(match.cod)

    # Scan upwards from match.code looking for nodes that return back to anchor
    scan = match.cod    # node -> partner
    while True:
      scan = cl.get_parent(scan)    # in other
      if scan == cl.forest_tnu: break
      back = cross_mrca(scan)
      if back != anchor:
        break
      matches.append(art.monotypic(node, scan, rel.eq))
  return matches

# Guaranteed invertible, except for monotypic node chains
# This code is derived from 
#   reference-taxonomy/org/opentreeoflife/conflict/ConflictAnalysis.java

def extensional_match(node, particles, other):
  part = particles.get(node)
  if part:
    return part
  partner = cross_mrca(node)      # node in other checklist; 'conode'
  if not partner:
    return None
  back = cross_mrca(partner)    # 'bounce'
  if not back:
    return None    # Not sure how this can happen but it does (NCBI vs. GBIF)
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
      pchild_back = cross_mrca(pchild)
      if pchild_back == None:
        # pchild ! node
        pass
      else:
        (x, y) = cross_compare(node, pchild)
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

def cross_compare(node, conode):
  back = cross_mrca(conode)
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
    saw = cross_compare(node, child)
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
  def half_analyze_cross_mrcas(checklist, other, checkp):
    def subanalyze_cross_mrcas(tnu, other):
      result = None
      probe = particles.get(tnu)
      if probe:
        assert probe.dom == tnu
        if not EOL:
          assert cl.is_accepted(probe.cod)
        if debug:
          print("#   particle(%s) = %s" %\
                (cl.get_unique(tnu), cl.get_unique(probe.cod)))
        result = probe.cod
      else:
        children = cl.get_children(tnu)
        if children:
          m = None      # None is the identity for mrca
          for child in children:
            if debug:
              print("# Child %s of %s" % (cl.get_unique(child), cl.get_unique(tnu)))
            m2 = subanalyze_cross_mrcas(child, other)
            if m2 != None:
              if debug:
                print("#  Folding %s into %s" % (cl.get_unique(m2), cl.get_unique(m)))
              m = cl.mrca(m, m2) if m != None else m2
              # ?????
              if debug:
                print("#   -> %s" % cl.get_unique(m))
          if m != None:
            result = m
            if debug:
              print("#   cm(%s) = %s)" % \
                    (cl.get_unique(tnu), cl.get_unique(m)))
      cross_mrcas[tnu] = result
      if debug and checkp and result:
        probe = cross_mrcas.get(result)
        if not probe:
          print("# ** No return cross-MRCA for %s -> %s -> ..." %\
                (cl.get_unique(tnu), cl.get_unique(result)))
      return result
    for root in cl.get_roots(checklist):
      subanalyze_cross_mrcas(root, other)
  half_analyze_cross_mrcas(A, B, False)
  half_analyze_cross_mrcas(B, A, True)
  return cross_mrcas

# Returns an accepted/accepted articulation

def cross_mrca(node):
  global cross_mrcas
  assert node > 0
  if not EOL:
    assert cl.is_accepted(node)
  return cross_mrcas.get(node)

# ---------- Particles

# A particle is mutual =-articulation of accepted nodes deriving from
# sameness of 'intrinsic' node properties: name, id, rank, parent,
# etc.

# This function returns a partial map from nodes to articulations.

def find_particles(here, other):
  intension.clear_cache()
  particles = {}
  count = [0]
  best = intension.best_intensional_matches(here, other)
  def log(node, message):
    if count[0] < 0:
      if debug:
       print("# fp(%s): %s" % (cl.get_unique(node), message))
      count[0] += 1
  def subanalyze(node, other):
    log(node, "subanalyze")
    if not cl.is_accepted(node):
      print("# ** Child %s of %s has an accepted name" %
            (cl.get_unique(node), cl.get_unique(cl.get_parent(node))))
      return False
    found_match = False
    for inf in cl.get_children(node):
      log(node, "child %s" % inf)
      if subanalyze(inf, other):
        found_match = True
    if found_match:    # Some descendant is a particle
      return True
    candidate = best.get(node)
    if candidate:
      rematch = best.get(candidate.cod)
      if rematch:
        if rematch.cod == node:
          if not EOL and not cl.is_accepted(candidate.cod):
            print("# ** Candidate is synonymlike: %s" % cl.get_unique(candidate.cod))
          particles[node] = candidate    # here -> other
          particles[candidate.cod] = art.reverse(candidate)  # other -> here
          return True
        else:
          # This situation probably reflects a split!
          log(node,
              "Round trip fail:\n  %s\n  %s\n" %
              (art.express(candidate),
               art.express(art.reverse(rematch))))
      else:
        log(node, "No rematch")
    return False
  log(0, "top")
  for root in cl.get_roots(here):
    log(root, "root")
    subanalyze(root, other)
  print ("# Number of particles:", len(particles)>>1)
  return particles

