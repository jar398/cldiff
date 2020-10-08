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
  particles = intension.find_particles(B, A)

  cross_mrcas = analyze_cross_mrcas(B, A, particles)
  print ("# Number of cross-mrcas:", len(cross_mrcas))

  the_alignment = finish_alignment(B, A, particles, cross_mrcas)
  print ("# Number of articulations in alignment:", len(the_alignment))

  return (the_alignment, cross_mrcas)

def finish_alignment(B, A, particles, xmrcas):
  global cross_mrcas
  cross_mrcas = xmrcas
  alignment = {}
  def half_finish(B, A):
    def process(node):
      m = articulate(node, A, particles, xmrcas)
      if m:
        if not EOL:
          assert cl.is_accepted(m.cod)
        alignment[node] = m
      else:
        # Otherwise, what?
        pass
      for child in cl.get_children(node):
        process(child)
    for root in cl.get_roots(B):
      process(root)
  half_finish(B, A)
  half_finish(A, B)
  fix_alignment(alignment)
  return alignment

# This is a kludge.

def fix_alignment(alignment):
  errants = []
  for node in alignment:
    m1 = alignment[node]
    if m1.relation == rel.eq:
      m2 = alignment.get(m1.cod)
      if (m2 == None or \
          m2.relation != rel.eq or \
          m2.cod != node):
        errants.append((m1, m2))
  for (m1, m2) in errants:
    node = m1.dom
    if m2 == None:
      ar = art.change_relation(m1, rel.lt, "foo?")
      alignment[node] = ar
      alignment[m1.cod] = art.reverse(ar)
    elif m2.relation != rel.eq:
      alignment[node] = art.change_relation(m1, rel.reverse(m2.relation), "reverse")
    else:    # m2.cod != node
      ca = cl.mrca(m2.cod, node)
      if ca == node:
        alignment[node] = art.change_relation(m1, rel.lt, "chain")
      elif ca == m2.cod:
        alignment[node] = art.change_relation(m1, rel.gt, "chain")
      else:
        # Split?  What to do?
        # e.g. ** B.Callosciurus_caniceps ∨ A.Petaurista_caniceps = B.Sciuridae
        print("# ** %s ∨ %s <= %s" %
              (cl.get_unique(node), cl.get_unique(m1.cod), cl.get_unique(ca)),
              file=dribble.dribble_file)
  return None

def is_mutual(m1, al):
  m2 = al[m1.cod]
  return m2 and art.inverses(m1, m2)

# ---------- One-sided best match

def articulate(node, other, particles, xmrcas):     # B-node to A
  assert node > 0
  assert cl.is_accepted(node)
  assert other
  assert cl.get_checklist(node) != other

  if debug:
    print("# matching B-node", cl.get_unique(node))
  exties = extensional_matches(node, particles, other)    #monotypic chain

  if len(exties) == 1:
    match = exties[0]
  else:
    # Filter out intensionals that aren't extensional (or v.v.)
    namies = intension.intensional_matches(node, other)
    targets = [namie.cod for namie in namies]
    #         [extie.cod for extie in exties]
    #bothies = [match for match in namies if match.cod in targets]
    bothies = [match for match in exties if match.cod in targets]

    if len(bothies) > 0:
      match = choose_best_match(bothies)
    elif len(exties) > 0:
      if debug: print("# %s exties to %s" % (len(exties), cl.get_unique(node)))
      match = choose_best_match(exties)
    else:
      # should not happen
      if debug: print("# no extie")
      match = choose_best_match(namies)

  return match

# ---------- EXTENSIONALITY / by split

# Starting with one match, extend to a set of matches by adding
# 'monotypic' ancestors

def extensional_matches(tnu, particles, other):       # B-node to A
  assert cl.get_checklist(tnu) != other

  match = extensional_match(tnu, particles, other)    # Single B/A extensional match
  if not match:
    return []
  matches = [match]
  if debug: print("# %s has initial extensional match %s" %
                  (cl.get_unique(tnu), cl.get_unique(match.cod)))
  if rel.is_variant(match.relation, rel.eq):

    here = cl.get_checklist(tnu)
    anchor = cross_mrca(match.cod)

    # Scan upwards from match.code looking for nodes that return back to anchor
    scan = match.cod    # tnu -> partner
    if debug: print("# scanning up from %s for returns to %s" %
                    (cl.get_unique(scan), cl.get_unique(anchor)))

    while True:
      scan = cl.get_parent(scan)    # in other
      if scan == None:
        print("# bailing out, no parent %s" % cl.get_unique(tnu))
        break
      if scan == cl.forest_tnu:
        print("# bailing out, %s --> forest" % cl.get_unique(tnu))
        break
      if debug:
        print("# considering %s" % cl.get_unique(scan))
      back = cross_mrca(scan)
      if back != anchor:
        if debug: print("# goes back to %s not %s, breaking out" %
                        (cl.get_unique(back), cl.get_unique(anchor)))
        break
      if debug: print("# adding return match %s" %
                      (cl.get_unique(scan)))
      matches.append(art.monotypic(tnu, scan, rel.eq))
    if debug: print("# %s matches %s nodes extensionally" %
                    (cl.get_unique(tnu), len(matches)))
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
  how = cl.how_related(node, back)
  # assert how != rel.disjoint - oddly, not always true
  reason = "particle set " + how.name
  if how == rel.eq:
    # Could be part of a 'monotypic' chain; fix later
    re = how
  elif how != rel.lt:
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
  re = cl.how_related(node, back)
  assert re != rel.conflict
  if re == rel.eq:
    return (conode, None)
  elif re == rel.gt:
    return (conode, None)
  elif re == rel.disjoint:
    return (None, conode)
  assert re == rel.lt
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

def cross_mrca(tnu):
  global cross_mrcas
  assert tnu > 0
  if not EOL:
    assert cl.is_accepted(tnu)
  return cross_mrcas.get(tnu)
