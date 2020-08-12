debug = False

# Temporary hack for experimenting with poorly formed EOL checklists
EOL = False

import sys, csv

import checklist as cl
import relation as rel
import articulation as art

# For each B-record, we choose an articulation with the closest
# A-record that it matches (preferably but not necessarily an '='
# articulation).

def align(B, A, dribfile=None):
  global intensional_matches_cache
  global dribble_file
  dribble_file = dribfile if dribfile else sys.stdout
  intensional_matches_cache = {}
  particles = find_particles(B, A)
  print ("# Number of particles:", len(particles)>>1)

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
      for child in get_children(node):
        process(child)
    for root in cl.get_roots(B):
      process(root)
  half_finish(B, A)
  half_finish(A, B)
  return alignment

def mutual_match(tnu, alignment):
  m1 = alignment[tnu]
  if m1:
    m2 = alignment[m1.cod]
    if m2 and m2.cod == tnu:
      return m1
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
    namies = intensional_matches(node, other)
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
        print("# bailing out, no parent %s" % scan)
        break
      if scan == cl.forest_tnu:
        print("# bailing out, forest")
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
  if how == rel.eq:
    # Could be part of a 'monotypic' chain; fix later
    re = how
    reason = "same particle set"
  elif how != rel.lt:
    re = how
    reason = "particle set comparison"
  else:               # must be rel.lt
    # Assume resolution (<) until conflict is proven
    re = rel.lt
    reason = "particle set containment"

    # Look for an intersection between any partner-child and node
    x_seen = None
    y_seen = None
    for pchild in get_children(partner):
      pchild_back = cross_mrca(pchild)
      if pchild_back == None:
        graft = pchild
      else:
        (x, y) = cross_compare(node, pchild)
        if x and y:
          re = rel.conflict
          reason = ("%s is in; %s is outside" %
                    (cl.get_unique(x), cl.get_unique(y)))
          break
        elif y:
          y_seen = y
        elif x:
          x_seen = x
    if not reason:
      reason = ("%s is in; %s is outside" % (cl.get_unique(x), cl.get_unique(y)))

  return art.extensional(node, partner, re, reason)

# ---------- Determine disjointness across checklists

# node bears any relationship to conode.
# We're looking for two things
#  1. whether they intersect (if so then parent conflict can be detected) -
#      i.e. a conode descendant that *is* a node descendant.
#      There may not be one.
#  2. why they are not equal (goes to conflict diagnosis) - example of
#      conode descendant that is *not* a node descendant.
#      There is always one of these.

# On each recursive call, we can provide one or the other of these, or both.
# Or maybe none.

# Returns (x, y) where x (under conode) is contained in node,
# and y (also under conode) is disjoint with node.  Either can be None.

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
  for child in get_children(conode):
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
        children = get_children(tnu)
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

# ---------- Particles

# A particle is mutual =-articulation of accepted nodes deriving from
# sameness of 'intrinsic' node properties: name, id, rank, parent,
# etc.

def find_particles(here, other):
  particles = {}
  count = [0]
  def log(tnu, message):
    if count[0] < 0:
      if debug:
       print("# fp(%s): %s" % (cl.get_unique(tnu), message))
      count[0] += 1
  def subanalyze(tnu, other):
    log(tnu, "subanalyze")
    if not cl.is_accepted(tnu):
      print("# ** Child %s of %s has an accepted name" %
            (cl.get_unique(tnu), cl.get_unique(cl.get_parent(tnu))))
      return False
    found_match = False
    for inf in get_children(tnu):
      log(tnu, "child %s" % inf)
      if subanalyze(inf, other):
        found_match = True
    if found_match:    # Some descendant is a particle
      return True
    candidate = best_intensional_match(tnu, other)
    if candidate:
      rematch = best_intensional_match(candidate.cod, here)
      if rematch:
        if rematch.cod == tnu:
          if not EOL and not cl.is_accepted(candidate.cod):
            print("# ** Candidate is synonymlike: %s" % cl.get_unique(candidate.cod))
          particles[tnu] = candidate    # here -> other
          particles[candidate.cod] = art.reverse(candidate)  # other -> here
          return True
        else:
          # This situation probably reflects a split!
          log(tnu,
              "Round trip fail:\n  %s\n  %s\n" %
              (art.express(candidate),
               art.express(art.reverse(rematch))))
      else:
        log(tnu, "No rematch")
    return False
  log(0, "top")
  for root in cl.get_roots(here):
    log(root, "root")
    subanalyze(root, other)
  return particles

# ---------- 'Intensional' matches to accepted nodes

# The source node ('tnu') may be accepted or a synonym.

# Three components: synonym-or-self o intensional o synonym-of-or-self
# from_accepted_articulations o intensional_matches o [to_accepted_articulation]

def best_intensional_match(tnu, other):
  matches = intensional_matches(tnu, other)
  return choose_best_match(matches)

# Matches where the target is accepted.

def intensional_matches(tnu, other):
  global intensional_matches_cache
  probe = intensional_matches_cache.get(tnu)
  if probe: return probe
  matches = weak_intensional_matches(tnu, other)
  matches = art.collapse_matches(matches)
  intensional_matches_cache[tnu] = matches
  return matches

# Like strong matches, but adds matches via source synonyms

def weak_intensional_matches(tnu, other):
  matches = strong_intensional_matches(tnu, other)
  syns = synonyms_locally(tnu)
  assert not syns or syns[0].dom == tnu
  matches = matches + [art.compose(syn, match)
                       for syn in syns
                       for match in strong_intensional_matches(syn.cod, other)
                       if art.composable(syn, match)]
  assert is_matches(matches)
  return matches

def strong_intensional_matches(tnu, other):
  matches = direct_matches(tnu, other)
  assert is_matches(matches)
  matches = [match_to_accepted(match) for match in matches]
  assert is_matches(matches)
  return matches

# Convert match a->b to a->b' where b' is accepted (maybe b=b')

def match_to_accepted(m):
  if cl.is_accepted(m.cod):
    return m
  else:
    return art.compose(m, has_accepted_locally(m.cod))

# Intensional matches by name (no synonym following)
# This ought to be cached I think?

def direct_matches(node, other):
  hits = cl.get_similar_records(other, node)
  return [art.intensional(node, hit) for hit in hits]

# ---------- Within-checklist articulations

# Handy for composing paths.

# Synonym-or-self = to-accepted
# Matches that involve going from an accepted name to a synonym are weak

# An accepted name can have many synonyms

def synonyms_locally(tnu):
  if cl.is_accepted(tnu):
    hases = [has_accepted_locally(syn) for syn in cl.get_synonyms(tnu)]
    return art.collapse_matches([art.reverse(ar) for ar in hases if ar])
  else:
    return []

# A synonym has only one accepted name

def has_accepted_locally(maybe_syn):   # goes from synonym to accepted
  assert maybe_syn > 0
  accepted = cl.get_accepted(maybe_syn)
  if accepted:
    if not EOL:
      assert accepted != maybe_syn
    return art.synonymy(maybe_syn, accepted)
  else:
    return None

# ---------- Best match among a set with common domain

def choose_best_match(arts):     # => art
  assert is_matches(arts)
  if len(arts) == 0: return None
  arts = choose_best_matches(arts)
  b = arts[0]
  if len(arts) == 1: return b
  print("# ** Multiple least-bad matches. Need to find more tie-breakers.",
        file=dribble_file)
  print("   %s -> %s" %
        (cl.get_unique(b.dom),
         [cl.get_unique(a.cod) for a in arts]),
        file=dribble_file)
  return None

# There can be multiple best matches

def choose_best_matches(arts):
  if len(arts) == 0:
    return []
  else:
    sorted_matches = art.sort_matches(art.collapse_matches(arts))
    badness = art.badness(sorted_matches[0])
    matches = []
    for match in sorted_matches:
      if art.badness(match) == badness:
        matches.append(match)
      else:
        break
    return matches

# Random

def get_children(node):
  assert node > 0
  return [node
          for node in cl.get_children(node) + cl.get_synonyms(node)
          if cl.is_accepted(node)]

# ---------- For debugging

def is_match(ar):
  c1 = cl.get_checklist(ar.dom)
  c2 = cl.get_checklist(ar.cod)
  if c1 == c2:
    print("%s and %s both in %s" %
          (cl.get_name(c1),
           c2.get_name(c2),
           c1.prefix))
    return False
  return True

def is_matches(matches):
  if len(matches) == 0: return True
  return is_match(matches[0])

def cache_it(fun, tnu, other, dict):
  not_cached = "not cached"
  probe = dict.get(tnu, not_cached)
  if probe != not_cached:
    return probe
  else:
    result = fun(tnu, other)
    dict[x] = result
    return result

# ----

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('left', help='A checklist')
  parser.add_argument('right', help='B checklist')
  parser.add_argument('--left-tag', default="A")
  parser.add_argument('--right-tag', default="B")
  parser.add_argument('--idspace', default=False)
  parser.add_argument('--out', help='file name for report', default='diff.csv')
  parser.add_argument('--format', help='report format', default='ad-hoc')
  args = parser.parse_args()
  main(args.left, args.left_tag, args.right, args.right_tag, args.out, args.format)

