debug = False

import sys, csv
import argparse

import checklist as cl
import relation as rel
import articulation as art
import diff

# ---------- OVERALL

def start(A, B):
  global particles
  global cross_mrcas
  global inverse_good_candidates
  global grafts

  particles = find_particles(A, B)
  print ("# Number of particles:", len(particles)>>1)

  cross_mrcas = analyze_cross_mrcas(A, B)
  print ("# Number of cross-mrcas:", len(cross_mrcas))

  alignment = finish_alignment(B, A)
  print ("# Number of alignment articulations:", len(alignment))
  inverse_good_candidates = invert_dict_by_cod(alignment)

  grafts = analyze_unmatched(A, B)
  print ("# Number of grafts: %s\n" % len(grafts))

  return alignment

# Partners list for reporting (A->B articulations)

def good_candidate(tnu, other):  # for children sort order
  return choose_best_match(good_candidates(tnu, other))

def good_candidates(tnu, other):
  matches = inverse_good_candidates.get(tnu) or []
  if len(matches) == 0:
    investigate = good_match(tnu, other)
    if investigate: matches = [investigate]
  return [match for match in matches if not cl.get_accepted(match.cod)]

# ---------- Alignments

# Fill the cache

def invert_dict_by_cod(d):
  inv = {}
  for (node, ar) in d.items():
    if ar:
      ar = art.reverse(ar)
      if ar.dom in inv:
        inv[ar.dom].append(ar)
      else:
        inv[ar.dom] = [ar]
  return inv

def finish_alignment(B, A):
  alignment = {}
  def process(node):
    m = good_match(node, A)
    if m:
      alignment[node] = m
    for child in get_children(node):
      process(child)
  for root in cl.get_roots(B):
    process(root)
  return alignment

# ---------- One-sided best match

def good_match(node, other = None):
  assert node > 0
  assert other
  assert cl.get_checklist(node) != other
  matches = good_matches(node, other)
  match = choose_best_match(matches)
  if match:
    if cl.get_accepted(match.cod) and not cl.get_accepted(match.dom):
      print("# ** Match is supposed to be terminal:\n  %s" %
            (art.express(match)))
    # Happens way often in GBIF
    if False and not cl.is_accepted(match.cod):
      print("# ** Match has taxonomic status %s\n  %s" %
            (cl.get_value(match.cod, cl.taxonomic_status),
             art.express(match)))
  return match

good_matches_cache = {}

def good_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in good_matches_cache: return good_matches_cache[node]

  if debug:
    print("# matching", cl.get_unique(node))
  exties = extensional_matches(node, other)

  if len(exties) > 0:
    if debug: print("# >0 exties")
    if rel.is_variant(exties[0].relation, rel.eq):
      if debug: print("# at least one eq")
      matches = exties + matches_to_accepted(node, other)
    else:
      if debug: print("# no eq")
      matches = exties
  else:
    if debug: print("# no extie")
    matches = matches_to_accepted(node, other)
  matches = collapse_matches(exties)

  good_matches_cache[node] = matches
  return matches

# ---------- EXTENSIONALITY

# Starting with one match, extend to a set of matches by adding
# extensionally equivalent ancestors

def extensional_matches(tnu, other):
  assert cl.get_checklist(tnu) != other

  match = extensional_match(tnu, other)    # Single extensional match
  if not match:
    if debug: print("# no ext match")
    return []
  matches = [match]
  if debug: print("# an ext match")

  if rel.is_variant(match.relation, rel.eq):
    # Scan upwards looking for nodes that return to us...
    scan = match.cod    # tnu -> partner
    here = cl.get_checklist(scan)
    while True:
      scan = cl.get_parent(scan)    # in other
      if not scan: break
      rev = extensional_match(scan, here)
      # rev: other -> here
      if not rev: break
      if not rel.is_variant(rev.relation, rel.eq): break
      if rev.cod != tnu: break
      matches.append(rel.reverse(rev))
  matches.reverse()    # hmm. for choose ?
  if debug: print("# %s ext matches to %s" % (len(matches), cl.get_unique(tnu)))
  return matches

# Guaranteed invertible, except for monotypic node chains

def extensional_match(tnu, other):
  assert tnu > 0
  here = cl.get_checklist(tnu)
  cross = cross_mrca(tnu, other)      # TNU in other checklist
  if not cross:
    return None
  match = how_related_extensionally(cross)
  return match_to_accepted(match)

# Returns a match

def how_related_extensionally(cross):
  tnu = cross.dom
  partner = cross.cod
  here = cl.get_checklist(tnu)
  atcha = cross_mrca(partner, here)
  assert atcha
  back = atcha.cod
  diffs = diff.compose(cross.differences, atcha.differences)
  if cl.are_disjoint(tnu, back):
    re = rel.disjoint
  elif tnu == back:
    particle = particle_match(tnu, partner)
    if particle:
      return art.conjoin(particle,
                         art.extensional(tnu, partner, rel.particle, particle.differences))
    elif diffs == 0:
      return art.extensional(tnu, partner, rel.identical, diffs)
    else:
      return art.extensional(tnu, partner, rel.same_extension, diffs)
  elif cl.mrca(tnu, back) == tnu:
    # Monotypic: back < tnu.
    re = rel.monotypic_over
  else:
    re = rel.lt
    for sub in get_children(partner):
      cross_back = cross_mrca(sub, here)
      back = cross_back.cod
      if back != None:
        if cl.mrca(tnu, back) == tnu and cross_disjoint(tnu, partner):
          re = rel.conflict
  re = extensionally(re)
  return art.extensional(tnu, partner, re, diffs)

def extensionally(re):
  return rel.justify(re, rel.extensionally, "particles " + rel.rcc5_name(re))

# ---------- Determine disjointness across checklists

def cross_disjoint(tnu, partner):
  assert tnu > 0
  assert partner > 0
  assert cl.get_checklist(tnu) != cl.get_checklist(partner)

  cross = cross_mrca(partner, cl.get_checklist(tnu))
  back = cross.cod
  if back == None:
    return True
  assert back > 0
  assert cl.get_checklist(back) == cl.get_checklist(tnu)
  if cl.are_disjoint(tnu, back):
    return True
  for child in get_children(partner):
    assert child > 0
    if not cross_disjoint(tnu, child):
      print("# ** Test %s ! %s failed because not ! %s" %\
            (cl.get_unique(tnu), cl.get_unique(partner), cl.get_unique(child)))
      return False
  return True

# ---------- Cross-MRCAs

# Returns a tnu - we never care about the reason for the match here

def cross_mrca(tnu, other):
  assert tnu > 0
  probe = cross_mrcas.get(tnu, 17)
  if probe != 17:
    return probe
  elif cl.is_accepted(tnu):
    if debug:
      print("# ** No cross-MRCA set for %s" % cl.get_unique(tnu))
  return (None, 1)

# Returns match with diffs = dropped taxa (including
# those of particles)

def analyze_cross_mrcas(A, B):
  cross_mrcas = {}
  def half_analyze_cross_mrcas(checklist, other, checkp):
    def subanalyze_cross_mrcas(tnu, other):
      pmatch = particle_match(tnu, other)
      diffs = diff.no_diffs     # Cumulative diffs for all descendants?
      m = None      # is the identity for mrca
      if pmatch:
        assert pmatch.dom == tnu
        if rel.is_variant(pmatch.relation, rel.eq):
          m = pmatch.cod
          diffs = diff.note_dropped_children(diffs)
          if debug:
           print("#   particle(%s) = %s" %\
                 (cl.get_unique(tnu), cl.get_unique(pmatch.cod)))
        else:
          # Can't happen
          print("#** Non-eq articulation from particle_match\n  %s" +
                art.express(pmatch))
          assert False
      else:
        children = get_children(tnu)
        if children:
          for child in children:
            if debug:
              print("# Child %s of %s" % (cl.get_unique(child), cl.get_unique(tnu)))
            cross = subanalyze_cross_mrcas(child, other)
            if cross:
              m2 = cross.cod
              if debug:
                print("#  Folding %s into %s" % (cl.get_unique(m2), cl.get_unique(m)))
              m = cl.mrca(m, m2) if m != None else m2
              diffs = diff.conjoin(diffs, cross.differences)
              if debug:
                print("#   -> %s" % cl.get_unique(m))
            else:
              diffs = diff.note_dropped_children(diffs)
      if debug:
        print("#   cm(%s) = %s, diffs %o)" % \
              (cl.get_unique(tnu), cl.get_unique(m), diffs))
      result = art.cross_mrca(tnu, m, diffs)
      if result:
        cross_mrcas[tnu] = result
      if debug and checkp:
        probe = cross_mrcas.get(m)
        if not probe:
          print("# ** No return cross-MRCA for %s -> %s -> ..." %\
                (cl.get_unique(tnu), cl.get_unique(m)))
      return result
    for root in cl.get_roots(checklist):
      subanalyze_cross_mrcas(root, other)
  half_analyze_cross_mrcas(A, B, False)
  half_analyze_cross_mrcas(B, A, True)
  return cross_mrcas

cross_mrcas = {}

# ---------- Particles

# A particle is mutual articulation deriving from sameness of 'intrinsic'
# node properties: name, id, rank, parent, etc.

def particle_match(tnu, other):
  return particles.get(tnu)

def find_particles(here, other):
  particles = {}
  count = [0]
  def log(tnu, message):
    if count[0] < 0:
      if debug:
       print("# fp(%s): %s" % (cl.get_unique(tnu), message))
      count[0] += 1
  def subanalyze(tnu, other):
    diffs = diff.no_diffs
    log(tnu, "subanalyze")
    if cl.get_accepted(tnu):
      print("# ** Child %s of %s has an accepted name" %
            (cl.get_unique(tnu), cl.get_unique(cl.get_parent(tnu))))
      return False
    found_match = False
    for inf in get_children(tnu):
      log(tnu, "child %s" % inf)
      if subanalyze(inf, other):
        found_match = True
      else:
        diffs = diff.note_dropped_children(diffs)
    if found_match:    # Some descendant is a particle
      return True
    candidate = choose_best_match(matches_to_accepted(tnu, other))
    if candidate:
      rematch = choose_best_match(matches_to_accepted(candidate.cod, here))
      if rematch:
        if rematch.cod == tnu:
          if cl.get_accepted(candidate.cod):
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

# ---------- Matches to accepted nodes

# Three components: synonym-or-self o intensional o synonym-of-or-self
# from_accepted_articulations o intensional_matches o [to_accepted_articulation]

# Matches where the target is accepted.

def matches_to_accepted(tnu, other):
  return weak_matches_to_accepted(tnu, other)

def weak_matches_to_accepted(tnu, other):
  matches = strong_matches_to_accepted(tnu, other)
  syns = synonyms_locally(tnu)
  assert not syns or syns[0].dom == tnu
  matches = matches + [art.compose(syn, match)
                       for syn in syns
                       for match in strong_matches_to_accepted(syn.cod, other)
                       if art.composable(syn, match)]
  assert is_matches(matches)
  return collapse_matches(matches)

def strong_matches_to_accepted(tnu, other):
  # (a) Go intensional to accepted
  # (b) Go intensional from synonym to synonym, then up to accepted
  # (c) Go from synonym to accepted, then intensional to accepted

  matches = intensional_matches(tnu, other)
  assert is_matches(matches)
  matches = [match_to_accepted(match) for match in matches]
  assert is_matches(matches)
  if cl.get_accepted(tnu):
    accept = has_accepted_locally(tnu)
    assert accept.dom == tnu
    more = [art.compose(accept, match) for match in intensional_matches(accept.cod, other)]
    matches = matches + more
  return collapse_matches(matches)

# Convert match a->b to a->b' where b' is accepted (maybe b=b')

def match_to_accepted(m):
  if cl.get_accepted(m.cod):
    return art.compose(m, has_accepted_locally(m.cod))
  else:
    return m


# Intensional matches by name (no synonym following)
# This ought to be cached I think?

def intensional_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  hits = cl.get_tnus_with_value(other,
                                cl.canonical_name,
                                cl.get_name(node))
  if shared_idspace:
    id_hit = cl.get_tnu_with_id(other, cl.get_tnu_id(node))
    if id_hit and not id_hit in hits:
      hits = hits + [id_hit]
  matches = [art.intensional(node, hit, 0) for hit in hits]
  return sort_by_badness(matches)   # ? is sort needed ?

# ---------- UNMATCHED

# Unmatched
# Find nodes in B that are not mutually matched to nodes in A

# A B_node that is not the best_match of any A_node ...

def analyze_unmatched(A, B):
  graft_points = {}
  def process(B_tnu):
    if not good_match(B_tnu, A) and not cl.get_accepted(B_tnu):
      point = get_graft_point(B_tnu, A)    # in A
      if point:
        graft_points[B_tnu] = point    # in A
    else:
      for child in cl.get_inferiors(B_tnu):
        process(child)
  for root in cl.get_roots(B):
    process(root)
  return cl.invert_dict(graft_points)

def get_graftees(A_node):
  return (grafts.get(A_node) or [])

def get_graft_point(B_tnu, A = None):
  B_parent = cl.get_parent(B_tnu)
  if B_parent:
    B_parent_match = good_match(B_parent, A)
    if B_parent_match:
      return B_parent_match.cod
  return None

# ---------- Within-checklist articulations

# Handy for composing paths.

# Synonym-or-self = to-accepted
# Matches that involve going from an accepted name to a synonym are weak

# An accepted name can have many synonyms

def synonyms_locally(tnu):
  if cl.get_accepted(tnu):
    return []
  else:
    hases = [has_accepted_locally(syn) for syn in cl.get_synonyms(tnu)]
    return collapse_matches([art.reverse(ar) for ar in hases if ar])

# A synonym has only one accepted name

def has_accepted_locally(maybe_syn):   # goes from synonym to accepted
  assert maybe_syn > 0
  if cl.get_accepted(maybe_syn):
    accepted = cl.get_accepted(maybe_syn)
    if accepted:
      if accepted != maybe_syn:
        return art.synonymy(maybe_syn, accepted)
      else:
        print("# ** Shouldn't happen", cl.get_unique(maybe_syn))
        return art.identity(maybe_syn)
    else:
      print("# ** Synonym %s has no accepted name" % cl.get_unique(maybe_syn))
      return None
  else:
    return None

# ---------- Best match among a set with common domain

# Assumes already sorted by art.badness

def choose_best_match(arts):     # => art
  assert is_matches(arts)
  if len(arts) == 0: return None
  arts = choose_best_matches(sort_by_badness(collapse_matches(arts)))
  b = arts[0]
  if len(arts) == 1: return b
  print("# ** Multiple least-bad matches. Need to find more tie-breakers.")
  print("   %s -> %s" %
        (cl.get_unique(b.dom),
         [cl.get_unique(a.cod) for a in arts]))
  print(', '.join(["%s" % ar.relation.name for ar in arts]))
  return None

# There can be multiple best matches

def choose_best_matches(sorted_matches):
  if len(sorted_matches) == 0:
    return []
  else:
    badness = art.badness(sorted_matches[0])
    matches = []
    for match in sorted_matches:
      if art.badness(match) == badness:
        matches.append(match)
      else:
        break
    return matches

def sort_by_badness(arts):
  return sorted(arts, key=art.badness)

# Random

def get_children(node):
  return [node
          for node in cl.get_children(node) + cl.get_synonyms(node)
          if not cl.get_accepted(node)]

def parent_changed(match):
  parent = cl.get_parent(match.dom)
  coparent = cl.get_parent(match.cod)
  if not parent and not coparent:
    return False
  if not parent or not coparent:
    return True
  other = cl.get_checklist(coparent)
  match = good_match(parent, other)
  if not match: return True
  return match.cod != coparent

# ---------- Utility: collapsing a set of matches

# Reduce a set of articulations so that all the codomains are
# different

def collapse_matches(arts):
  if len(arts) <= 1: return arts
  arts = sorted(arts, key=art.conjoin_sort_key)
  previous = None
  matches = []
  for ar in arts:
    assert ar
    if previous and art.conjoinable(previous, ar):
      previous = art.conjoin(previous, ar)
    else:
      if previous:
        matches.append(previous)
      previous = ar
  if previous:
    matches.append(previous)
  assert len(matches) <= len(arts)
  return matches

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
  shared_idspace = args.idspace
  main(args.left, args.left_tag, args.right, args.right_tag, args.out, args.format)

