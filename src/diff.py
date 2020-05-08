debug = False

import sys, csv
import argparse

import checklist as cl
import relation as rel
import articulation as art
import eulerx

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

# Partners list for reporting (A->B articulations)

def good_candidate(tnu, other):  # for children sort order
  return best_match(good_candidates(tnu, other))

def good_candidates(tnu, other):
  matches = inverse_good_candidates.get(tnu) or []
  if len(matches) == 0:
    investigate = good_match(tnu, other)
    if investigate: matches = [investigate]
  return [match for match in matches if not cl.get_accepted(match.cod)]


# Fill the cache

def assign_matches(here, other):
  global inverse_good_candidates
  good_match_map = {}
  def process(tnu):
    m = good_match(tnu, other)  # Fill the cache
    if m: good_match_map[tnu] = m    # here -> other
    for child in get_children(tnu):
      process(child)
  for root in cl.get_roots(here):
    process(root)
  inverse_good_candidates = invert_dict_by_cod(good_match_map)

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

# ---------- OVERALL

def start(A, B):
  global particles
  global cross_mrcas

  particles = find_particles(A, B)
  print ("# Number of particles:", len(particles)>>1)

  cross_mrcas = analyze_cross_mrcas(A, B)
  print ("# Number of cross-mrcas:", len(cross_mrcas))

  assign_matches(B, A)
  print ("# Number of besties:", len(inverse_good_candidates))

  analyze_unmatched(A, B)
  print ("# Number of grafts: %s\n" % len(grafts))

  # return finish_alignment(B, A)

# ---------- UNMATCHED

# Unmatched
# Find nodes in B that are not mutually matched to nodes in A

# A B_node that is not the best_match of any A_node

def analyze_unmatched(A, B):
  global grafts
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
  grafts = cl.invert_dict(graft_points)

def get_graftees(A_node):
  return (grafts.get(A_node) or [])

def get_graft_point(B_tnu, A = None):
  B_parent = cl.get_parent(B_tnu)
  if B_parent:
    B_parent_match = good_match(B_parent, A)
    if B_parent_match:
      return B_parent_match.cod
  return None

# ---------- Alignments

def finish_alignment(B, A):
  alignment = {}
  def process(node):
    m = good_match(node, A)
    if m:
      alignment[node] = m
    for child in cl.get_children(node):
      process(child)
    for synonym in cl.get_synonyms(node):
      process(synonym)
  for root in cl.get_roots(B):
    process(root)
  return alignment

# ---------- One-sided best match

alignment = {}

def good_match(node, other = None):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in alignment:
    return alignment[node]

  assert other
  matches = good_matches(node, other)
  match = best_match(matches)
  if match:
    if cl.get_accepted(match.cod) and not cl.get_accepted(match.dom):
      print("# ** Match is supposed to be terminal:\n  %s" %
            (art.express(match)))
    # Happens way often in GBIF
    if False and not cl.is_accepted(match.cod):
      print("# ** Match has taxonomic status %s\n  %s" %
            (cl.get_value(match.cod, cl.taxonomic_status_field),
             art.express(match)))

  alignment[node] = match
  return match

good_matches_cache = {}

def good_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  if node in good_matches_cache: return good_matches_cache[node]

  exties = extensional_matches(node, other)
  exties = [match_to_accepted(m) for m in exties]
  nameys = matches_to_accepted(node, other)
  matches = collapse_matches(exties + nameys)

  good_matches_cache[node] = matches
  return matches

# ---------- EXTENSIONALITY

# Starting with one match, extend to a set of matches by adding
# extensionally equivalent ancestors

def extensional_matches(tnu, other):
  assert cl.get_checklist(tnu) != other

  match = extensional_match(tnu, other)    # Single extensional match
  if not match: return []
  matches = [match]

  if rel.is_variant(match.relation, rel.extension_eq):
    # Scan upwards looking for nodes that return to us...
    scan = match.cod    # tnu -> partner
    here = cl.get_checklist(scan)
    while True:
      scan = cl.get_superior(scan)    # in other
      if not scan: break
      rev = extensional_match(scan, here)
      # rev: other -> here
      if not rev: break
      if not rel.is_variant(rev.relation, rel.extension_eq): break
      if rev.cod != tnu: break
      matches.append(rel.reverse(rev))
  matches.reverse()    # hmm. for choose ?
  return matches

# Guaranteed invertible, except for monotypic node chains

def extensional_match(tnu, other):
  assert tnu > 0
  here = cl.get_checklist(tnu)
  match = particle_match(tnu, other)
  if match:
    return match
  else:
    (partner, _, _) = cross_mrca(tnu, other)     # TNU in other checklist
    if not partner:
      return None
    match = bridge(tnu, partner, how_related_extensionally(tnu, partner))
    return match_to_accepted(match)

def how_related_extensionally(tnu, partner):
  here = cl.get_checklist(tnu)
  i = particle_match(tnu, partner)
  if i:
    return i.relation
  else:
    (back, _, _) = cross_mrca(partner, here)
    if back == None:                # shouldn't happen
      print("# ** No return cross-MRCA\n  %s -> %s -> ??" %\
            (cl.get_unique(tnu), cl.get_unique(partner)))
      assert False
      return rel.extension_disjoint
    elif cl.are_disjoint(tnu, back):
      return rel.extension_disjoint
    elif tnu == back:
      return rel.eq
    elif cl.mrca(tnu, back) == tnu:
      # Monotypic: back < tnu.
      return rel.monotypic_in
    else:
      for sub in get_children(partner):
        (back, _, _) = cross_mrca(sub, here)
        if back != None:
          if cl.mrca(tnu, back) == tnu and cross_disjoint(tnu, partner):
            return rel.extension_conflict
      return rel.extension_lt

# ---------- Determine disjointness across checklists

def cross_disjoint(tnu, partner):
  assert tnu > 0
  assert partner > 0
  assert cl.get_checklist(tnu) != cl.get_checklist(partner)

  (back, _, _) = cross_mrca(partner, cl.get_checklist(tnu))
  if back == None:
    return True
  assert back > 0
  assert cl.get_checklist(back) == cl.get_checklist(tnu)
  if cl.are_disjoint(tnu, back):
    return True
  for child in get_children(partner):
    assert child > 0
    if not cross_disjoint(tnu, child):
      print("# Test %s ! %s failed because not ! %s" %\
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
    print("# ** No cross-MRCA set for %s" % cl.get_unique(tnu))
  return (cl.forest, 1, 0)    # shouldn't happen

# Returns (the-mrca, number-unmatched)

def analyze_cross_mrcas(A, B):
  cross_mrcas = {cl.forest: (cl.forest, 0, 0)}
  def half_analyze_cross_mrcas(checklist, other, checkp):
    def subanalyze_cross_mrcas(tnu, other):
      pmatch = particle_match(tnu, other)
      drops = 0     # These get upgraded
      adds = 0
      m = 0      # is the identity for mrca
      if pmatch:
        assert pmatch.dom == tnu
        if rel.is_variant(pmatch.relation, rel.eq):
          m = pmatch.cod
          drops = len(get_children(pmatch.dom))
          adds = len(get_children(pmatch.cod))
          if debug:
           print("#   particle(%s) = %s" % (cl.get_unique(tnu), cl.get_unique(pmatch.cod)))
        else:
          # Can't happen
          print("# ** Non-eq articulation from particle_match\n  %s" +
                art.express(pmatch))
          assert False
      else:
        children = get_children(tnu)
        if children:
          for child in children:
            if debug:
              print("# Child %s of %s" % (cl.get_unique(child), cl.get_unique(tnu)))
            (m2, d2, a2) = subanalyze_cross_mrcas(child, other)
            if m2 != None:
              drops += d2
              adds += a2
              if debug:
                print("#  Folding %s into %s" % (cl.get_unique(m2), cl.get_unique(m)))
              m = cl.mrca(m, m2)
              if debug:
                print("#   -> %s" % cl.get_unique(m))
        else:
          drops = 1
      if debug:
        print("#   cm(%s) = (%s, %s, %s)" % \
              (cl.get_unique(tnu), cl.get_unique(m), drops, adds))
      result = (m, drops, adds)
      cross_mrcas[tnu] = result
      if checkp:
        (back, _, _) = cross_mrcas[m]
        if not back in cross_mrcas:
          print("# ** No return cross-MRCA for %s -> %s -> %s" %\
                (cl.get_unique(tnu), cl.get_unique(m), cl.get_unique(back)))
          assert False
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
    if found_match:    # Some descendant is an particle
      return True
    candidate = best_match(matches_to_accepted(tnu, other))
    if candidate:
      rematch = best_match(matches_to_accepted(candidate.cod, here))
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
                                cl.canonical_name_field,
                                cl.get_name(node))
  matches = [bridge(node, hit, rel.same_name) for hit in hits]
  assert is_matches(matches) and True
  if shared_idspace:
    id_hit = cl.get_tnu_with_id(other, cl.get_tnu_id(node))
  else: id_hit = None
  if id_hit and not id_hit in hits:
    matches = [bridge(node, id_hit, rel.same_id)] + matches
  assert is_matches(matches)
  rank_matches = [bridge(node, match.cod, rel.same_rank)
                  for match in matches
                  if cl.get_nominal_rank(match.cod) == cl.get_nominal_rank(node)]
  return sort_by_badness(collapse_matches(matches + rank_matches))    # Should be cached.

# ---------- Within-checklist articulations

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
        status = rel.synonym_relation(maybe_syn_status(maybe_syn))
        return art.art(maybe_syn, accepted, status)
      else:
        print("# ** Shouldn't happen", cl.get_unique(maybe_syn))
        return art.identity(maybe_syn)
    else:
      print("# ** Synonym %s has no accepted name" % cl.get_unique(maybe_syn))
      return None
  else:
    return None

def maybe_syn_status(synonym):
  return cl.get_value(synonym, cl.nomenclatural_status_field) or \
         cl.get_value(synonym, cl.taxonomic_status_field) or \
         "synonym"

# ---------- Best match among a set with common domain

# Assumes already sorted by art.badness

def best_match(arts):     # => art
  assert is_matches(arts)
  if len(arts) == 0: return None
  arts = best_matches(sort_by_badness(collapse_matches(arts)))
  b = arts[0]
  if len(arts) == 1: return b
  print("# ** Multiple least-bad matches. Need to find more tie-breakers.")
  print("   %s -> %s" %
        (cl.get_unique(b.dom),
         [cl.get_unique(a.cod) for a in arts]))
  print(', '.join(["%s" % ar.relation.name for ar in arts]))
  return None

# There can be multiple best matches

def best_matches(sorted_matches):
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

def bridge(dom, cod, re):
  assert cl.get_checklist(dom) != cl.get_checklist(cod)
  return art.art(dom, cod, re)

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

