debug = False

import sys, csv

import checklist as cl
import relation as rel
import articulation as art

# For each B-record, we choose an articulation with the closest
# A-record that it matches (preferably but not necessarily an '='
# articulation).

def align(B, A):
  particles = find_particles(B, A)
  print ("# Number of particles:", len(particles)>>1)
  cross_mrcas = analyze_cross_mrcas(B, A, particles)
  print ("# Number of cross-mrcas:", len(cross_mrcas))
  the_alignment = finish_alignment(B, A, cross_mrcas)
  print ("# Number of articulations in alignment:", len(the_alignment))
  grafts = analyze_unmatched(B, A, cross_mrcas)
  print ("# Number of grafts: %s\n" % len(grafts))
  return (the_alignment, grafts)

def finish_alignment(B, A, xmrcas):
  global cross_mrcas
  cross_mrcas = xmrcas
  alignment = {}
  def process(node):
    m = articulate(node, A, xmrcas)
    if m:
      alignment[node] = m
    for child in get_children(node):
      process(child)
  for root in cl.get_roots(B):
    process(root)
  return alignment

# ---------- One-sided best match

def articulate(node, other, xmrcas):     # B-node to A
  assert node > 0
  assert not cl.get_accepted(node)
  assert other
  assert cl.get_checklist(node) != other
  matches = good_matches(node, other)
  match = choose_best_match(matches)
  if match:
    if cl.get_accepted(match.cod) and not cl.get_accepted(match.dom):
      print("# ** Match is supposed to be terminal:\n  %s" %
            (art.express(match)))
  return match

good_matches_cache = {}

def good_matches(node, other):
  assert node > 0
  assert not cl.get_accepted(node)
  assert cl.get_checklist(node) != other
  if node in good_matches_cache: return good_matches_cache[node]

  if debug:
    print("# matching B-node", cl.get_unique(node))
  exties = extensional_matches(node, other)    #monotypic chain

  namies = matches_to_accepted(node, other)
  targets = [extie.cod for extie in exties]
  namies = [match for match in namies if match.cod in targets]

  if len(namies) > 0:
    matches = namies
  elif len(exties) > 0:
    if debug: print("# %s exties to %s" % (len(exties), cl.get_unique(node)))
    matches = exties
  else:
    if debug: print("# no extie")
    matches = matches_to_accepted(node, other)

  good_matches_cache[node] = matches
  return matches

# ---------- EXTENSIONALITY / by split

# Starting with one match, extend to a set of matches by adding
# 'monotypic' ancestors

def extensional_matches(tnu, other):       # B-node to A
  assert cl.get_checklist(tnu) != other

  match = extensional_match(tnu, other)    # Single B/A extensional match
  if not match:
    return []
  matches = [match]
  if debug: print("# %s has initial extensional match %s" %
                  (cl.get_unique(tnu), cl.get_unique(match.cod)))
  if rel.is_variant(match.relation, rel.eq):

    here = cl.get_checklist(tnu)
    anchor = cross_mrca(match.cod, here).cod

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
      back = cross_mrca(scan, here)
      if back.cod != anchor:
        if debug: print("# goes back to %s not %s, breaking out" %
                        (cl.get_unique(back.cod), cl.get_unique(anchor)))
        break
      if debug: print("# adding return match %s" %
                      (cl.get_unique(scan)))
      matches.append(art.extensional(tnu, scan, rel.same_particles))
    if debug: print("# %s matches %s nodes extensionally" %
                    (cl.get_unique(tnu), len(matches)))
  return matches

# Guaranteed invertible, except for monotypic node chains

def extensional_match(tnu, other):
  cross = cross_mrca(tnu, other)      # TNU in other checklist
  if not cross:
    return None
  partner = cross.cod
  here = cl.get_checklist(tnu)
  atcha = cross_mrca(partner, here)
  back = atcha.cod
  if cl.are_disjoint(tnu, back):
    re = rel.disjoint
  elif cl.mrca(tnu, back) == tnu:
    re = rel.same_particles
  else:
    re = rel.lt
    for sub in get_children(partner):
      cross_back = cross_mrca(sub, here)
      if cross_back != None:
        back = cross_back.cod
        if cl.mrca(tnu, back) == tnu and cross_disjoint(tnu, partner):
          re = rel.conflict
  re = extensionally(re)
  return art.extensional(tnu, partner, re)

def extensionally(re):
  return rel.justify(re, rel.extensionally, "extension " + rel.rcc5_name(re))

# ---------- Determine disjointness across checklists

def cross_disjoint(tnu, partner):
  assert tnu > 0
  assert partner > 0
  assert cl.get_checklist(tnu) != cl.get_checklist(partner)

  cross = cross_mrca(partner, cl.get_checklist(tnu))
  if cross == None:
    return True
  back = cross.cod
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

# Returns an accepted/accepted articulation

def cross_mrca(tnu, other):
  assert tnu > 0
  assert not cl.get_accepted(tnu)
  probe = cross_mrcas.get(tnu, 17)
  if probe != 17:
    return probe
  return None

def analyze_cross_mrcas(A, B, particles):
  cross_mrcas = {}
  def half_analyze_cross_mrcas(checklist, other, checkp):
    def subanalyze_cross_mrcas(tnu, other):
      result = particles.get(tnu)
      if result:
        assert result.dom == tnu
        assert not cl.get_accepted(result.cod)
        if debug:
          print("#   particle(%s) = %s" %\
                (cl.get_unique(tnu), cl.get_unique(result.cod)))
      else:
        children = get_children(tnu)
        if children:
          m = None      # is the identity for mrca
          for child in children:
            if debug:
              print("# Child %s of %s" % (cl.get_unique(child), cl.get_unique(tnu)))
            cross = subanalyze_cross_mrcas(child, other)
            if cross:
              m2 = cross.cod
              if debug:
                print("#  Folding %s into %s" % (cl.get_unique(m2), cl.get_unique(m)))
              m = cl.mrca(m, m2) if m != None else m2
              # ?????
              if debug:
                print("#   -> %s" % cl.get_unique(m))
          if m != None:
            result = art.cross_mrca(tnu, m)
            if debug:
              print("#   cm(%s) = %s)" % \
                    (cl.get_unique(tnu), cl.get_unique(m)))
      if result:
        cross_mrcas[tnu] = result
        if debug and checkp:
          m = result.cod
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
    if cl.get_accepted(tnu):
      print("# ** Child %s of %s has an accepted taxname" %
            (cl.get_unique(tnu), cl.get_unique(cl.get_parent(tnu))))
      return False
    found_match = False
    for inf in get_children(tnu):
      log(tnu, "child %s" % inf)
      if subanalyze(inf, other):
        found_match = True
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
  return art.collapse_matches(matches)

def strong_matches_to_accepted(tnu, other):
  # (a) Go intensional (by name) to accepted
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
  return art.collapse_matches(matches)

# Convert match a->b to a->b' where b' is accepted (maybe b=b')

def match_to_accepted(m):
  if cl.get_accepted(m.cod):
    return art.compose(m, has_accepted_locally(m.cod))
  else:
    return m

# Intensional matches by name (no synonym following)
# This ought to be cached I think?

def intensional_matches(node, other):
  hits = cl.get_similar_records(other, node, shared_idspace)
  return [art.intensional(node, hit) for hit in hits]

# ---------- UNMATCHED (grafts)

# Unmatched - TBD: ought to be folded into alignment
# (but we need to distinguish grafts from insertions)
# Find nodes in B (?) that are not mutually matched to nodes in A

# A B_node that is not the best_match of any A_node ...

def analyze_unmatched(B, A, xmrcas):
  graft_points = {}
  def process(B_tnu):
    if not cl.get_accepted(B_tnu) and not articulate(B_tnu, A, xmrcas):
      point = get_graft_point(B_tnu, A, xmrcas)    # in A
      if point:
        graft_points[B_tnu] = point    # in A
    else:
      for child in cl.get_inferiors(B_tnu):
        process(child)
  for root in cl.get_roots(B):
    process(root)
  return cl.invert_dict(graft_points)

def get_graft_point(B_tnu, A, xmrcas):
  B_parent = cl.get_parent(B_tnu)
  if B_parent:
    B_parent_match = articulate(B_parent, A, xmrcas)
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
    return art.collapse_matches([art.reverse(ar) for ar in hases if ar])

# A synonym has only one accepted name

def has_accepted_locally(maybe_syn):   # goes from synonym to accepted
  assert maybe_syn > 0
  accepted = cl.get_accepted(maybe_syn)
  if accepted:
    if accepted != maybe_syn:
      return art.synonymy(maybe_syn, accepted)
    else:
      print("# ** Shouldn't happen", cl.get_unique(maybe_syn))
      return art.identity(maybe_syn)
  else:
    return None

# ---------- Best match among a set with common domain

def choose_best_match(arts):     # => art
  assert is_matches(arts)
  if len(arts) == 0: return None
  arts = choose_best_matches(arts)
  b = arts[0]
  if len(arts) == 1: return b
  print("# ** Multiple least-bad matches. Need to find more tie-breakers.")
  print("   %s -> %s" %
        (cl.get_unique(b.dom),
         [cl.get_unique(a.cod) for a in arts]))
  print(', '.join(["%s" % ar.relation.name for ar in arts]))
  return None

# There can be multiple best matches

def choose_best_matches(arts):
  if len(arts) == 0:
    return []
  else:
    sorted_matches = art.sort_by_badness(art.collapse_matches(arts))
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
  return [node
          for node in cl.get_children(node) + cl.get_synonyms(node)
          if not cl.get_accepted(node)]

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
  shared_idspace = args.idspace # global
  main(args.left, args.left_tag, args.right, args.right_tag, args.out, args.format)

