debug = False

import sys, csv

import checklist as cl
import relation as rel
import articulation as art

# For each B-record, we choose an articulation with the closest
# A-record that it matches (preferably but not necessarily an '='
# articulation).

def align(B, A):
  global intensional_matches_cache
  intensional_matches_cache = {}
  particles = find_particles(B, A)
  print ("# Number of particles:", len(particles)>>1)

  cross_mrcas = analyze_cross_mrcas(B, A, particles)
  print ("# Number of cross-mrcas:", len(cross_mrcas))

  the_alignment = finish_alignment(B, A, cross_mrcas)
  print ("# Number of articulations in alignment:", len(the_alignment))

  return (the_alignment, cross_mrcas)

def finish_alignment(B, A, xmrcas):
  global cross_mrcas
  cross_mrcas = xmrcas
  alignment = {}
  def half_finish(B, A):
    def process(node):
      m = articulate(node, A, xmrcas)
      if m:
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

def articulate(node, other, xmrcas):     # B-node to A
  assert node > 0
  assert cl.is_accepted(node)
  assert other
  assert cl.get_checklist(node) != other

  if debug:
    print("# matching B-node", cl.get_unique(node))
  exties = extensional_matches(node, other)    #monotypic chain

  # Filter out intentionals that aren't extensional matches
  namies = intensional_matches(node, other)
  targets = [extie.cod for extie in exties]
  bothies = [match for match in namies if match.cod in targets]

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
      matches.append(art.extensional(tnu, scan, rel.eq))
    if debug: print("# %s matches %s nodes extensionally" %
                    (cl.get_unique(tnu), len(matches)))
  return matches

# Guaranteed invertible, except for monotypic node chains

def extensional_match(tnu, other):
  partner = cross_mrca(tnu)      # TNU in other checklist
  if not partner:
    return None
  here = cl.get_checklist(tnu)
  back = cross_mrca(partner)
  if cl.are_disjoint(tnu, back):
    re = rel.disjoint
  elif cl.mrca(tnu, back) == tnu:
    re = rel.eq
  else:
    re = rel.lt
    for sub in get_children(partner):
      cross_back = cross_mrca(sub)
      if cross_back != None:
        if cl.mrca(tnu, cross_back) == tnu and cross_disjoint(tnu, partner):
          re = rel.conflict
  return art.extensional(tnu, partner, re)

# ---------- Determine disjointness across checklists

def cross_disjoint(tnu, partner):
  assert tnu > 0
  assert partner > 0
  assert cl.get_checklist(tnu) != cl.get_checklist(partner)

  back = cross_mrca(partner)
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

def analyze_cross_mrcas(A, B, particles):
  cross_mrcas = {}
  def half_analyze_cross_mrcas(checklist, other, checkp):
    def subanalyze_cross_mrcas(tnu, other):
      result = None
      probe = particles.get(tnu)
      if probe:
        assert probe.dom == tnu
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
    candidate = best_intensional_match(tnu, other)
    if candidate:
      rematch = best_intensional_match(candidate.cod, here)
      if rematch:
        if rematch.cod == tnu:
          if not cl.is_accepted(candidate.cod):
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
  hits = cl.get_similar_records(other, node, shared_idspace)
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
  shared_idspace = args.idspace # global
  main(args.left, args.left_tag, args.right, args.right_tag, args.out, args.format)

