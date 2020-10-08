import checklist as cl
import articulation as art
import dribble

# Temporary hack for experimenting with poorly formed EOL checklists
EOL = False

# ---------- Particles

# A particle is mutual =-articulation of accepted nodes deriving from
# sameness of 'intrinsic' node properties: name, id, rank, parent,
# etc.

def find_particles(here, other):
  clear_cache()
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
    for inf in cl.get_children(tnu):
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
  print ("# Number of particles:", len(particles)>>1)
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
  matches = art.direct_matches(tnu, other)
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
        file=dribble.dribble_file)
  print("   %s -> %s" %
        (cl.get_unique(b.dom),
         [cl.get_unique(a.cod) for a in arts]),
        file=dribble.dribble_file)
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

def clear_cache():
  global intensional_matches_cache
  intensional_matches_cache = {}


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

