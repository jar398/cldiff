import checklist as cl
import articulation as art
import dribble

# Temporary hack for experimenting with poorly formed EOL checklists
EOL = False

# ---------- 'Intensional' matches to accepted nodes

# The source node ('tnu') may be accepted or a synonym.

def best_intensional_match_map(A, B, captured):
  best = {}
  for node in captured:
    ar = captured[node]
    if is_variant(ar.relation, rel.eq):
      best[ar.dom] = ar
      best[ar.cod] = art.reverse(ar)
  def process(here, there):
    for node in here.get_all_nodes():
      if cl.is_accepted(node) and node not in best:
        ar = best_intensional_match(node, there)
        if ar:
          assert cl.is_accepted(ar.cod)
          best[node] = ar
  process(A, B)
  process(B, A)
  dribble.log("%s best matches" % len(best))
  return best

# Three components: synonym-or-self o direct o synonym-of-or-self
# from_accepted_articulations o direct_matches o [to_accepted_articulation]

def best_intensional_match(tnu, other):
  matches = intensional_matches(tnu, other)
  return choose_best_match(matches)

# Matches where the target is accepted.  This is going away soon...

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
  accepted = cl.to_accepted(maybe_syn)
  if accepted != maybe_syn:
    return art.synonymy(maybe_syn, accepted)
  else:
    return None

# ---------- Best match among a set with common domain

def choose_best_match(arts):     # => art
  assert is_matches(arts)
  if len(arts) == 0: return None
  arts = skim_best_matches(arts)
  b = arts[0]
  if len(arts) == 1: return b
  dribble.log("** Multiple least-bad matches. Need to find tie-breakers.")
  dribble.log("   %s -> %s" %
              (cl.get_unique(b.dom),
               [cl.get_unique(a.cod) for a in arts]))
  return None

# There can be multiple best matches

def skim_best_matches(arts):
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

