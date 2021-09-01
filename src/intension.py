import checklist as cl
import articulation as art
import relation as rel
import dribble

# Temporary hack for experimenting with poorly formed EOL checklists
EOL = False

# ---------- 'Intensional' matches to accepted nodes

# The source node ('node') may be accepted or a synonym.

def best_intensional_match_map(A, B):
  best = {}
  def process(here, there):
    for node in here.get_all_nodes():
      if cl.is_accepted(node) and not node in best:
        ar = best_intensional_match(node, there)
        if dribble.watch(node):
          dribble.log("# Best: %s" % art.express(ar))
        if ar:
          assert ar.dom == node
          assert cl.is_accepted(ar.cod)
          best[node] = ar   # art.proclaim(best, ar)
  process(A, B)
  process(B, A)
  dribble.log("%s best matches" % len(best))
  return best

# Three components: synonym-or-self o direct o synonym-of-or-self
# from_accepted_articulations o direct_matches o [to_accepted_articulation]

def best_intensional_match(node, other):
  matches = intensional_matches(node, other)
  return choose_best_match(matches)

# Matches where the target is accepted.  This is going away soon...

def intensional_matches(node, other):
  matches = weak_intensional_matches(node, other)
  matches = art.collapse_matches(matches)
  return matches

# Like strong matches, but adds matches via source synonyms

def weak_intensional_matches(node, other):
  matches = strong_intensional_matches(node, other)
  syns = synonyms_locally(node)
  assert not syns or syns[0].dom == node
  matches = matches + [art.compose(syn, match)
                       for syn in syns
                       for match in strong_intensional_matches(syn.cod, other)
                       if art.composable(syn, match)]
  assert is_matches(matches)
  return matches

def strong_intensional_matches(node, other):
  matches = art.direct_matches(node, other)
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

def synonyms_locally(node):
  if cl.is_accepted(node):
    hases = [has_accepted_locally(syn) for syn in cl.get_synonyms(node)]
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

# ---- Find ad hoc splits and merges based on multimatches.

def intensional_proposal(best, A, B):
  matches = tipward(best, A, B)

  result = {}

  incoming = index_by_target(matches)

  # Suppose `node` x comes from the A checklist, and there is a split
  # such that x matches multiple nodes y1, y2 in the B checklist.

  # Modify the relation for all approximate-match nodes.
  for y in incoming:            # Many x's, one y
    arts = incoming[y]

    # Canonical.  back.cod will be among the incoming, by construction.
    back = matches.get(y)       # back : y -> x
    if not back: continue
    x0 = back.cod               # Back match y -> x -> y

    revarts = incoming[x0]

    if len(arts) > 1:     # multiple x's
      if len(revarts) > 1:
        art.proclaim_eq(result, art.set_relation(back, rel.eq))
        dribble.log("** Tangle:\n   %s\n   %s" %
                    ("\n   ".join(map(art.express, arts)),
                    ("\n   ".join(map(art.express, revarts)))))
      else:

        # OK.  We're going to just throw away all non-sibling matches.

        x_rent = cl.get_parent(x0)
        sibs = [ar for ar in arts if cl.get_parent(ar.dom) == x_rent]
        # e.g. ar: x2 -> y
        # Don't even try to do anything with N->M node tangles.
        if len(sibs) == 1:
          art.proclaim_eq(result, art.set_relation(back, rel.eq))
        else:
          for sib in sibs:
            # change x2 ~ y to x2 < y
            ar = art.change_relation(sib, rel.lt, "merge", "split")
            art.proclaim(result, ar)
          # y < parent(x2)
          # art.proclaim(result, art.bridge(y, x_rent, rel.lt, "split", "merge"))
          # Report!
          dribble.log("# Split/lump %s < %s < %s" %
                      (" + ".join(map(lambda e:cl.get_unique(e.dom), sibs)),
                       cl.get_unique(y),
                       cl.get_unique(x_rent)))

    elif len(revarts) > 1:   # multiple y's
      pass
    else:
      # n.b. arts[0] is reverse of back
      art.proclaim_eq(result, art.set_relation(back, rel.eq))

  return result

def index_by_target(al):
  incoming = {}
  for node in al:            # y node
    ar = al[node]      # ar : y -> x
    y = ar.cod
    if y in incoming:
      incoming[y].append(ar) # x -> (y->x, y2->x)
    else:
      incoming[y] = [ar]     # x -> (y->x)
  return incoming


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

# ---------- Tipwards

# A particle is a mutual =-articulation of tipward accepted nodes
# deriving from sameness of 'intrinsic' node properties: name, id,
# rank, parent, etc.

# This function returns a partial map from nodes to articulations.

# Filter out internal nodes (those having a matched descendant)

def tipward(amap, A, B):
  tw = {}
  def filter(node):
    watch = dribble.watch(node)
    found_match = None
    for child in cl.get_children(node):
      ar = filter(child)
      if ar:
        found_match = ar
    if found_match:    # Some descendant is a particle
      if watch: dribble.log("# %s: descendant matches, not keeping: %s" %
                            (cl.get_unique(node), art.express(found_match)))
      return found_match
    else:
      ar = amap.get(node)
      if ar:
        tw[node] = ar
        if watch:
          dribble.log("# %s is a tipward match, keeping: %s" %
                      (cl.get_unique(node), art.express(ar)))
      else:
        if watch: dribble.log("# %s is unmatched" % cl.get_unique(node))
      return ar
  for root in cl.get_roots(A): filter(root)
  for root in cl.get_roots(B): filter(root)
  return tw

