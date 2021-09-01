# Articulations must support:
#   The basics: source and destination taxa (nodes), and an RCC-5 (or similar) 
#     relation between them
#   An explanation for why the articulation might be true
#   Reversal
#   Composition
#   Disjunction ??

import collections
import relation as rel
import checklist as cl
import property
import changes
import dribble

# Articulations

# reason and factors are mutually exclusive.  reason is only for
# non-composed articulations.

Articulation = \
  collections.namedtuple('Articulation',
                         ['dom', 'cod', 'relation', 'factors',
                          'reason', 'revreason', 'diff'])

def _articulation(dom, cod, re,
                  reason = None, revreason = None, factors = None):
  assert dom > 0
  assert cod > 0
  assert re
  assert re.name
  dif = changes.all_diffs
  if cl.is_accepted(dom) and cl.is_accepted(cod):
    dif = changes.differences(dom, cod)
  assert reason or factors
  if reason and revreason == None: revreason = reason + " of"
  ar = Articulation(dom, cod, re, factors, reason, revreason, dif)
  return ar

def express(ar):
  if ar == None:
    return "none"
  else:
    return "%s %s %s" % (cl.get_unique(ar.dom),
                         ar.relation.name,
                         cl.get_unique(ar.cod))

def compose(p, q):
  if not composable(p, q):
    print("** Not composable:\n  %s &\n  %s" %
          (express(p), express(q)))
    assert False
  if is_identity(p): return q
  if is_identity(q): return p
  return _articulation(p.dom,
                       q.cod,
                       rel.compose(p.relation, q.relation),
                       reason = None,
                       factors = (p.factors or [p]) + (q.factors or [q]))

def reason(p):
  if p.factors:
    return "+".join(map(reason, p.factors))
  else:
    return p.reason

def composable(p, q):
  return (p.cod == q.dom and
          rel.composable(p.relation, q.relation))

def conjoin(p, q):
  if not conjoinable(p, q):
    print("** Not conjoinable:\n  %s &\n  %s" %
          (express(p), express(q)))
    assert False
  return p                      # ???

def conjoinable(p, q):
  return (p.dom == q.dom and
          p.cod == q.cod and
          rel.conjoinable(p.relation, q.relation))

def get_comment(ar):
  return ar.relation.name

def reverse(ar):
  if ar.factors:
    f = list(reversed(ar.factors))
  else:
    f = None
  return _articulation(ar.cod, ar.dom, rel.reverse(ar.relation),
                       reason = ar.reason,
                       revreason = ar.revreason,
                       factors = f)

def is_identity(ar):
  return ar.dom == ar.cod and ar.relation == rel.eq

def inverses(ar1, ar2):
  return (ar1.cod == ar2.dom and 
          ar1.dom == ar2.cod and 
          rel.inverses(ar1.relation, ar2.relation))

# Foo.  Phase out

def set_relation(ar, re):      # re = rel.eq
  return _articulation(ar.dom, ar.cod, re, ar.reason, ar.revreason, ar.factors)

def change_relation(ar, re, reason, revreason = reason):  # re = rel.gt
  assert reason
  return compose(set_relation(ar, re),
                 _articulation(ar.cod, ar.cod, rel.eq, reason, revreason))

# ---------- Synonymy relationship within one tree

def synonymy(synonym, accepted):
  assert synonym > 0
  assert accepted > 0
  status = (cl.get_nomenclatural_status(synonym) or \
            cl.get_taxonomic_status(synonym) or \
            "synonym")
  return _articulation(synonym, accepted, rel.matches,
                       reason = status, revreason = status + "-of")

# Some NCBI nomenclatural status values.  This is just a comment.

def declare_synonym_relations(b):
  b("homotypic synonym")    # GBIF
  b("authority")
  b("scientific name")        # (actually canonical) exactly one per node
  b("equivalent name")        # synonym but not nomenclaturally
  b("misspelling")
  b("unpublished name")    # non-code synonym
  b("genbank synonym")        # at most one per node; first among equals
  b("anamorph")
  b("genbank anamorph")    # at most one per node
  b("teleomorph")
  b("acronym")
  b("blast name")             # large well-known taxa
  b("genbank acronym")      # at most one per node
  b("BOLD id")

  # More dubious
  synonym = b("synonym")
  b("heterotypic synonym")      # GBIF
  b("misnomer")
  b("type material")
  b("merged id", revname = "split id")    # ?
  b("accepted")    # EOL
  b("invalid")     # EOL

  # Really dubious
  b("genbank common name")    # at most one per node
  b("common name")

  b("includes", relation=rel.gt, name="part-of", revname="included-in")
  b("in-part",  relation=rel.lt, name="included-in", revname="part-of")  # part of a polyphyly
  b("proparte synonym", relation=rel.lt)

# ---------- Different kinds of articulation

def intensional(dom, cod, reason):
  return bridge(dom, cod, rel.matches, reason)

def extensional(dom, cod, re, reason, revreason = None):
  ar = bridge(dom, cod, re, reason, revreason)
  if dribble.watch(dom):
    dribble.log("# Extensional articulation %s" % art.express(ar))
  return ar

def bridge(dom, cod, re, reason, revreason = None):
  assert cl.get_checklist(dom) != cl.get_checklist(cod)
  assert reason
  return _articulation(dom, cod, re, reason=reason, revreason=revreason)

# Intensional matches by name (no synonym following)
# This ought to be cached I think?

def direct_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  seen = []
  arts = []
  for prop in [cl.ncbi_id,
               cl.eol_page_id,
               cl.scientific_name,
               cl.canonical_name,
               cl.gbif_id]:
    val = cl.get_value(node, prop)
    if val != None:
      more = cl.get_nodes_with_value(other, prop, val)
      for hit in more:
        if hit and not hit in seen:
          seen.append(hit)
          arts.append(intensional(node, hit, prop.pet_name))
  return arts

# ---------- Utility: collapsing a set of matches

# Reduce a set of articulations grouped first by RCC5 relation, then
# within each group, collapsed (conjoined) so that the codomains are
# all different.

def collapse_matches(arts):
  if len(arts) <= 1: return arts
  arts = sorted(arts, key=conjoin_sort_key)
  previous = None
  matches = []
  for ar in arts:
    if not previous:
      previous = ar
    elif conjoinable(ar, previous):
      previous = conjoin(ar, previous)
    else:
      matches.append(previous)
      previous = None
  if previous:
    matches.append(previous)
  assert len(matches) <= len(arts)
  return matches

# This one is for deduplication (grouping by codomain)

def conjoin_sort_key(ar):
  assert ar.dom
  return (rel.sort_key(ar.relation),
          ar.cod)

# ---------- This one is for tie breaking (when codomains differ)

# Less-bad articulations first.

def badness(ar):
  (drop, change, add) = ar.diff
  return(
         rel.sort_key(ar.relation),     # '=' sorts earliest
         # Changes are bad
         # Low-bit changes are better than high-bit changes
         # Additions don't matter
         change,
         # Using synonym is bad, using two is worse
         len(ar.factors) if ar.factors else 1,
         drop,
         # Added fields are benign
         # Dropped fields are so-so
         # What is this about?
         cl.get_mutex(ar.cod))

def sort_matches(arts):
  return sorted(arts, key=badness)

# -----

def proclaim_eq(draft, ar):
  assert ar.relation == rel.eq
  proclaim(draft, ar)
  proclaim(draft, reverse(ar))

def proclaim(draft, ar):
  assert cl.get_checklist(ar.dom) != cl.get_checklist(ar.cod)
  assert ar.relation != rel.matches
  p = proclaimable(draft, ar)
  if p == MEH:
    if dribble.watch(ar.dom):
      dribble.log("# Meh: %s" % express(ar))
    pass
  elif p:
    if dribble.watch(ar.dom):
      dribble.log("# storing")
    draft[ar.dom] = ar
  else:
    print("** Not OK to replace %s\n   with %s" %
          (express(draft.get(ar.dom)), express(ar)))
    assert False

  if dribble.watch(ar.dom):
    dribble.log("# Proclaim %s\n  yields %s" %
                (express(ar), express(draft.get(ar.dom))))

def proclaimable(draft, ar):
  if dribble.watch(ar.dom): print("# Proclaimable? %s" % express(ar))
  before = draft.get(ar.dom)
  if not before:
    if dribble.watch(ar.dom): print("# New: %s" % express(ar))
    return True                 # Improvement

  elif before.relation == rel.matches:
    if dribble.watch(ar.dom): print("# Upgrade ~: %s" % express(ar))
    return True                 # Improvement

  elif ar.cod == before.cod:
    if ar.relation == before.relation:
      # Update reason
      if dribble.watch(ar.dom): print("# Change reason only: %s" % express(ar))
      return True                 # Improvement
    else:
      # Relations differ, so this is an incompatible change.
      if dribble.watch(ar.dom): print("# Wrong RCC5: %s" % express(ar))
      return False                # Inconsistent

  # Codomains are different, relations may or may not be different

  elif ar.relation == rel.eq and before.relation == rel.eq:
    if dribble.watch(ar.dom): print("# UNA: %s" % express(ar))
    return False                # Inconsistent (unique name assumption)
  elif before.relation == rel.gt:
    if dribble.watch(ar.dom): print("# Upgrade: %s" % express(ar))
    return True

  else:
    if dribble.watch(ar.dom):
      print("# Meh: before %s,\n  after        %s" %
            (express(before), express(ar)))
    return MEH

MEH = 2
