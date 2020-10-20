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
import diff

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
  dif = diff.all_diffs
  if cl.is_accepted(dom) and cl.is_accepted(cod):
    dif = diff.differences(dom, cod)
  ar = Articulation(dom, cod, re, factors, reason, revreason, dif)
  return ar

def express(ar):
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
  return _articulation(ar.dom, ar.cod, re, ar.reason, ar.revreason)

def change_relation(ar, re, reason, revreason):  # re = rel.gt
  return compose(set_relation(ar, re),
                 _articulation(ar.cod, ar.cod, rel.eq, reason, revreason))

# ---------- Synonymy relationship within one tree

def synonymy(synonym, accepted):
  assert synonym > 0
  assert accepted > 0
  status = (cl.get_nomenclatural_status(synonym) or \
            cl.get_taxonomic_status(synonym) or \
            "synonym")
  re = synonym_relation(status)
  return _articulation(synonym, accepted, re,
                       reason="synonym")

# I don't understand this

def synonym_relation(nom_status):
  re = synonym_relations.get(nom_status)
  if re: return re
  print("Unrecognized nomenclatural status: %s" % nom_status)
  return rel.intensional

# These relations go from synonym to accepted (the "has x" form)
# TBD: Put these back into the articulation somehow

synonym_relations = {}

def declare_synonym_relations():

  def b(nstatus, relation = rel.intensional, name = None, revname = None):
    if False:
      if name == None: name = "has-" + nstatus.replace(" ", "-")
      if revname == None: revname = nstatus.replace(" ", "-") + "-of"
    re = rel.reverse(relation)
    synonym_relations[nstatus] = re
    return re

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

declare_synonym_relations()

# ---------- Different kinds of articulation

def intensional(dom, cod, how):
  return bridge(dom, cod, rel.intensional, how)

def cross_mrca(dom, cod):
  return bridge(dom, cod, rel.cross_mrca, "cross-mrca")

def extensional(dom, cod, re, reason):
  return bridge(dom, cod, re, reason)

def monotypic(dom, cod, re):
  return bridge(dom, cod, re, "monotypic")

def bridge(dom, cod, re, reason):
  assert cl.get_checklist(dom) != cl.get_checklist(cod)
  return _articulation(dom, cod, re, reason=reason)

# Intensional matches by name (no synonym following)
# This ought to be cached I think?

def direct_matches(node, other):
  assert node > 0
  assert cl.get_checklist(node) != other
  seen = []
  arts = []
  for prop in [cl.scientific_name,
               cl.canonical_name,
               cl.ncbi_id,
               cl.gbif_id,
               cl.eol_id]:
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
    elif conjoinable(previous, ar):
      previous = conjoin(previous, ar)
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
         drop,
         # Using synonym is bad, using two is worse
         len(ar.factors) if ar.factors else 1,
         # Added fields are benign
         # Dropped fields are so-so
         # What is this about?
         cl.get_mutex(ar.cod))

def sort_matches(arts):
  return sorted(arts, key=badness)
