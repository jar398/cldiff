# Articulations must support:
#   The basics: source and destination taxa, and an RCC-5 (or similar) relation
#   A comment explaining the reason the articulation might be true
#   Reversal
#   Composition
#   Disjunction ??

import collections
import relation as rel
import checklist as cl

# Articulations

Articulation = \
  collections.namedtuple('Articulation',
                         ['dom', 'cod', 'relation', 'differences'])

def _articulation(dom, cod, re, differences):
  assert dom > 0
  assert cod > 0
  assert re
  assert re.name
  assert differences >= 0
  return Articulation(dom, cod, re, differences)

def express(ar):
  return "%s %s %s (%s)" % (cl.get_unique(ar.dom),
                            ar.relation.name,
                            cl.get_unique(ar.cod),
                            ar.differences)

def identity(node):
  return _articulation(node, node, rel.eq, 0)

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
                       p.differences + q.differences)

def composable(p, q):
  return (p.cod == q.dom and
          rel.composable(p.relation, q.relation))

def conjoin(p, q):
  if not conjoinable(p, q):
    print("** Not conjoinable:\n  %s &\n  %s" %
          (express(p), express(q)))
    assert False
  re = rel.conjoin(p.relation, q.relation)
  return Articulation(p.dom, p.cod, re,
                      min(p.differences, q.differences))

def conjoinable(p, q):
  return (p.dom == q.dom and
          q.cod == q.cod and
          rel.conjoinable(p.relation, q.relation))

def get_comment(art):
  return art.relation.name

def reverse(art):
  return Articulation(art.cod, art.dom, rel.reverse(art.relation), art.differences)

def is_identity(art):
  return art.dom == art.cod and art.relation == rel.eq

# This one is for deduplication (grouping by codomain)

def conjoin_sort_key(ar):
  assert ar.dom
  return (ar.cod,
          rel.rcc5_key(ar.relation))

# This one is for tie breaking (when codomains differ)

def badness(ar):
  return(rel.sort_key(ar.relation),
         ar.differences,
         cl.get_mutex(ar.cod))

# Synonymy relationship within one tree

def synonymy(syn, accepted):
  status = (cl.get_value(synonym, cl.nomenclatural_status_field) or \
            cl.get_value(synonym, cl.taxonomic_status_field) or \
            "synonym")
  re = rel.synonym_relation(status)
  return _articulation(syn, accepted, re, 1)

different_subtrees = 1 << cl.number_of_fields

def extensional(dom, cod, re, same_subtree):
  diffs = 0
  if not same_subtree:
    diffs = different_subtrees
  return bridge(dom, cod, re, diffs)

def intensional(dom, cod, diffs):
  re = rel.identical_subtrees if diffs == 0 else rel.similar_subtrees
  return bridge(dom, cod, re, diffs)

def bridge(dom, cod, re, diffs):
  assert cl.get_checklist(dom) != cl.get_checklist(cod)
  diffs |= cl.differences(dom, cod)
  return _articulation(dom, cod, re, diffs)

def cross_mrca(dom, cod, diffs):
  return _articulation(dom, cod, rel.le, diffs)
