
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
                         ['dom', 'cod', 'relation'])

def art(dom, cod, re):
  assert dom > 0
  assert cod > 0
  assert re
  assert re.name
  return Articulation(dom, cod, re)

def express(ar):
  return "%s %s %s" % (cl.get_unique(ar.dom),
                       ar.relation.name,
                       cl.get_unique(ar.cod))

def identity(node):
  return art(node, node, rel.eq)

def compose(p, q):
  if not composable(p, q):
    print("** Not composable:\n  %s &\n  %s" %
          (express(p), express(q)))
    assert False
  if is_identity(p): return q
  if is_identity(q): return p
  return Articulation(p.dom, q.cod, rel.compose(p.relation, q.relation))

def composable(p, q):
  return (p.cod == q.dom and
          rel.composable(p.relation, q.relation))

def conjoin(p, q):
  if not conjoinable(p, q):
    print("** Not conjoinable:\n  %s &\n  %s" %
          (express(p), express(q)))
    assert False
  re = rel.conjoin(p.relation, q.relation)
  return Articulation(p.dom, p.cod, re)

def conjoinable(p, q):
  return (p.dom == q.dom and
          q.cod == q.cod and
          rel.conjoinable(p.relation, q.relation))

def get_comment(art):
  return art.relation.name

def reverse(art):
  return Articulation(art.cod, art.dom, rel.reverse(art.relation))

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
         cl.get_rank(ar.cod))
