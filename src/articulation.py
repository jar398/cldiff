
# Articulations must support:
#   The basics: source and destination taxa, and an RCC-5 (or similar) relation
#   A 'badness' score so that they can be compared and sorted
#   A comment explaining the reason the articulation might be true
#   Reversal
#   Composition
#   Disjunction ??

import collections
import relation as rel

# Articulations

Articulation = \
  collections.namedtuple('Articulation',
                         ['dom', 'cod', 'relation', 'relations'])

def art(dom, cod, re):
  assert dom > 0
  assert cod > 0
  assert re
  assert re.name
  return Articulation(dom, cod, re, [re])

def identity(node):
  return art(node, node, rel.eq)

def compose(p, q):
  assert p.cod == q.dom
  if is_identity(p): return q
  if is_identity(q): return p
  return Articulation(p.dom,
                      q.cod,
                      rel.compose(p.relation, q.relation),
                      p.relations + q.relations)

def get_comment(art):
  return "; ".join([re.name for re in art.relations])

def reverse(art):
  return Articulation(art.cod, art.dom, rel.reverse(art.relation),
                      [rel.reverse(re) for re in art.relations[::-1]])

def is_identity(art):
  return art.dom == art.dom and art.relation == rel.eq
