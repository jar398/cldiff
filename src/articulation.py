
# Articulations must support:
#   The basics: source and destination taxa, and an RCC-5 (or similar) relation
#   A comment explaining the reason the articulation might be true
#   Reversal
#   Composition
#   Disjunction ??

import collections
import relation as rel

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

def identity(node):
  return art(node, node, rel.eq)

def compose(p, q):
  assert p.cod == q.dom
  if is_identity(p): return q
  if is_identity(q): return p
  return Articulation(p.dom,
                      q.cod,
                      rel.compose(p.relation, q.relation))

def conjoin(p, q):
  assert p.dom == q.dom
  if p.cod == q.cod and rel.is_variant(p.relation, q.relation):
    re = rel.conjoin(p.relation, q.relation)
    return Articulation(p.dom, p.cod, re)
  else:
    return None

def get_comment(art):
  return art.relation.name

def reverse(art):
  return Articulation(art.cod, art.dom, rel.reverse(art.relation))

def is_identity(art):
  return art.dom == art.dom and art.relation == rel.eq
