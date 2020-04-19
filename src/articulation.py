
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
                         ['dom', 'cod', 'rel', 'badness', 'comment'])

def art(dom, cod, re, badness, comment):
  assert dom
  assert cod
  assert badness >= 0
  return Articulation(dom, cod, re, badness, comment)

def reverse(art):
  return Articulation(art.cod, art.dom, rel.reverse(art.rel), art.badness, art.comment)

def compose(p, q):
  assert p.cod == q.dom
  if p.comment == q.comment:
    comment = p.comment
  elif p.comment == None:
    comment = q.comment
  elif q.comment == None:
    comment = p.comment
  else:
    comment = "%s; %s" % (p.comment, q.comment)
  return Articulation(p.dom,
                      q.cod,
                      rel.compose(p.rel, q.rel),
                      max(p.badness, q.badness),
                      comment)
