# Relations
# Characterized by P(B|A), P(A|B)

import collections

Relation = \
  collections.namedtuple('Relation', ['b_given_a', 'a_given_b', 'name'])

defined_relations = {}
def get_relation(b_given_a, a_given_b, name = None):
  key = (b_given_a, a_given_b)
  if key in defined_relations:
    return defined_relations[key]
  if name:
    rel = Relation(b_given_a, a_given_b, name)
    defined_relations[key] = rel
    return rel
  print("unrecognized relationship", b_given_a, a_given_b)
  return None

def reverse(rel, name = None):
  return get_relation(rel.a_given_b, rel.b_given_a, name)

def compose(rel1, rel2, name = None):
  b_given_a = min(rel1.b_given_a, rel2.b_given_a)
  a_given_b = min(rel1.a_given_b, rel2.a_given_b)
  return get_relation(b_given_a, a_given_b, name)
  
def fuzz(rel):
  return get_relation(min(rel.b_given_a, 0.9),
                      min(rel.a_given_b, 0.9))

# RCC5: find a representation that makes composition possible

# compose = < > with min
# but the others don't compose

eq        = get_relation(1, 1, '=')
lt        = get_relation(1, 0.5, '<')
gt        = reverse(lt, '>')
conflict  = get_relation(0.5, 0.5, '⟂')    # careful with this
disjoint  = get_relation(0, 0, '||')       # ∥ careful with this

# Non-RCC5 options

similar  = get_relation(0.9, 0.9, '~=')  # ≈
le       = get_relation(0.9, 0.5, '<=')  # ≤
ge       = reverse(le, '>=')  # ≥
intersect = get_relation(0.5, 0.5, '∩')

assert reverse(eq) == eq
assert compose(eq, eq) == eq
assert compose(lt, lt) == lt
assert compose(disjoint, gt) == disjoint
assert compose(eq, similar) == similar
# assert compose(lt, le) == lt  - doesn't work
