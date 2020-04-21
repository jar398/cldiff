# Relations
# Characterized by P(B|A), P(A|B)

import sys
import collections

Relation = \
  collections.namedtuple('Relation',
                         ['b_given_a', 'a_given_b', 'badness', 'name', 'revname'])

defined_relations = {}
def get_relation(b_given_a, a_given_b, badness, name = None, revname = None):
  key = (b_given_a, a_given_b, badness)
  if key in defined_relations:
    return defined_relations[key]
  if name:
    if not revname:
      if b_given_a == a_given_b:
        revname = name
      else:
        revname = name + " of"
    def establish(b_given_a, a_given_b, badness, name, revname):
      re = Relation(b_given_a, a_given_b, badness, name, revname)
      defined_relations[key] = re
      return re
    establish(a_given_b, b_given_a, badness, revname, name)
    return establish(b_given_a, a_given_b, badness, name, revname)
  print("unrecognized relationship", name, b_given_a, a_given_b, badness, file=sys.stderr)
  assert False

def variant(re, badness, name = None, revname = None):
  assert re.name
  assert badness >=0
  return get_relation(re.b_given_a, re.a_given_b,
                      badness, name, revname)

def is_variant(rel1, rel2):
  return rel1.b_given_a == rel2.b_given_a and \
         rel1.a_given_b == rel2.a_given_b

def reverse(re):
  rre = get_relation(re.a_given_b, re.b_given_a, re.badness,
                     re.revname, re.name)
  assert rre
  return rre

def compose(rel1, rel2):
  if rel1.b_given_a < 0.8 and rel2.a_given_b < 0.8:
    print("cannot compose",
          rel1.name, rel1.b_given_a,
          rel2.a_given_b, rel2.name,
          file=sys.stderr)
    assert False
  b_given_a = min(rel1.b_given_a, rel2.b_given_a)
  a_given_b = min(rel1.a_given_b, rel2.a_given_b)
  badness = max(rel1.badness, rel2.badness)
  return get_relation(b_given_a, a_given_b, badness,
                      "%s; %s" % (rel1.name, rel2.name),
                      "%s; %s" % (rel2.revname, rel1.revname))
  
# RCC5: find a representation that makes composition possible

eq        = get_relation(1, 1,     0, '=')
lt        = get_relation(1, 0.5,   0, '<', '>')
gt        = reverse(lt)
conflict  = get_relation(0.5, 0.5, 0, '⟂') 
disjoint  = get_relation(0,   0,   0, '||')

# Not really RCC5

le = get_relation(0.9, 0.5, 1, '<=', '>=')

# Non-RCC5 options

intersect = get_relation(0.5, 0.5, 0, '∩')

child_parent = variant(lt, 0, 'parent', 'child')
parent_child = reverse(child_parent)

def self_tests():
  assert reverse(eq) == eq
  assert compose(eq, eq) == eq
  assert compose(lt, lt) == lt
  assert compose(disjoint, gt) == disjoint
  assert compose(eq, similar) == similar
