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
  b_given_a = min(rel1.b_given_a, rel2.b_given_a)
  a_given_b = min(rel1.a_given_b, rel2.a_given_b)
  if a_given_b == 0 and b_given_a == 0:
    pass                        # Composing with no_info
  elif rel1.b_given_a < 0.8 and rel2.a_given_b < 0.8:
    print("losing information",
          rel1.name, rel1.b_given_a,
          rel2.a_given_b, rel2.name,
          file=sys.stderr)
    return compose(compose(rel1, no_info), rel2)
  badness = max(rel1.badness, rel2.badness)
  name = "%s; %s" % (rel1.name, rel2.name)
  revname = "%s; %s" % (rel2.revname, rel1.revname)
  return Relation(b_given_a, a_given_b, badness,
                  name, revname)
  
# RCC5: find a representation that makes composition possible

goodness = 0.01

eq        = get_relation(1, 1,     goodness, '=')
lt        = get_relation(1, 0.5,   goodness, '<', '>')
gt        = reverse(lt)
conflict  = get_relation(0.5, 0.5, goodness, '⟂') 
disjoint  = get_relation(0.1, 0.1, goodness, '||')

# Non-RCC5 options

no_info   = get_relation(0,   0,   0.05, '?')
le        = get_relation(1,   0.7, 0.02, '<=', '>=')
intersect = get_relation(0.3, 0.3, 0.04, '∩')

child     = variant(lt, 0, 'parent', 'child')

def self_tests():
  assert reverse(eq) == eq
  assert compose(eq, eq) == eq
  assert compose(lt, lt) == lt
  assert compose(disjoint, gt) == disjoint
  assert compose(eq, similar) == similar
