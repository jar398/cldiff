# Relations
# Characterized by P(B|A), P(A|B)

import sys
import collections

# RCC5 relations

Relation = \
  collections.namedtuple('Rcc5', ['b_given_a', 'a_given_b', 'name', 'revname'])

relations_by_params = {}
relations_by_name = {}

def declare_relation(ba, ab, name, revname):
  re  = Relation(ba, ab, name, revname)
  relations_by_params[(ba, ab)] = re
  relations_by_name[name] = re
  if ab != ba:
    rev = Relation(ab, ba, revname, name)
    relations_by_params[(ab, ba)] = rev
    relations_by_name[revname] = rev
  return re

def _relation(b_given_a, a_given_b):
  return relations_by_params[(b_given_a, a_given_b)]

# Reverse relation

def reverse(re):
  return _relation(re.a_given_b, re.b_given_a)

def inverses(r1, r2):
  return (r1.b_given_a == r2.a_given_b and
          r1.a_given_b == r2.b_given_a)

# RCC-5 relations

eq = declare_relation(1,   1,   '=',  '=')
refines = declare_relation(1,   0.6, '<1', '>1')    # potential child
refined_by = reverse(refines)
lt = declare_relation(1,   0.5, '<',  '>')    # not potential child (inconsistency somewhere)
gt = reverse(lt)
conflict = declare_relation(0.5, 0.5, '><', '><')
matches  = declare_relation(0.9, 0.9, '~', '~')     # approximate
disjoint = declare_relation(0,   0,   '!',  '!')

# Conjuction of two relations.  The 5 RCC-5 relations are mutually exclusive

def conjoinable(rel1, rel2):
  if rel1 == rel2:
    return rel1
  else:
    return False

def conjoin(rel1, rel2):
  re = conjoinable(rel1, rel2)
  assert re
  return re

# Composition of two relations

def compose(r1, r2):
  re = composable(r1, r2)
  if not re:
    print("** losing information", r1, r2,
          file=sys.stderr)
    assert False
  return re

def composable(r1, r2):
  if r2 == eq: return r1
  elif r1 == eq: return r2
  elif r1 == matches:
    return r1 if r2 == matches else None
  elif r1 == refines:
    # a <1 b and a < c imply a < c
    if is_lt_like(r2):
      return rel.lt
    else:
      return r2 if r2 == disjoint else None
  elif r1 == lt:
    # a < b and b < c imply a < c
    # a < b and b ! c imply a ! c
    return r2 if (is_lt_like(r2) or r2 == disjoint) else None
  elif r1 == gt:
    # a > b and b > c imply a > c
    return r2 if r1 == gt else None
  elif r1 == disjoint:
    # a ! b and b > c imply a ! c
    return r1 if r2 == gt else None
  else:                         # >< or <1
    return None

# Main thing is to sort '=' ahead of everything else

def sort_key(re):
  return (-re.b_given_a, -re.a_given_b) # distinguish < from >

def is_lt_like(re):
  return re == lt or re == refines

def self_tests():
  assert reverse(eq) == eq
  assert compose(eq, eq) == eq
  assert compose(lt, lt) == lt
  assert compose(disjoint, gt) == disjoint

if __name__ == '__main__':
  self_tests()
