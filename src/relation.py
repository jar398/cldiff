# Relations
# Characterized by P(B|A), P(A|B)

import sys
import collections

# RCC5 relations

Relation = \
  collections.namedtuple('Rcc5', ['b_given_a', 'a_given_b', 'name', 'revname'])

relations_by_params = {}
relations_by_name = {}

for (ba, ab, name, revname) in \
    [(0,   0,   '!',  '!'),
     (0.5, 0.5, '><', '><'),
     (1,   0.5, '<',  '>'),
     (1,   1,   '=',  '='),
     (0.9, 0.9, '~', '~')]:     # approximate
  re  = Relation(ba, ab, name, revname)
  relations_by_params[(ba, ab)] = re
  relations_by_name[name] = re
  if ab != ba:
    rev = Relation(ab, ba, revname, name)
    relations_by_params[(ab, ba)] = rev
    relations_by_name[revname] = rev

def _relation(b_given_a, a_given_b):
  return relations_by_params[(b_given_a, a_given_b)]

# RCC-5 relations

disjoint    = relations_by_name['!']
conflict    = relations_by_name['><'] 
lt          = relations_by_name['<']
gt          = relations_by_name['>']
eq          = relations_by_name['=']
matches     = relations_by_name['~']

# Reverse relation

def reverse(re):
  return _relation(re.a_given_b, re.b_given_a)

def inverses(r1, r2):
  return (r1.b_given_a == r2.a_given_b and
          r1.a_given_b == r2.b_given_a)

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
  if r1 == matches and r2 == matches: return r1
  b_given_a = min(r1.b_given_a, r2.b_given_a)
  a_given_b = min(r1.a_given_b, r2.a_given_b)
  if r1.b_given_a < 1 and r2.a_given_b < 1:
    return False
  return _relation(b_given_a, a_given_b)

# Main thing is to sort '=' ahead of everything else

def sort_key(re):
  return (-re.b_given_a, -re.a_given_b) # distinguish < from >

def self_tests():
  assert reverse(eq) == eq
  assert compose(eq, eq) == eq
  assert compose(lt, lt) == lt
  assert compose(disjoint, gt) == disjoint

if __name__ == '__main__':
  self_tests()
