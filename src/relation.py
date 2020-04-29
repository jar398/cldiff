# Relations
# Characterized by P(B|A), P(A|B)

import sys
import collections

Relation = \
  collections.namedtuple('Relation',
                         ['b_given_a', 'a_given_b', 'badness', 'name', 'revname'])

defined_relations = {}
def intern_relation(b_given_a, a_given_b, badness, name = None, revname = None):
  key = (b_given_a, a_given_b, badness)
  if key in defined_relations:
    return defined_relations[key]
  revname = default_revname(name, revname, b_given_a, a_given_b)
  if name:
    def establish(b_given_a, a_given_b, badness, name, revname):
      re = Relation(b_given_a, a_given_b, badness, name, revname)
      defined_relations[key] = re
      return re
    establish(a_given_b, b_given_a, badness, revname, name)
    return establish(b_given_a, a_given_b, badness, name, revname)
  print("unrecognized relationship", name, b_given_a, a_given_b, badness, file=sys.stderr)
  assert False

def get_relation(b_given_a, a_given_b, badness, name = None, revname = None):
  assert badness >= 0
  revname = default_revname(name, revname, b_given_a, a_given_b)
  if badness == 0 or name == None:
    return intern_relation(b_given_a, a_given_b, badness, name, revname)
  else:
    return Relation(b_given_a, a_given_b, badness, name, revname)

def default_revname(name, revname, b_given_a, a_given_b):
  if name != None and revname == None:
    if b_given_a == a_given_b:
      revname = name
    else:
      revname = name + " of"
  return revname

def variant(re, badness, name = None, revname = None):
  assert re.name
  assert badness >=0
  return intern_relation(re.b_given_a, re.a_given_b, badness, name, revname)

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

eq        = get_relation(1, 1,     0, '=')
lt        = get_relation(1, 0.5,   0, '<', '>')
gt        = reverse(lt)
conflict  = get_relation(0.5, 0.5, 0, '><') 
disjoint  = get_relation(0.1, 0.1, 0, '!')

# Non-RCC5 options

le        = get_relation(1,   0.7, 10, '<=', '>=')
intersect = get_relation(0.3, 0.3, 10, '∩')
no_info   = get_relation(0,   0,   10, '?')

child     = variant(lt, 0, 'parent', 'child')

def self_tests():
  assert reverse(eq) == eq
  assert compose(eq, eq) == eq
  assert compose(lt, lt) == lt
  assert compose(disjoint, gt) == disjoint
  assert compose(eq, similar) == similar

# -------------------- Synonyms

def synonym_relation(nomenclatural_status):
  if nomenclatural_status == None: return badnesses["synonym"]
  status = badnesses.get(nomenclatural_status)
  if status:
    return status
  print("Unrecognized nomenclatural status: %s" % status)
  return badnesses["synonym"]

badnesses = {}
badness = 100

def declare_badnesses():
  def b(revname, re, name = None):
    assert re
    global badness
    name = name or (revname + " of")
    badnesses[revname] = variant(re, badness, name, revname)
    badness += 1

  b("homotypic synonym", eq)    # GBIF
  b("authority", eq)
  b("scientific name", eq)        # (actually canonical) exactly one per node
  b("equivalent name", eq)        # synonym but not nomenclaturally
  b("misspelling", eq)
  b("genbank synonym", eq)        # at most one per node; first among equals
  b("anamorph", eq)
  b("genbank anamorph", eq)    # at most one per node
  b("teleomorph", eq)
  b("unpublished name", eq)    # non-code synonym
  b("merged id", eq)
  b("acronym", eq)

  # above here: equivalence implied. below here: acc>=syn implied.
  # except in the case if 'in-part' which is acc<syn.

  b("synonym", eq)
  b("misnomer", eq)
  b("includes", gt, "included in")
  b("in-part", lt, "part of")      # this node is part of a polyphyly
  b("type material", eq)
  b("blast name", eq)             # large well-known taxa
  b("genbank common name", eq)    # at most one per node
  b("genbank acronym", eq)      # at most one per node
  b("common name", eq)

declare_badnesses()

