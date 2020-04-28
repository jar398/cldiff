# Relations
# Characterized by P(B|A), P(A|B)

import sys
import collections

Relation = \
  collections.namedtuple('Relation',
                         ['b_given_a', 'a_given_b', 'goodness', 'name', 'revname'])

def _relation(b_given_a, a_given_b, goodness, name = None, revname = None):
  assert name or revname
  if name == None and revname != None:
    if b_given_a == a_given_b:
      name = revname
    else:
      name = revname + " of"
  if revname == None:
    if b_given_a == a_given_b:
      revname = name
    else:
      revname = name + " of"
  return Relation(b_given_a, a_given_b, goodness, name, revname)

def variant(re, goodness, name, revname = None):
  assert name or revname
  return _relation(re.b_given_a, re.a_given_b, goodness, name, revname)

def is_variant(rel1, rel2):
  return rel1.b_given_a == rel2.b_given_a and \
         rel1.a_given_b == rel2.a_given_b

def reverse(re):
  rre = _relation(re.a_given_b, re.b_given_a, re.goodness,
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
  goodness = rel1.goodness & rel2.goodness
  name = "%s; %s" % (rel1.name, rel2.name)
  revname = "%s; %s" % (rel2.revname, rel1.revname)
  return _relation(b_given_a, a_given_b, goodness,
                              name, revname)
  
def conjoin(rel1, rel2, name = None, revname = None):
  assert is_variant(rel1, rel2)
  if name == None:
    if rel1.name == rel2.name: name = rel1.name
    else: name = "%s & %s" % (rel1.name, rel2.name)
  if revname == None:
    if rel1.revname == rel2.revname: name = rel1.revname
    else: revname = "%s & %s" % (rel1.revname, rel2.revname)
  return variant(rel1,
                 rel1.goodness | rel2.goodness,
                 name,
                 revname)

def sort_key(re):
  return -1 - re.goodness

# Goodness represented as bit manipulation

def bit(b): return (1 << b)

no_info   = _relation(1, 1, bit(16)-1, '?')

# Direct matches
same_rank          = variant(no_info, bit(0), "rank=")
same_parent        = variant(no_info, bit(1), "parent=")
vernacular         = variant(no_info, bit(2), "vernacular")
synonym            = variant(no_info, bit(3), "synonym")
same_name          = variant(no_info, bit(4), "name=")
same_id            = variant(no_info, bit(5), "id=")

same_name_and_id = conjoin(same_name, same_id)

# RCC5
topo_disjoint  = _relation(0,   0,   bit(9), 'tips ||')
topo_conflict  = _relation(0.5, 0.5, bit(10), 'tips ⟂') 
topo_lt        = _relation(1,   0.5, bit(11), 'tips<', 'tips>')
topo_eq        = _relation(1,   1,   bit(12), 'tips=')           # =
topo_gt        = reverse(topo_lt)

same_checklist     = variant(no_info, bit(13), 'checklist=')

def rcc5_name(re):
  if is_variant(re, eq): return eq.name
  if is_variant(re, lt): return lt.name
  if is_variant(re, gt): return gt.name
  if is_variant(re, conflict): return conflict.name
  if is_variant(re, disjoint): return disjoint.name
  else: return re.name

def rcc5(topo, name, revname = None):
  return variant(topo, (same_checklist.goodness | topo.goodness), name, revname)

disjoint = rcc5(topo_disjoint, '||')
conflict = rcc5(topo_conflict, '⟂')
lt       = rcc5(topo_lt,       '<', '>')
eq       = rcc5(topo_eq,       '=')
gt       = reverse(lt)

fringe_and_name = conjoin(eq, same_name)

def self_tests():
  assert reverse(eq) == eq
  assert compose(eq, eq) == eq
  assert compose(lt, lt) == lt
  assert compose(disjoint, gt) == disjoint
  assert compose(eq, similar) == similar

# -------------------- Synonyms 

def synonym_relation(nomenclatural_status):
  if nomenclatural_status == None:
    return synonym
  re = synonym_types.get(nomenclatural_status)
  if re: return re
  print("Unrecognized nomenclatural status: %s" % status)
  return synonym

synonym_types = {}

def declare_synonym_types():

  def b(nstatus, rcc5 = eq, name = None, relation = synonym):
    revname = nstatus
    synonym_types[nstatus] = variant(rcc5, relation.goodness, name, revname)

  b("homotypic synonym")    # GBIF
  b("authority")
  b("scientific name")        # (actually canonical) exactly one per node
  b("equivalent name")        # synonym but not nomenclaturally
  b("misspelling")
  b("unpublished name")    # non-code synonym
  b("genbank synonym")        # at most one per node; first among equals
  b("anamorph")
  b("genbank anamorph")    # at most one per node
  b("teleomorph")
  b("acronym")
  b("blast name")             # large well-known taxa
  b("genbank acronym")      # at most one per node

  # More dubious
  b("synonym")
  b("misnomer")
  b("type material")
  b("merged id")    # ?

  # Really dubious
  b("genbank common name", relation = vernacular)    # at most one per node
  b("common name", relation = vernacular)

  b("includes", rcc5=gt, name="part of")
  b("in-part",  rcc5=lt, name="included in")  # part of a polyphyly

declare_synonym_types()

