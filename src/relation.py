# Relations
# Characterized by P(B|A), P(A|B)

import sys
import collections

Rcc5 = \
  collections.namedtuple('Rcc5', ['b_given_a', 'a_given_b'])

Relation = \
  collections.namedtuple('Relation',
                         ['rcc5', 'goodness', 'name', 'revname'])

weak = 100

def _relation(rcc5, goodness, name, revname = None):
  assert name
  assert isinstance(rcc5, Rcc5)
  if revname == None:
    if rcc5.b_given_a == rcc5.a_given_b:
      revname = name
    else:
      revname = name + "-of"
  return Relation(rcc5, goodness, name, revname)

identity = _relation(Rcc5(1, 1), 0, 'identical', 'identical')

def variant(re, goodness, name, revname = None):
  assert name or revname
  return _relation(re.rcc5, goodness, name, revname)

def is_variant(rel1, rel2):
  return rel1.rcc5 == rel2.rcc5

def justify(re, bit, name, revname = None):
  return _relation(re.rcc5, re.goodness | bit,
                   name, revname)

def reverse(re):
  return _relation(Rcc5(re.rcc5.a_given_b, re.rcc5.b_given_a),
                   re.goodness, re.revname, re.name)

def compose(rel1, rel2):
  goodness = rel1.goodness & rel2.goodness
  name    = _compose_names(rel1.name, rel2.name)
  revname = _compose_names(rel2.revname, rel1.revname)
  return _relation(compose_rcc5(rel1.rcc5, rel2.rcc5), goodness,
                   name, revname)
  
def _compose_names(name1, name2):
  if name1.startswith("["): name1 = name1[1:-1]
  if name2.startswith("["): name2 = name2[1:-1]
  return "[%s → %s]" % (name1, name2)

def composable(rel1, rel2):
  return composable_rcc5s(rel1.rcc5, rel2.rcc5) != None

def composable_rcc5s(r1, r2):
  b_given_a = min(r1.b_given_a, r2.b_given_a)
  a_given_b = min(r1.a_given_b, r2.a_given_b)
  if r1.b_given_a < 1 and r2.a_given_b < 1:
    return None
  return Rcc5(b_given_a, a_given_b)

def compose_rcc5(r1, r2):
  rcc5 = composable_rcc5s(r1, r2)
  if not rcc5:
    print("** losing information", r1, r2,
          file=sys.stderr)
    return noinfo_rcc5
  return rcc5

def conjoin(rel1, rel2):
  assert conjoinable(rel1, rel2)
  return _relation(rel1.rcc5,
                   rel1.goodness | rel2.goodness,
                   _conjoin_names(rel1.name, rel2.name),
                   _conjoin_names(rel1.revname, rel2.revname))

def _conjoin_names(name1, name2):
  if name1 == name2: return name1
  if not name1: return name2
  if True: return name1

  # Clever stuff that now I think I don't want
  if not name2: return name1
  elif name1.endswith(" ...}"):
    return name1                # Discard lesser reasons
  elif name1.endswith("}"):
    return name1[0:-1] + " ...}"
  elif name2.startswith("{"):
    return "{%s & ...}" % name1
  elif name2.startswith("["):
    return "{%s & ...}" % name1
  else:
    return "{%s & %s}" % (name1, name2)

def conjoinable(rel1, rel2):
  return is_variant(rel1, rel2)

# re1 < re2 in sort order iff sort_key(re1) < sort_key(re2)

def sort_key(re):
  return (rcc5_key(re),    # distinguish < from >
          -1-re.goodness)

def rcc5_key(re):
  assert re.name
  return (-(re.rcc5.a_given_b + re.rcc5.b_given_a),
          re.rcc5.b_given_a - re.rcc5.a_given_b)

def better(re1, re2):
  return sort_key(re1) < sort_key(re2)

# Hack

def rcc5_name(re):
  if is_variant(re, eq): return eq.name
  if is_variant(re, lt): return lt.name
  if is_variant(re, gt): return gt.name
  if is_variant(re, conflict): return conflict.name
  if is_variant(re, disjoint): return disjoint.name
  else: return re.name

# Basic RCC-5 relations

noinfo_rcc5 = _relation(Rcc5(-1, -1),   0, 'noinfo')    # error
disjoint    = _relation(Rcc5(0,   0),   0, '!')
conflict    = _relation(Rcc5(0.5, 0.5), 0, '><') 
lt          = _relation(Rcc5(1,   0.5), 0, '<', '>')
eq          = _relation(Rcc5(1,   1),   0, '=')
gt          = reverse(lt)

# For cross-mrcas

le        = _relation(Rcc5(1,   0.7), 0, '<=', '>=')  # never composed

# ----------
# Goodness represented as bit manipulation.  Higher numbered bits are 'better'.

def bit(b): return (1 << b)

# Intensional matches

has_vernacular     = justify(eq, bit(2), "vernacular")
has_synonym        = justify(eq, bit(3), "synonym")

# Justified by particle set (split?)

extensionally = bit(10)
same_particle  = justify(eq, extensionally, "particle")
same_particles = justify(eq, extensionally, "same particles")

# Share all fields, transitively to all descendants

similar_subtrees   = justify(eq, bit(11), "similar subtrees")
identical_subtrees = justify(eq, bit(12), "identical subtrees")
identical          = justify(eq, bit(16), "= identically")

# -------------------- Synonyms 

def synonym_relation(nom_status):
  if nom_status == None:
    return synonym
  re = synonym_relations.get(nom_status)
  if re: return re
  print("Unrecognized nomenclatural status: %s" % nom_status)
  return reverse(synonym)

# These relations go from synonym to accepted (the "has x" form)

synonym_relations = {}

def declare_synonym_relations():
  global synonym

  def b(nstatus, rcc5 = eq, name = None, revname = None, relation = has_synonym):
    if name == None: name = "has-" + nstatus.replace(" ", "-")
    if revname == None: revname = nstatus.replace(" ", "-") + "-of"
    re = reverse(variant(rcc5, relation.goodness, name, revname))
    synonym_relations[nstatus] = re
    return re

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
  b("BOLD id")

  # More dubious
  synonym = b("synonym")
  b("heterotypic synonym")      # GBIF
  b("misnomer")
  b("type material")
  b("merged id", revname = "split id")    # ?

  # Really dubious
  b("genbank common name", relation = has_vernacular)    # at most one per node
  b("common name", relation = has_vernacular)

  b("includes", rcc5=gt, name="part-of", revname="included-in")
  b("in-part",  rcc5=lt, name="included-in", revname="part-of")  # part of a polyphyly
  b("proparte synonym", rcc5=lt)

declare_synonym_relations()

def self_tests():
  assert reverse(eq) == eq
  assert compose(eq, eq) == eq
  assert compose(lt, lt) == lt
  assert compose(disjoint, gt) == disjoint
  assert compose(eq, similar) == similar
