# Relations
# Characterized by P(B|A), P(A|B)

import sys
import collections

Rcc5 = \
  collections.namedtuple('Rcc5', ['b_given_a', 'a_given_b'])

Relation = \
  collections.namedtuple('Relation',
                         ['rcc5', 'goodness',
                          'path_length', 'name', 'revname'])

strong = 1
weak = 100

def _relation(rcc5, goodness,
              name, revname = None, path_length = weak):
  assert path_length >= 0
  assert name
  assert isinstance(rcc5, Rcc5)
  if revname == None:
    if rcc5.b_given_a == rcc5.a_given_b:
      revname = name
    else:
      revname = name + "-of"
  return Relation(rcc5, goodness, path_length, name, revname)

identity = _relation(Rcc5(1, 1), 0, 'identical', 'identical', 0)

def variant(re, goodness, name, revname = None):
  assert name or revname
  return _relation(re.rcc5, goodness, name,
                   revname, re.path_length)

def is_variant(rel1, rel2):
  return rel1.rcc5 == rel2.rcc5

def justify(re, bit, name, revname = None):
  return _relation(re.rcc5, re.goodness | bit,
                   name, revname, re.path_length)

def adjacently(re, name, revname = None):
  return _relation(re.rcc5, re.goodness,
                   name, revname, strong)

def reverse(re):
  rre = _relation(Rcc5(re.rcc5.a_given_b, re.rcc5.b_given_a), re.goodness,
                  re.revname, re.name, re.path_length)
  assert rre
  return rre

def compose(rel1, rel2):
  goodness = rel1.goodness & rel2.goodness
  name    = _compose_names(rel1.name, rel2.name)
  revname = _compose_names(rel2.revname, rel1.revname)
  len = rel1.path_length + rel2.path_length
  return _relation(compose_rcc5(rel1.rcc5, rel2.rcc5), goodness,
                   name, revname, len)
  
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
                   _conjoin_names(rel1.revname, rel2.revname),
                   min(rel1.path_length, rel2.path_length))

def _conjoin_names(name1, name2):
  if name1 == name2: return name1
  elif name1.endswith(" ...}"):
    return name1                # Discard lesser reasons
  elif name1.endswith("}"):
    return name1[0:-1] + " ...}"
  elif name2.startswith("{"):
    return "{%s ...}" % name1
  elif name2.startswith("["):
    return "{%s ...}" % name1
  else:
    return "{%s & %s}" % (name1, name2)

def conjoinable(rel1, rel2):
  return is_variant(rel1, rel2)

# re1 < re2 in sort order iff sort_key(re1) < sort_key(re2)

def sort_key(re):
  return (rcc5_key(re),    # distinguish < from >
          -1-re.goodness,
          re.path_length)

def rcc5_key(re):
  assert re.name
  return (-(re.rcc5.a_given_b + re.rcc5.b_given_a),
          re.rcc5.b_given_a - re.rcc5.a_given_b)

def better(re1, re2):
  return sort_key(re1) < sort_key(re2)

# Goodness represented as bit manipulation

def bit(b): return (1 << b)

# Unjustified

noinfo_species = Rcc5(-1, -1)
noinfo    = _relation(noinfo_species,   0, 'noinfo')
disjoint  = _relation(Rcc5(0,   0),   0, '!')
conflict  = _relation(Rcc5(0.5, 0.5), 0, '><') 
lt        = _relation(Rcc5(1,   0.5), 0, '<', '>')
eq        = _relation(Rcc5(1,   1),   0, '=')
gt        = reverse(lt)

child_of  = adjacently(lt, 'child of', 'parent of')
parent_of = reverse(child_of)
sibling_of = adjacently(disjoint, 'sibling of')  # ?

# Intensional half-matches ?

# Intensional matches

same_rank          = justify(eq, bit(0), "same-rank")
same_parent        = justify(eq, bit(1), "same-parent")
has_vernacular     = justify(eq, bit(2), "vernacular")
has_synonym        = justify(eq, bit(3), "synonym")
same_name          = justify(eq, bit(4), "same-name")
same_id            = justify(eq, bit(5), "same-id")
same_name_and_id = conjoin(same_name, same_id)

# Individual half-matches?
# Individual matches?
# Cross-mrca half-matches ??
# Cross-mrca matches ???

# Justified extensionally (i.e. assuming nonentities don't matter)

extensionally = bit(12)

extension_disjoint  = justify(disjoint, extensionally, "x !")
extension_conflict  = justify(conflict, extensionally, "x ><") 
extension_lt        = justify(lt, extensionally, "x <", "x >")
extension_eq        = justify(eq, extensionally, "x !")
extension_gt = reverse(extension_lt)

# Monotypy is a tough one.  The extension is the same, but the
# intension differs.

monotypic_over = justify(gt, bit(14), "monotypic-over", "monotypic-in")
monotypic_in   = reverse(monotypic_over)

def rcc5_name(re):
  if is_variant(re, eq): return eq.name
  if is_variant(re, lt): return lt.name
  if is_variant(re, gt): return gt.name
  if is_variant(re, conflict): return conflict.name
  if is_variant(re, disjoint): return disjoint.name
  else: return re.name

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
  re = synonym_relations.get(nomenclatural_status)
  if re: return re
  print("Unrecognized nomenclatural status: %s" % status)
  return reverse(synonym)

# These relations go from synonym to accepted (the "has x" form)

synonym_relations = {}

def declare_synonym_relations():

  def b(nstatus, rcc5 = eq, name = None, revname = None, relation = has_synonym):
    if name == None: name = "has-" + nstatus.replace(" ", "-")
    if revname == None: revname = nstatus.replace(" ", "-") + "-of"
    synonym_relations[nstatus] = \
      reverse(variant(rcc5, relation.goodness,
                      name, revname))

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
  b("genbank common name", relation = has_vernacular)    # at most one per node
  b("common name", relation = has_vernacular)

  b("includes", rcc5=gt, name="part-of", revname="included-in")
  b("in-part",  rcc5=lt, name="included-in", revname="part-of")  # part of a polyphyly

declare_synonym_relations()

