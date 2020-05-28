debug = False

import os, csv

import relation as rel
import rank
import chaitin
import property as prop

# ---------- Fields (columns) in taxon table

# Fields (columns in taxon table) are represented by selectors

field_selectors = []    # list of (Selector, position)
highest_priority = 0

def field_position(field):
  (_, col) = field
  return col

def field_selector(field):
  (sel, _) = field
  return sel

# Returns a field (selector, position)

def field_by_name(name):
  sel = prop.by_name(name)
  if sel in field_selectors:
    return (sel, field_selectors.index(sel))
  else:
    return None

# Particular fields of interest here (not all possible DwC fields)

def reserve(name):
  global highest_specificity
  sel = prop.by_name(name)
  assert sel
  field = (sel, len(field_selectors))
  field_selectors.append(sel)
  highest_specificity = max(highest_priority, sel.specificity)
  return field

checklist_field      = reserve("source")
nomenclatural_status = reserve("nomenclaturalStatus")
taxonomic_status     = reserve("taxonomicStatus")
taxon_rank           = reserve("taxonRank")
parent_taxon_id      = reserve("parentNameUsageID")
taxon_id             = reserve("taxonID")
accepted_taxon_id    = reserve("acceptedNameUsageID")
canonical_name       = reserve("canonicalName")
scientific_name      = reserve("scientificName")

number_of_fields = len(field_selectors)

# ---------- Taxon registry and taxa

# Registry = tnu uid -> value vector
#  where value vector is a vector that parallels the `field_selectors` list.
#  (it contains a pointer to the originating checklist.)

def is_tnu(x):
  # change to be the list itself?... no
  return isinstance(x, int) and x >= 0

forest_tnu = 0
registry = ["no forest record"]
sequence_numbers = {}

# Get the value of a field of a TNU record (via global registry)

def _get_record(uid):
  assert is_tnu(uid)
  return registry[uid]     # checklist is in record

def get_value(uid, field):
  if uid == forest_tnu: return None
  record = _get_record(uid)
  return record[field_position(field)]

def set_value(uid, field, value):
  record = _get_record(uid)
  record[field_position(field)] = value
         
def get_checklist(uid):
  assert is_tnu(uid)               # Forest has no checklist
  # get_value(uid, prop.source)
  if uid == forest_tnu: return None
  return get_value(uid, checklist_field)

def get_sequence_number(uid):
  assert uid                    # Forest has no checklist
  return sequence_numbers[uid]

def fresh_tnu():
  record = [None] * number_of_fields
  uid = len(registry)
  registry.append(record)
  return uid

def differences(tnu1, tnu2):  # mask
  r1 = _get_record(tnu1)
  r2 = _get_record(tnu2)
  comp = 0
  for position in range(number_of_fields):
    sel = field_selectors[position]
    v1 = r1[position]
    v2 = r2[position]
    if v1 and v2:
      if v1 != v2: comp |= 1 << position
  if debug:
    print("# Differences(%s, %s) = %o (octal)" %\
          (get_unique(tnu1), get_unique(tnu2), comp))
  return comp

# ---------- Checklists

class Checklist:
  def __init__(self, prefix, name):
    assert prefix
    self.prefix = prefix
    self.tnus = []
    self.indexes = [None] * number_of_fields
    self.metadata = None
    self.name = name
  def get_all_tnus(self):
    return self.tnus
  def tnu_count(self):
    return len(self.tnus)
  def new_tnu(self):
    tnu = fresh_tnu()
    set_value(tnu, checklist_field, self)
    self.tnus.append(tnu)
    return tnu
  def get_index(self, field):
    position = field_position(field)
    # create the index if necessary, on demand
    if self.indexes[position] == None:
      index = {}
      for tnu in self.tnus:
        # Check each record
        record = registry[tnu]
        value = record[position]
        if value:
          if value in index:
            index[value].append(tnu)
          else:
            index[value] = [tnu]
      self.indexes[position] = index
    return self.indexes[position]
  def populate_from_directory(self, indirs):
    self.populate_from_file(get_taxon_file_path(indir))
  def populate_from_file(self, inpath):
    if '(' in inpath:
      rows = chaitin.parse(inpath)
      self.populate_from_generator((row for row in rows))
    else:
      with open(inpath, "r") as infile:
        reader = csv.reader(infile, delimiter="\t", quotechar='\a', quoting=csv.QUOTE_NONE)
        self.populate_from_generator(reader)
  def populate_from_generator(self, row_generator):
    head = next(row_generator)
    guide = [None] * number_of_fields
    for position_in_row in range(len(head)):
      key = head[position_in_row]
      f = field_by_name(key)
      if f:
        guide[field_position(f)] = position_in_row
      else:
        print("# ** Column label %s isn't recognized" % key)
    if not "taxonID" in head:
      print("# ** Didn't find a column for taxonID")
      print("%s" % guide)
      assert False
    for row in row_generator:
      uid = self.new_tnu()
      record = _get_record(uid)
      # Fill in fields in record in order left to right
      for i in range(number_of_fields):
        sel = field_selectors[i]
        label = sel.pet_name
        position_in_row = guide[i]
        # Idiot python treats 0 as false
        if position_in_row != None:
          value = row[position_in_row]
          if value != None and value != '':
            record[i] = value
    self.assign_sequence_numbers()
  def assign_sequence_numbers(self):
    n = len(sequence_numbers)
    def process(tnu, n):
      sequence_numbers[tnu] = n
      n = n + 1
      for inf in get_inferiors(tnu):
        n = process(inf, n)
      return n
    for root in get_roots(self):
      n = process(root, n)

# A checklist's TNU list

def get_all_tnus(checklist):
  return checklist.get_all_tnus()

# Read a checklist from a file

def read_checklist(specifier, prefix, name):
  assert prefix
  checklist = Checklist(prefix, name)
  checklist.populate_from_file(specifier)
  return checklist

# Utility - copied from another file - really ought to be shared

def get_taxon_file_path(dwca_dir):
  for name in ["taxon.tsv",
               "Taxon.tsv",
               "taxon.tab",
               "Taxon.tab",
               "taxa.txt",
               "taxon.txt",
               "Taxon.txt"]:
    path = os.path.join(dwca_dir, name)
    if os.path.exists(path):
      return path
  raise ValueError("cannot find taxon file in this directory", dwca_dir)

# -------------------- indexing

# Get TNUs in checklist having a given value

canonical_empty_list = []

def get_tnus_with_value(checklist, field, value):
  return checklist.get_index(field).get(value, canonical_empty_list)

# Get unique (we hope) TNU possessing a given identifier

def get_tnu_id(tnu):
  id = get_value(tnu, taxon_id)
  if id == None:
    print("** No taxon id?? for %s" % _get_record(tnu))
    assert False
  return id

def get_tnu_with_id(checklist, id):
  tnus = checklist.get_index(taxon_id).get(id, None)
  if tnus:
    return tnus[0]
  else:
    return None
  return get_value(tnu, taxon_id)

# ----------------------------------------

# Logic for particular fields

# The following is for display purposes

def get_name(tnu):
  name = get_value(tnu, canonical_name)
  if name != None: return name
  name = get_value(tnu, scientific_name)
  if name != None: return name  
  return get_tnu_id(tnu)

def get_nominal_rank(tnu):
  if is_container(tnu): return None
  return get_value(tnu, taxon_rank)

# Unique name of the sort Nico likes

def get_spaceless(tnu):
  if tnu == None: return "none"
  assert is_tnu(tnu)
  if tnu == forest_tnu: return "forest"
  checklist = get_checklist(tnu)
  name = get_name(tnu)

  tnus_with_this_name = \
    get_tnus_with_value(checklist, canonical_name, name)
  if len(tnus_with_this_name) > 1:
    name = name + "#" + get_tnu_id(tnu)

  status = get_value(tnu, taxonomic_status)
  if status == "accepted" or not status:
    pass
  elif status == "synonym":
    name = "?" + name
  else:
    name = "%s (%s)" % (name, status)

  name = name.replace(" ", "_")
  return name

def get_unique(tnu):
  if tnu:
    return get_checklist(tnu).prefix + get_spaceless(tnu)
  else:
    return get_spaceless(tnu)

# Roots - accepted tnus without parents

def get_roots(checklist):
  roots = []
  for tnu in get_all_tnus(checklist):
    assert tnu > 0
    if not get_accepted(tnu) and not get_parent(tnu):
      roots.append(tnu)
  return roots

# Superior/inferior

def get_inferiors(tnu):
  assert is_tnu(tnu)
  assert tnu != forest_tnu
  return get_synonyms(tnu) + get_children(tnu)

# ----------
# Parent/children and accepted/synonyms

superiors_cache = {}

def get_superior(tnu):
  (accepted, parent) = get_superiors(tnu)
  return parent or accepted

def get_parent(tnu):
  assert is_tnu(tnu)
  assert tnu
  (_, parent) = get_superiors(tnu)
  return parent

def get_accepted(tnu):
  assert is_tnu(tnu)
  (accepted, _) = get_superiors(tnu)
  assert accepted != forest_tnu
  return accepted

# Returns (accepted, parent)

def get_superiors(tnu):
  assert is_tnu(tnu)
  probe = superiors_cache.get(tnu)
  if probe: return probe
  result = get_superiors_really(tnu)
  superiors_cache[tnu] = result
  if debug:
    print("# Cached gs(%s) = (%s, %s)" % \
          (get_unique(tnu), get_unique(result[0]), get_unique(result[1])))
  return result

def get_superiors_really(tnu):
  assert is_tnu(tnu)
  if tnu == forest_tnu: return (None, None)
  parent_id = get_value(tnu, parent_taxon_id)
  if parent_id == None:
    parent = forest_tnu
  else:
    parent_uid = get_tnu_with_id(get_checklist(tnu), parent_id)
    if parent_uid == None:
      parent = forest_tnu
    else:
      parent = to_accepted(parent_uid)
  
  # If id doesn't resolve, it's a root

  accepted_id = get_value(tnu, accepted_taxon_id)
  if accepted_id:
    accepted = to_accepted(get_tnu_with_id(get_checklist(tnu), accepted_id))
    if accepted:
      accepted_taxo_status = get_taxonomic_status(accepted)
      if accepted_taxo_status and accepted_taxo_status != "accepted":
        # Example: GBIF backbone:
        # C.?Hylobates sericus#8955169 -> C.Bunopithecus_sericus_(doubtful)
        print("# ** Putative accepted has wrong taxonomic status %s\n  %s -> %s" % \
              (accepted_taxo_status, get_unique(tnu), get_unique(accepted)))
      # Iterate???
    else:
      print("# ** Id %s for accepted of %s doesn't resolve" % \
            (accepted_id, get_unique(tnu)))
  else:
    n_status = get_taxonomic_status(tnu)
    if n_status == "synonym":
      print("# ** Node with status = 'synonym' has no accepted node\n  %s" % \
            (get_unique(tnu)))

    accepted = None

  if parent and accepted:
    parent_taxo_status = get_taxonomic_status(parent)
    if parent_taxo_status == "accepted" and accepted_taxo_status == "accepted":
      if False:
        print("# ** Taxon has both accepted and parent links\n  %s -> %s, %s" % \
              (get_unique(tnu), get_unique(accepted), get_unique(parent)))
        print("# ** Preferring parent to accepted")
      parent = None
    else:
      # This case is common in GBIF files
      accepted = None
  if parent != None: parent = to_accepted(parent)
  if accepted != None: accepted = to_accepted(accepted)
  return (accepted, parent)

def get_children(parent):
  return get_tnus_with_value(get_checklist(parent),
                             parent_taxon_id,
                             get_tnu_id(parent))

def to_accepted(tnu):
  assert is_tnu(tnu)
  probe = get_accepted(tnu)     # cached
  if probe:
    return probe
  else:
    return tnu

# ----------
# Accepted/synonyms

accepteds_cache = {}

def get_accepted_really(tnu):
  accepted = None
  if is_accepted(tnu):
    pass                        # This is normal
  else:
    accepted_id = get_value(tnu, accepted_taxon_id)
    if not accepted_id:
      pass                      # This is normal
    else:
      probe = get_tnu_with_id(get_checklist(tnu), accepted_id)
      if probe:
        if is_accepted(probe):
          accepted = probe      # Normal case
        else:
          print("# ** Accepted taxon %s for %s is not accepted" %
                (get_unique(probe), get_unique(tnu)))
      else:
        print("# ** Accepted id %s for %s doesn't resolve" %
              (accepted_id, get_unique(tnu)))

  return accepted

def get_synonyms(tnu):
  return [syn
          for syn in get_tnus_with_value(get_checklist(tnu),
                                         accepted_taxon_id,
                                         get_tnu_id(tnu))
          if not get_parent(tnu)]

def get_taxonomic_status(tnu):
  return get_value(tnu, taxonomic_status)

def is_accepted(tnu):
  return get_taxonomic_status(tnu) == "accepted"

def is_synonym(tnu):
  return get_taxonomic_status(tnu) == "synonym"

# ----------
# Totally general utilities from here down... I guess...

def invert_dict(d):
  inv = {}
  for (key, val) in d.items():
    if val in inv:
      inv[val].append(key)
    else:
      inv[val] = [key]
  return inv

# ---------- Hierarchy analyzers

def how_related(tnu1, tnu2):
  assert is_tnu(tnu1)
  assert is_tnu(tnu2)
  if tnu1 == tnu2:
    # If in differently checklists, this could be an incompatibility
    return rel.eq
  (peer1, peer2) = find_peers(tnu1, tnu2)  
  if peer1 != peer2:
    return rel.disjoint
  if peer1 == tnu1:
    return rel.gt
  if peer2 == tnu2:
    return rel.lt
  assert False

def are_disjoint(tnu1, tnu2):
  assert is_tnu(tnu1)
  assert is_tnu(tnu2)
  if tnu1 == forest_tnu: return False
  if tnu2 == forest_tnu: return False
  if tnu1 == tnu2: return False
  (tnu1, tnu2) = find_peers(tnu1, tnu2)
  return tnu1 != tnu2

# Find ancestor(s) of tnu1 and/or tnu2 that are in the same mutex: either
# disjoint or equal.

def find_peers(tnu_1, tnu_2):
  assert is_tnu(tnu_1)
  assert is_tnu(tnu_2)

  tnu1 = to_accepted(tnu_1)
  tnu2 = to_accepted(tnu_2)

  if tnu1 == forest_tnu or tnu2 == forest_tnu:
    return (forest_tnu, forest_tnu)
  assert get_checklist(tnu1) == get_checklist(tnu2)  #?

  mutex1 = get_mutex(tnu1)
  mutex2 = get_mutex(tnu2)

  if mutex1 == mutex2:
    if debug:
      print ("# Kludge %s %s" % (get_unique(tnu1), get_unique(tnu2)))
    tnu1 = get_parent(tnu1)
    mutex1 = get_mutex(tnu1)    

  #print("# Mutexes are %s %s" % (mutex1, mutex2))
  # Mutex of the forest is 0.  Going in 0-ward direction.

  while mutex1 != mutex2:
    assert mutex1 >= 0
    assert mutex2 >= 0
    if mutex1 > mutex2:
      # If p2 is closer to the root, try going rootward from p1
      if tnu1 == forest_tnu: return (forest_tnu, forest_tnu)
      tnu1 = get_parent(tnu1)
      mutex1 = get_mutex(tnu1)
    else: # mutex1 < mutex2:
      # If p1 is closer to the root, try going rootward from p2
      if tnu2 == forest_tnu: return (forest_tnu, forest_tnu)
      tnu2 = get_parent(tnu2)
      mutex2 = get_mutex(tnu2)

  if debug:
   print("# find_peers(%s, %s) = %s, %s" % \
        (get_unique(tnu_1), get_unique(tnu_2), 
         get_unique(tnu1), get_unique(tnu2)))
  return (tnu1, tnu2)

# Common ancestor - utility
# Also computes number of matched tips
# None (not 0) is the identity for mrca

def mrca(tnu1, tnu2):
  while True:
    assert is_tnu(tnu1)
    assert is_tnu(tnu2)
    # to_accepted ??
    if tnu1 == tnu2: return tnu1
    (tnu1, tnu2) = find_peers(tnu1, tnu2)
    assert get_mutex(tnu1) == get_mutex(tnu2)

mutex_table = {}

def set_mutex(tnu, mutex):
  have = mutex_table.get(tnu, mutex)
  if have != mutex:
    verb = "Promoting" if have > mutex else "Demoting"
    print("# ** %s %s, %s -> %s" % \
          (verb, get_unique(tnu),
           rank.mutex_to_name(have),
           rank.mutex_to_name(mutex)))
  mutex_table[tnu] = mutex

def get_mutex(tnu):
  if not tnu:
    # Above root of tree = forest_tnu
    return rank.forest
  probe = mutex_table.get(tnu)
  if probe: return probe
  mutex = get_mutex_really(tnu)
  assert mutex >= 0
  mutex_table[tnu] = mutex    # Perhaps amended later
  return mutex

# Higher numbers are more tipward.
# Parent's level > child's level for all children.

def get_mutex_really(tnu):

  probe = get_accepted(tnu)
  if probe: return get_mutex(probe)

  # Findest most rootward mutex of all the children
  children_mutex = rank.atom     # identity for min
  for child in get_children(tnu):
    children_mutex = min(children_mutex, get_mutex(child))

  # Treat given rank, if any, as normative
  if get_parent(tnu) == forest_tnu:
    mutex = rank.root
  else:
    nominal = get_nominal_mutex(tnu)
    mutex = nominal or (children_mutex - 10)
  set_mutex(tnu, mutex)

  # Demote children that have higher rank
  correct_children_mutexes(tnu, mutex)

  return mutex

def correct_children_mutexes(parent, parent_mutex):
  for child in get_children(parent):
    child_mutex = get_mutex(child)
    if child_mutex <= parent_mutex:
      if child_mutex == parent_mutex:
        print("# ** Child %s has same rank as parent %s" % \
              (get_unique(child), get_unique(parent)))
      else:
        print("# ** Child %s is of higher rank than parent %s" %\
              (get_unique(child), get_unique(parent)))
      if is_container(child):
        new_mutex = parent_mutex + 1 # demote!
        set_mutex(child, new_mutex)
        correct_children_mutexes(child, new_mutex) # ?
      else:
        set_mutex(child, parent_mutex + 10)  # demote!

def get_nominal_mutex(tnu):
  nominal = get_nominal_rank(tnu) # name of rank
  return rank.name_to_mutex(nominal)

def is_container(tnu):
  name = get_name(tnu).lower()
  return "unclassified" in name or \
         "incertae sedis" in name or \
         "unallocated" in name or \
         "unassigned" in name
  

# Test

if __name__ == '__main__':
  checklist = read_checklist("../../ncbi/2020-01-01/primates/")
  print (len(get_all_tnus(checklist)))
  tnus = checklist.get_index(taxon_id)
  tnu = get_tnu_with_id(checklist, '43777')
  print ("Specimen tnu:", tnu)
  print (get_value(tnu, checklist_field))
  print ("Name:", get_name(tnu))
  print ("Superior:", get_superior(tnu))
  synos = get_synonyms(tnu)
  print ("Synonyms:", list(map(get_name, synos)))
  print ("Back atcha:", [get_accepted(syno) for syno in synos])
  print ("Roots:", get_roots(checklist))
