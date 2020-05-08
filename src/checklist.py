debug = False

import os, csv

import relation as rel
import rank
import chaitin

# ---------- Fields

# Darwin Core field = (label, position)

the_fields = []
label_to_field = {}
def define_field(label):
  field = (label, len(the_fields))
  the_fields.append(field)
  label_to_field[label] = field
  return field

def field_label(field):
  (label, position) = field
  return label

def field_position(field):
  (label, position) = field
  return position

# Particular fields of interest here (not all possible DwC fields)

tnu_id_field               = define_field("taxonID")
accepted_tnu_id_field      = define_field("acceptedNameUsageID")
taxonomic_status_field     = define_field("taxonomicStatus")
nomenclatural_status_field = define_field("nomenclaturalStatus")
canonical_name_field       = define_field("canonicalName")    # without authority
scientific_name_field      = define_field("scientificName")  # with authority
parent_tnu_id_field        = define_field("parentNameUsageID")
taxon_rank_field           = define_field("taxonRank")

# ---------- Taxon registry and taxa

# Registry = tnu uid -> (value vector, checklist)
#  where value vector is a vector that parallels the `the_fields` list

def is_tnu(x):
  return isinstance(x, int) and x >= 0

forest_tnu = 0
registry = {forest_tnu: "no forest record"}
sequence_numbers = {}

# Get the value of a field of a TNU record (via global registry)

def _get_record(uid):
  assert is_tnu(uid)
  (record, _) = registry[uid]
  return record

def get_value(uid, field):
  if uid == forest_tnu: return None
  record = _get_record(uid)
  return record[field_position(field)]

def get_checklist(uid):
  assert is_tnu(uid)               # Forest has no checklist
  if uid == forest_tnu: return None
  (_, checklist) = registry[uid]
  return checklist

def get_sequence_number(uid):
  assert uid                    # Forest has no checklist
  return sequence_numbers[uid]

def fresh_tnu(checklist):
  record = [None] * len(the_fields)
  uid = len(registry)
  registry[uid] = (record, checklist)
  return uid

# ---------- Checklists

class Checklist:
  def __init__(self, prefix):
    assert prefix
    self.prefix = prefix
    self.tnus = []
    self.indexes = [None] * len(the_fields)
    self.metadata = None
  def get_all_tnus(self):
    return self.tnus
  def tnu_count(self):
    return len(self.tnus)
  def new_tnu(self):
    tnu = fresh_tnu(self)
    self.tnus.append(tnu)
    return tnu
  def get_index(self, field):
    (_, position) = field
    # create the index if necessary, on demand
    if self.indexes[position] == None:
      index = {}
      for tnu in self.tnus:
        # Check each record
        (record, _) = registry[tnu]
        value = record[position]
        if value:
          if value in index:
            index[value].append(tnu)
          else:
            index[value] = [tnu]
      self.indexes[position] = index
    return self.indexes[position]
  def populate_from_directory(self, indirs):
    self.populate_from_file(get_tnu_path(indir))
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
    guide = {key: position for (key, position) in zip(head, range(len(head)))}
    if guide.get("taxonID", None) == None:
      print("** Didn't find a column for taxonID")
      print("%s" % guide)
      assert False
    for row in row_generator:
      uid = self.new_tnu()
      record = _get_record(uid)
      for (label, (_, position)) in label_to_field.items():
        position_in_row = guide.get(label, None)
        # Idiot python treats 0 as false
        if position_in_row != None:
          value = row[position_in_row]
          if value != '':
            record[position] = value
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

def read_checklist(specifier, prefix):
  assert prefix
  checklist = Checklist(prefix)
  checklist.populate_from_file(specifier)
  return checklist

# Utility - copied from another file - really ought to be shared

def get_tnu_path(dwca_dir):
  for name in ["taxon.tsv",
               "Taxon.tsv",
               "taxa.txt",
               "taxon.txt",
               "Taxon.txt"]:
    path = os.path.join(dwca_dir, name)
    if os.path.exists(path):
      return path
  raise ValueError("cannot find TNU file in this directory", dwca_dir)

# -------------------- indexing

# Get TNUs in checklist having a given value

canonical_empty_list = []

def get_tnus_with_value(checklist, field, value):
  return checklist.get_index(field).get(value, canonical_empty_list)

# Get unique (we hope) TNU possessing a given identifier

def get_tnu_id(tnu):
  id = get_value(tnu, tnu_id_field)
  assert id != None
  return id

def get_tnu_with_id(checklist, id):
  tnus = checklist.get_index(tnu_id_field).get(id, None)
  if tnus:
    return tnus[0]
  else:
    return None

  return get_value(tnu, tnu_id_field)

# Create a dict mapping field values to lists of keys (uids, tnus)
#  - deprecated

def index_by_column(checklist, field):
  return checklist.get_index(field)

# You could use this for, say, taxonID, or maybe authority (scientificName)
# Sort of a kludge, we really shouldn't have to copy the dict
#  - deprecated

def index_unique_by_column(checklist, field):
  index = checklist.get_index(field)
  unique = {}
  for key, tnus in index.items():
    if len(tnus) != 1:
      raise ValueError("conflict in supposedly unique index: %s -> %s, %s" % (key, tnus, field))
    unique[key] = tnus[0]
  return unique

# ----------------------------------------

# Logic for particular fields

# The following is for display purposes

def get_name(tnu):
  name = get_value(tnu, canonical_name_field)
  if name != None: return name
  name = get_value(tnu, scientific_name_field)
  if name != None: return name  
  return get_tnu_id(tnu)

def get_nominal_rank(tnu):
  if is_container(tnu): return None
  return get_value(tnu, taxon_rank_field)

# Unique name of the sort Nico likes

def get_unique(tnu):
  if tnu == None: return "none"
  assert is_tnu(tnu)
  if tnu == forest: return "forest"
  checklist = get_checklist(tnu)
  name = get_name(tnu)
  tnus_with_this_name = \
    get_tnus_with_value(checklist, canonical_name_field, name)

  status = get_value(tnu, taxonomic_status_field)
  if status == "accepted" or not status:
    pass
  elif status == "synonym":
    name = "?" + name
  else:
    name = "%s (%s)" % (name, status)

  better = checklist.prefix + name.replace(" ", "_")
  if len(tnus_with_this_name) <= 1:
    return better
  else:
    return checklist.prefix + name + "#" + get_tnu_id(tnu)

# Roots - accepted tnus without parents

def get_roots(checklist):
  roots = []
  for tnu in get_all_tnus(checklist):
    if not get_accepted(tnu) and not get_parent(tnu):
      roots.append(tnu)
  return roots

# Superior/inferior

def get_inferiors(tnu):
  assert is_tnu(tnu)
  assert tnu != forest
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
  assert parent != None
  return parent

def get_accepted(tnu):
  assert is_tnu(tnu)
  (accepted, _) = get_superiors(tnu)
  return accepted

# Returns (accepted, parent)

def get_superiors(tnu):
  assert is_tnu(tnu)
  probe = superiors_cache.get(tnu)
  if probe: return probe
  result = get_superiors_really(tnu)
  superiors_cache[tnu] = result
  if False:
    print("# Cached gs(%s) = (%s, %s)" % \
          (get_unique(tnu), get_unique(result[0]), get_unique(result[1])))
  return result

def get_superiors_really(tnu):
  assert is_tnu(tnu)
  if tnu == forest_tnu: return (None, None)
  parent_id = get_value(tnu, parent_tnu_id_field)
  if parent_id == None:
    parent = forest_tnu
  else:
    parent_uid = get_tnu_with_id(get_checklist(tnu), parent_id)
    if parent_uid == None:
      parent = forest_tnu
    else:
      parent = to_accepted(parent_uid)
  
  # If id doesn't resolve, it's a root

  accepted_id = get_value(tnu, accepted_tnu_id_field)
  if accepted_id:
    accepted = to_accepted(get_tnu_with_id(get_checklist(tnu), accepted_id))
    if accepted:
      accepted_taxo_status = get_taxonomic_status(accepted)
      if accepted_taxo_status and accepted_taxo_status != "accepted":
        # Example: GBIF backbone:
        # C.?Hylobates sericus#8955169 -> C.Bunopithecus_sericus_(doubtful)
        print("** Putative accepted has wrong taxonomic status %s\n  %s -> %s" % \
              (accepted_taxo_status, get_unique(tnu), get_unique(accepted)))
      # Iterate???
    else:
      print("** Id %s for accepted of %s doesn't resolve" % \
            (accepted_id, get_unique(tnu)))
  else:
    n_status = get_taxonomic_status(tnu)
    if n_status == "synonym":
      print("** Node with status = 'synonym' has no accepted node\n  %s" % \
            (get_unique(tnu)))

    accepted = None

  if parent and accepted:
    parent_taxo_status = get_taxonomic_status(parent)
    if parent_taxo_status == "accepted" and accepted_taxo_status == "accepted":
      if False:
        print("** Taxon has both accepted and parent links\n  %s -> %s, %s" % \
              (get_unique(tnu), get_unique(accepted), get_unique(parent)))
        print("** Preferring parent to accepted")
      parent = None
    else:
      # This case is common in GBIF files
      accepted = None
  if parent != None: parent = to_accepted(parent)
  if accepted != None: accepted = to_accepted(accepted)
  return (accepted, parent)

def get_children(parent):
  return get_tnus_with_value(get_checklist(parent),
                             parent_tnu_id_field,
                             get_tnu_id(parent))

def to_accepted(tnu):
  assert is_tnu(tnu)
  probe = get_accepted(tnu)     # cached
  if probe == None:
    return tnu
  else:
    return probe

# ----------
# Accepted/synonyms

accepteds_cache = {}

def get_accepted_really(tnu):
  accepted = None
  if is_accepted(tnu):
    pass                        # This is normal
  else:
    accepted_id = get_value(tnu, accepted_tnu_id_field)
    if not accepted_id:
      pass                      # This is normal
    else:
      probe = get_tnu_with_id(get_checklist(tnu), accepted_id)
      if probe:
        if is_accepted(probe):
          accepted = probe      # Normal case
        else:
          print("** Accepted taxon %s for %s is not accepted" %
                (get_unique(probe), get_unique(tnu)))
      else:
        print("** Accepted id %s for %s doesn't resolve" %
              (accepted_id, get_unique(tnu)))

  return accepted

def get_synonyms(tnu):
  return [syn
          for syn in get_tnus_with_value(get_checklist(tnu),
                                         accepted_tnu_id_field,
                                         get_tnu_id(tnu))
          if not get_parent(tnu)]

def get_taxonomic_status(tnu):
  return get_value(tnu, taxonomic_status_field)

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
  assert is_tnu(tnu)
  if tnu1 == forest: return False
  if tnu2 == forest: return False
  (tnu1, tnu2) = find_peers(tnu1, tnu2)
  return tnu1 != tnu2

# Find ancestor(s) of tnu1 and tnu2 that are in the same mutex: either
# disjoint or equal.

def find_peers(tnu1, tnu2):
  assert is_tnu(tnu1)
  assert is_tnu(tnu2)

  tnu1 = to_accepted(tnu1)
  tnu2 = to_accepted(tnu2)

  if tnu1 == forest_tnu or tnu2 == forest_tnu:
    return (forest_tnu, forest_tnu)
  assert get_checklist(tnu1) == get_checklist(tnu2)  #?

  mutex1 = get_mutex(tnu1)
  mutex2 = get_mutex(tnu2)

  #print("# Mutexes are %s %s" % (mutex1, mutex1))

  while mutex1 != mutex2:
    assert mutex1 >= 0
    assert mutex2 >= 0
    if mutex1 > mutex2:
      # If p1 is closer to the root, try going rootward from p2
      tnu2 = get_parent(tnu2)
      if tnu2 == forest_tnu: return (forest_tnu, forest_tnu)
      mutex2 = get_mutex(tnu2)
    else: # mutex1 < mutex2:
      # If p2 is closer to the root, try going rootward from p1
      tnu1 = get_parent(tnu1)
      if tnu1 == forest_tnu: return (forest_tnu, forest_tnu)
      mutex1 = get_mutex(tnu1)
    #print("# Adjusted mutexes are %s %s" % (mutex1, mutex1))
  return (tnu1, tnu2)

forest = 0

# Common ancestor - utility
# Also computes number of matched tips
# None (not 0) is the identity for mrca

def mrca(tnu1, tnu2):
  if tnu1 == None: return tnu2
  if tnu2 == None: return tnu1
  assert is_tnu(tnu1)
  assert is_tnu(tnu2)
  while True:
    (tnu1, tnu2) = find_peers(tnu1, tnu2)
    if tnu1 == tnu2: return tnu1
    tnu1 = get_parent(tnu1)
    tnu2 = get_parent(tnu2)

mutex_table = {}

def set_mutex(tnu, mutex):
  have = mutex_table.get(tnu, mutex)
  if have != mutex:
    verb = "Promoting" if have > mutex else "Demoting"
    print("** %s %s, %s -> %s" % \
          (verb, get_unique(tnu),
           rank.mutex_to_name(have),
           rank.mutex_to_name(mutex)))
  mutex_table[tnu] = mutex

def get_mutex(tnu):
  if not tnu:
    # Above root of tree = forest
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
        print("** Child %s has same rank as parent %s" % \
              (get_unique(child), get_unique(parent)))
      else:
        print("** Child %s is of higher rank than parent %s" %\
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
  tnus = checklist.get_index(tnu_id_field)
  tnu = get_tnu_with_id(checklist, '43777')
  print ("Specimen tnu:", tnu)
  print (registry[tnu][0])
  print ("Name:", get_name(tnu))
  print ("Superior:", get_superior(tnu))
  synos = get_synonyms(tnu)
  print ("Synonyms:", list(map(get_name, synos)))
  print ("Back atcha:", [get_accepted(syno) for syno in synos])
  print ("Roots:", get_roots(checklist))
