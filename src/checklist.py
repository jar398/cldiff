import os, csv

import relation as rel
import rank

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

registry = {0: "unused"}
sequence_numbers = {}

# Get the value of a field of a TNU record (via global registry)

def _get_record(uid):
  assert uid > 0
  (record, _) = registry[uid]
  return record

def get_value(uid, field):
  record = _get_record(uid)
  return record[field_position(field)]

def get_checklist(uid):
  assert uid > 0
  (_, checklist) = registry[uid]
  return checklist

def get_sequence_number(uid):
  assert uid > 0
  return sequence_numbers[uid]

# ---------- Checklists

class Checklist:
  def __init__(self):
    self.tnus = []
    self.indexes = [None] * len(the_fields)
    self.metadata = None
  def get_all_tnus(self):
    return self.tnus
  def tnu_count(self):
    return len(self.tnus)
  def add_tnu(self, record):
    uid = len(registry)
    registry[uid] = (record, self)
    self.tnus.append(uid)
    return uid
  def new_tnu(self):
    record = [None] * len(the_fields)
    return self.add_tnu(record)
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
    with open(inpath, "r") as infile:
      reader = csv.reader(infile, delimiter="\t", quotechar='\a', quoting=csv.QUOTE_NONE)
      head = next(reader)
      guide = {key: position for (key, position) in zip(head, range(len(head)))}
      for row in reader:
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
  checklist = Checklist()
  checklist.prefix = prefix
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
  return get_value(tnu, tnu_id_field)

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
  return get_value(tnu, taxon_rank_field)

# Unique name of the sort Nico likes

def get_unique(tnu):
  if not tnu: return "none"
  checklist = get_checklist(tnu)
  name = get_name(tnu)
  tnus_with_this_name = \
    get_tnus_with_value(checklist, canonical_name_field, name)

  if not is_accepted(tnu):
    if is_synonym(tnu):
      name = "?" + name
    else:
      status = get_value(tnu, nomenclatural_status_field)
      if not status: status = "nom. status not given"
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
  assert tnu > 0
  return get_synonyms(tnu) + get_children(tnu)

# ----------
# Parent/children and accepted/synonyms

superiors_cache = {}

def get_superior(tnu):
  (parent, accepted) = get_superiors(tnu)
  return parent or accepted

def get_parent(tnu):
  (parent, _) = get_superiors(tnu)
  return parent

def get_accepted(tnu):
  (_, accepted) = get_superiors(tnu)
  return accepted

def get_superiors(tnu):
  probe = superiors_cache.get(tnu)
  if probe: return probe
  result = get_superiors_really(tnu)
  superiors_cache[tnu] = result
  return result

def get_superiors_really(tnu):
  parent_id = get_value(tnu, parent_tnu_id_field)
  if parent_id:
    parent = get_tnu_with_id(get_checklist(tnu), parent_id)
    # If id doesn't resolve, it's a root
  else:
    parent = None

  accepted_id = get_value(tnu, accepted_tnu_id_field)
  if accepted_id:
    accepted = get_tnu_with_id(get_checklist(tnu), accepted_id)
    if accepted:
      a_status = get_taxonomic_status(accepted)
      if a_status != "accepted":
        print("** Putative accepted has wrong taxonomic status %s\n  %s -> %s" % \
              (a_status, get_unique(tnu), get_unique(accepted)))
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
    print("** Taxon has both accepted and parent links\n  %s -> %s, %s" % \
          (get_unique(tnu), get_unique(accepted), get_unique(parent)))
    print("** Preferring parent to accepted")
    parent = None
  return (to_accepted(parent), to_accepted(accepted))

def get_children(parent):
  return get_tnus_with_value(get_checklist(parent),
                             parent_tnu_id_field,
                             get_tnu_id(parent))

def to_accepted(tnu):
  if not tnu: return tnu
  probe = get_accepted(tnu)
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

# Find ancestor(s) of tnu1 and tnu2 that are in the same mutex: either
# disjoint or equal.

def find_peers(tnu1, tnu2):
  if (not tnu1) or (not tnu2):
    return (None, None)
  assert tnu1 > 0
  assert tnu2 > 0
  assert get_checklist(tnu1) == get_checklist(tnu2)  #?
  tnu1 = to_accepted(tnu1)
  tnu2 = to_accepted(tnu2)
  mutex1 = get_mutex(tnu1)
  mutex2 = get_mutex(tnu2)
  while mutex1 != mutex2:
    if mutex1 < mutex2:
      # If p1 is closer to the root, try going rootward from p2
      tnu2 = get_superior(tnu2)
      mutex2 = get_mutex(tnu2)
    else:
      # If p2 is closer to the root, try going rootward from p1
      tnu1 = get_superior(tnu1)
      mutex1 = get_mutex(tnu1)
  return (tnu1, tnu2)

def how_related(tnu1, tnu2):
  assert tnu1 > 0
  assert tnu2 > 0
  tnu1 = get_accepted(tnu1)    # TBD: implement in-part 'synonyms'
  tnu2 = get_accepted(tnu2)
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
  assert tnu1 > 0
  assert tnu2 > 0
  (tnu1, tnu2) = find_peers(tnu1, tnu2)
  return tnu1 != tnu2

# Common ancestor - utility
# Also computes number of matched tips

def mrca(tnu1, tnu2):
  if not tnu1 or not tnu2: return tnu1
  tnu1 = to_accepted(tnu1)
  tnu2 = to_accepted(tnu2)
  while tnu1 != tnu2:
    tnu1 = get_superior(tnu1)
    tnu2 = get_superior(tnu2)
    (tnu1, tnu2) = find_peers(tnu1, tnu2)
  return tnu1

def get_mutex(tnu):
  return rank.height_to_depth(level_table.get(tnu, rank.root_height))

level_table = {}

def get_height(tnu):
  if not tnu:
    return anything_depth
  probe = level_table.get(tnu)
  if probe: return probe
  height = get_height_really(tnu)
  assert height >= 0
  level_table[tnu] = height    # Perhaps amended later
  return height

def get_height_really(tnu):
  # Must be higher than any child
  level = 0    # Can increase
  inferiors = get_inferiors(tnu)
  if len(inferiors) > 0:
    # The children are mutually exclusive
    unplaced_level = 0
    for inf in inferiors:
      h = get_height(inf)
      level = max(level, h)
      if is_container(h):
        unplaced_level = max(unplaced_level, h)
    if unplaced_level >= level:
      print("** Shouldn't happen: height(%s) >= height(%s)" % \
            (get_unique(inf), get_unique(tnu)))
      level += 1
    for inf in get_inferiors(tnu):
      if not is_container(inf):
        # This may change the value that's already there
        h = level_table[inf]
        if h != level:
          print("** Promotion: %s from %s to %s" % \
                (get_unique(inf), h, level))
          level_table[inf] = level
    if not is_container(tnu):
      level += 1
      rank = get_nominal_rank(tnu)
      if rank:
        r = rank.rank_to_height(rank)
        if r > level:
          print("** Promoting %s from %s to %s" % \
                (get_unique(tnu), level, r))
          level = r
  return level

def is_container(tnu):
  name = get_name(inf).lower()
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
