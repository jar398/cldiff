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
  print("Roots:", roots)
  return roots

# Superior/inferior

def get_superior(tnu):
  assert tnu > 0
  return get_accepted(tnu) or get_parent(tnu)

def get_inferiors(tnu):
  assert tnu > 0
  return get_synonyms(tnu) + get_children(tnu)

# ----------
# Parent/children and accepted/synonyms

parents_cache = {}

def get_parent(tnu):
  probe = parents_cache.get(tnu)
  if probe != None: return probe
  parent = get_parent_really(tnu)
  parents_cache[tnu] = parent or False
  return parent

def get_parent_really(tnu):
  parent = None
  probe = get_accepted(tnu)
  if probe:
    # This happens a lot in GBIF.  Don't bother.
    if False:
      print("** In get_parent: %s has accepted name %s," %
            (get_unique(tnu), get_unique(probe)))
    tnu = probe
  parent_id = get_value(tnu, parent_tnu_id_field)
  if parent_id != None:
    # No parent means root taxon
    probe = get_tnu_with_id(get_checklist(tnu), parent_id)
    if probe:
      if is_accepted(tnu) and not is_accepted(probe):
        # GBIF contains all sorts of junk
        print("** Parent %s of accepted %s is not accepted" %
              (get_unique(probe), get_unique(tnu)))
        probe2 = get_accepted(probe)
        if probe2:
          parent = probe2
        elif not probe2:
          print("** Unaccepted parent %s of %s has no accepted taxon" %
                (get_unique(probe), get_unique(tnu)))
      else:
        parent = probe
    else:
      # This is a root
      pass
  return parent

def get_children(parent):
  return get_tnus_with_value(get_checklist(parent),
                             parent_tnu_id_field,
                             get_tnu_id(parent))

# ----------
# Accepted/synonyms

accepteds_cache = {}

def get_accepted(tnu):
  probe = accepteds_cache.get(tnu)
  if probe != None: return probe

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

  accepteds_cache[tnu] = accepted or False
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

def to_accepted(tnu):
  probe = get_accepted(tnu)
  if probe:
    return probe
  else:
    return tnu

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

def find_peers(tnu1, tnu2):
  if (not tnu1) or (not tnu2):
    return (None, None)
  assert tnu1 > 0
  assert tnu2 > 0
  assert get_checklist(tnu1) == get_checklist(tnu2)
  d1 = get_rank(tnu1)
  d2 = get_rank(tnu2)
  while True:
    parent1 = get_superior(tnu1)
    p1 = get_rank(parent1)
    # No leapfrogging
    if p1 < d2: break
    tnu1 = parent1
    d1 = p1
  while True:
    parent2 = get_superior(tnu2)
    p2 = get_rank(parent2)
    if p2 < d1: break
    tnu2 = parent2
    d2 = p2
  return (tnu1, tnu2)

def how_related(tnu1, tnu2):
  assert tnu1 > 0
  assert tnu2 > 0
  if tnu1 == tnu2:
    return rel.eq
  (peer1, peer2) = find_peers(tnu1, tnu2)  
  if peer1 != peer2:
    return rel.disjoint
  if peer1 == tnu1:
    # TBD: check for synonymy
    return rel.gt
  if peer2 == tnu2:
    return rel.lt
  assert False

def are_disjoint(tnu1, tnu2):
  assert tnu1 > 0
  assert tnu2 > 0
  tnu1 = to_accepted(tnu1)
  tnu2 = to_accepted(tnu2)
  (tnu1, tnu2) = find_peers(tnu1, tnu2)
  return tnu1 != tnu2

# Common ancestor - utility
# Also computes number of matched tips

def mrca(tnu1, tnu2):
  if not tnu1: return tnu2
  if not tnu2: return tnu1
  if tnu1 == tnu2: return tnu1
  d1 = get_rank(tnu1)
  d2 = get_rank(tnu2)
  while tnu1 != tnu2:
    if d1 >= d2:
      tnu1 = get_superior(tnu1)
      d1 = get_rank(tnu1)
    else:
      tnu2 = get_superior(tnu2)
      d2 = get_rank(tnu2)
  return tnu1

depth_cache = {}

def get_rank(tnu):
  if not tnu: return 0    # None is OK, means we fell off top
  assert tnu > 0
  depth = depth_cache.get(tnu, None)
  if depth: return depth
  d = 1
  rankname = get_nominal_rank(tnu)
  if rankname:
    r = rank.name_to_rank.get(rankname)
    if r:
      d = max(d, r)
  sup = get_superior(tnu)
  if sup:
    d = max(d, get_rank(sup) + 1)
  depth_cache[tnu] = d
  return d

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
