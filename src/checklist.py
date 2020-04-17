import os, csv

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

# Registry = tnu uid -> (value vector, checklist)
#  where value vector is a vector that parallels the `the_fields` list

registry = {0: "unused"}

# Get the value of a field of a TNU record (via global registry)

def get_value(uid, field):
  (values, checklist) = registry[uid]
  (label, position) = field
  return values[position]

def get_checklist(uid):
  (_, checklist) = registry[uid]
  return checklist

# Checklist and registry

class Checklist:
  def __init__(self):
    self.tnus = []
    self.indexes = [None] * len(the_fields)
    self.metadata = None
  def add_tnu(self, values):
    uid = len(registry)
    registry[uid] = (values, self)
    self.tnus.append(uid)
  def get_all_tnus(self):
    return self.tnus
  def get_index(self, field):
    (_, position) = field
    # create the index if necessary, on demand
    if not self.indexes[position]:
      index = {}
      for tnu in self.tnus:
        # Check each row
        (values, _) = registry[tnu]
        value = values[position]
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
        values = [None] * len(the_fields)
        for (label, (_, position)) in label_to_field.items():
          position_in_row = guide.get(label, None)
          # Idiot python treats 0 as false
          if position_in_row != None:
            value = row[position_in_row]
            if value != '':
              values[position] = value
        self.add_tnu(values)

# A checklist's TNU list

def get_all_tnus(checklist):
  return checklist.get_all_tnus()

# Read a checklist from a file

def read_checklist(specifier, prefix):
  checklist = Checklist()
  checklist.populate_from_file(specifier)
  checklist.prefix = prefix
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

def get_tnus_with_value(checklist, field, value):
  return checklist.get_index(field).get(value, ())

# Get unique (we hope) TNU possessing a given identifier

def get_tnu_with_id(checklist, id):
  index = checklist.get_index(tnu_id_field)
  tnus = index.get(id, None)
  if tnus:
    return tnus[0]
  else:
    return None

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
  if name != None:
    return name
  return get_value(tnu, scientific_name_field)

# Unique name of the sort Nico likes

def get_unique(tnu):
  checklist = get_checklist(tnu)
  name = get_name(tnu)
  tnus_with_this_name = \
    get_tnus_with_value(checklist, canonical_name_field, name)

  better = checklist.prefix + name.replace(" ", "_")
  if len(tnus_with_this_name) <= 1:
    return better
  else:
    return checklist.prefix + name + "#" + get_value(tnu, tnu_id_field)

# All synonyms, sorted

def get_synonyms(tnu):
  id = get_value(tnu, tnu_id_field)
  if id is None: return ()
  index = get_checklist(tnu).get_index(accepted_tnu_id_field)
  return sorted(index.get(id, ()), key=badness)

def is_accepted(tnu):
  return get_value(tnu, taxonomic_status_field) == "accepted"

def get_accepted(A_candidate):
  if is_accepted(A_candidate):
    return A_candidate
  else:
    A_accepted_id = get_value(A_candidate, accepted_tnu_id_field)
    if A_accepted_id == None:
      print ("tnu is accepted, but it has no accepted id", A_candidate)
    return get_tnu_with_id(get_checklist(A_candidate), A_accepted_id)

# Roots - accepted tnus without parents

def get_roots(checklist):
  roots = []
  for tnu in get_all_tnus(checklist):
    if is_accepted(tnu):
      parent = get_parent(tnu)
      if parent is None:
        roots.append(tnu)
  print (len(roots), "roots")
  return roots

# Parent/children

def get_parent(tnu):
  parent_id = get_value(tnu, parent_tnu_id_field)
  if parent_id != None:
    return get_tnu_with_id(get_checklist(tnu), parent_id)
  else:
    return None

# List of child tnus, or () if none

def get_children(parent):
  parent_id = get_value(parent, tnu_id_field)
  return get_tnus_with_value(get_checklist(parent), parent_tnu_id_field, parent_id)

# For assigning priorities to synonyms

def badness(tnu):
  status = get_value(tnu, nomenclatural_status_field)
  if status is None:
    status = get_value(tnu, taxonomic_status_field)
    if status is None:
      return 99
  badness = badnesses.get(status, None)
  if badness is None: badness = 99
  return badness

# Name classes, best to worst

badnesses = {
  "authority": 0,
  "scientific name": 1,        # (actually canonical) exactly one per node
  "accepted": 1,
  "equivalent name": 2,        # synonym but not nomenclaturally
  "misspelling": 3,
  "genbank synonym": 4,        # at most one per node; first among equals
  "synonym": 5,
  "anamorph": 5.1,
  "teleomorph": 5.2,
  "misnomer": 5.5,
  "includes": 6,
  "in-part": 6.5,              # this node is part of a polyphyly
  "type material": 7,
  "blast name": 8,             # large well-known taxa
  "genbank common name": 9,    # at most one per node
  "genbank acronym": 9.2,      # at most one per node
  "genbank anamorph": 9.4,     # at most one per node
  "common name": 10,
  "acronym": 10.5,
  "unpublished name": 10.7,    # non-code synonym
  "id": 11,
  "merged id": 12,
}

# Totally general utilities from here down... I guess...

def invert_dict(d):
  inv = {}
  for (key, val) in d.items():
    if val in inv:
      inv[val].append(key)
    else:
      inv[val] = [key]
  return inv

# Test

if __name__ == '__main__':
  checklist = read_checklist("../../ncbi/2020-01-01/primates/")
  print (len(get_all_tnus(checklist)))
  tnus = checklist.get_index(tnu_id_field)
  tnu = get_tnu_with_id(checklist, '43777')
  print ("Specimen tnu:", tnu)
  print (registry[tnu][0])
  print ("Name:", get_name(tnu))
  print ("Parent:", get_parent(tnu))
  synos = get_synonyms(tnu)
  print ("Synonyms:", list(map(get_name, synos)))
  print ("Back atcha:", [get_accepted(syno) for syno in synos])
  print ("Roots:", get_roots(checklist))
