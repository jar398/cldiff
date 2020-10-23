import csv
import property

# A table can be read and/or written
# If a table is populated it can be indexed
# The columns of a table have labels

class Table:
  def __init__(self):
    self.record_uids = []
    pass

  def header(self):
    return self.header

  # Position in record of column for given property
  def get_position(self, prop):
    return self.position_index[prop.uid]

  def process_header(self, header):
    assert not "\t" in header[0]
    assert not "," in header[0]
    self.header = header
    self.position_index = [None] * property.number_of_properties
    self.methods = [not_present] * property.number_of_properties
    # TBD: If there is a meta.xml, get the properties that way.
    # NB: by_name returns None if label is unrecognized
    self.properties = [property.by_name(label) for label in header]
    self.indexes = [None] * property.number_of_properties
    for position in range(len(header)):
      # position is column position within record (table specific)
      prop = self.properties[position]
      if prop:
        self.position_index[prop.uid] = position
        self.methods[prop.uid] = self.fetcher(prop, position)

  def fetcher(self, prop, position):
    def fetch(r):
      if r[position] == '': return None
      return r[position]
    return fetch

  def populate_from_generator(self, record_generator):
    self.process_header(next(record_generator))
    for record in record_generator:
      id = _register(record, self)
      self.record_uids.append(id)

  def populate_from_file(self, inpath):
    # Look for a meta.xml file in same directory?
    (delim, qc, qu) = csv_parameters(inpath)
    # print("# Parameters %s %s %s" % (delim, qc, qu))
    with open(inpath, "r") as infile:
      reader = csv.reader(infile, delimiter=delim, quotechar=qc, quoting=qu)
      self.populate_from_generator(reader)

  # Create indexes on demand.  Position is column position specific to
  # this table, which can be determined using get_position.

  def get_index(self, prop):
    if self.methods[prop.uid] == not_present: return {}
    if self.indexes[prop.uid] == None:
      index = {}
      for id in self.record_uids:
        value = get_value(id, prop)
        if value != None:
          if value in index:
            index[value].append(id)
          else:
            index[value] = [id]
      self.indexes[prop.uid] = index
    return self.indexes[prop.uid]

def csv_parameters(path):
  if ".csv" in path:
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

def read_table(specifier):
  table = Table()
  table.populate_from_file(specifier)
  return table

# Record registry

def is_record(x):
  return isinstance(x, int) and x > 0

_registry = ["there is no record 0"]

# We store a pair (record, table) where table is the table that "owns"
# the record.

def record_and_table(record_uid):    # returns (record, table)
  return _registry[record_uid]

def _register(record, table):
  record_uid = len(_registry)
  _registry.append((record, table))
  return record_uid

def not_present(record):
  return None

def get_value(record_uid, prop):
  (r, t) = record_and_table(record_uid)
  val = t.methods[prop.uid](r)
  return val

def get_table(record_uid):
  (r, t) = record_and_table(record_uid)
  return t

# ---------- Self-test

def self_test():
  table = read_table("work/ncbi/2020-01-01/primates.csv")
  print ("Records read: %s" % len(table.record_uids))
  prop = property.by_name("taxonID")
  pos = table.get_position(prop)    # column number
  print ("taxonID position in table is %s" % pos)
  assert pos != None
  print ("taxonID position in properties is %s" % prop.uid)

  rec = table.record_uids[0]
  (r, _) = record_and_table(rec)
  print ("Sample record: %s" % r)
  print ("Taxon id of sample record: %s" % get_value(rec, prop))

  idx = table.get_index(prop)
  print ("Ids indexed: %s" % len(idx))
  print (idx["9443"])

if __name__ == '__main__':
  self_test()
