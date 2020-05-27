import csv

# A table can be read and/or written
# If a table is populated it can be indexed
# The columns of a table have labels

class Table:
  def __init__(self):
    self.records = []
    pass

  def header(self):
    return self.header

  # Position in record of column having given label
  def get_position(self, label):
    return self.position_index.get(label)

  def process_header(self, header):
    self.header = header
    self.position_index = {}
    for position in range(len(header)):
      label = header[position]
      self.position_index[label] = position
    print(self.position_index)
    self.indexes = [None] * len(header)

  def populate_from_generator(self, record_generator):
    self.process_header(next(record_generator))
    for record in record_generator:
      # rerepresent ?  normalize ?
      self.records.append(record)
      _register(record, self)
    # sort ??

  def populate_from_file(self, inpath):
    # Look for a meta.xml file in same directory?
    (delim, qc, qu) = csv_parameters(inpath)
    print ("delimiter: %s quote char: %s" % (delim, qc))
    with open(inpath, "r") as infile:
      reader = csv.reader(infile, delimiter=delim, quotechar=qc, quoting=qu)
      self.populate_from_generator(reader)

  # Create indexes on demand
  def get_index(self, position):
    if self.indexes[position] == None:
      index = {}
      for record in self.records:
        value = record[position]
        if value:
          if value in index:
            index[value].append(record)
          else:
            index[value] = [record]
      self.indexes[position] = index
    return self.indexes[position]

def csv_parameters(path):
  if path.endswith(".csv"):
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

def read_table(specifier):
  table = Table()
  table.populate_from_file(specifier)
  return table

# Record registry

def is_record_id(x):
  return isinstance(x, int) and x > 0

_registry = ["there is no record 0"]

# We store a pair (record, table) where table is the table that "owns"
# the record.

def record_and_table(uid):    # returns (record, table)
  return _registry[uid]

def _register(record, table):
  uid = len(_registry)
  _registry.append((record, table))
  return uid



if __name__ == '__main__':
  table = read_table("work/ncbi/2020-01-01/primates.csv")
  print (len(table.records))
  pos = table.get_position("taxonID")    # column number
  assert pos != None
  idx = table.get_index(pos)
  print (idx["9443"])
