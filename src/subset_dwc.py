"""
 Makes a subset of a checklist based on one subtree of a taxonmoy.

 python3 subset_dwc.py [--taxonomy tax_dwc] source_dwc out_dwc

 Assumption: every accepted record has a taxonID
"""

debug = False

import sys, os, csv, argparse

def main(checklist, tax_path, root_id, outpath):
  topo = read_topology(tax_path)
  print ("Nodes with children:", len(topo))
  all = closure(topo, root_id)
  print ("Nodes in subset:", len(all))
  write_subset(checklist, root_id, all, outpath)

def write_subset(checklist, root_id, all, outpath):

  (delimiter, quotechar, mode) = csv_parameters(checklist)
  with open(checklist, "r") as infile:
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
    head = next(reader)

    tid_column = head.index("taxonID") 
    aid_column = head.index("acceptedNameUsageID")
    pid_column = head.index("parentNameUsageID")
    sid_column = head.index("taxonomicStatus")

    with open(outpath, "w") as outfile:
      (delimiter, quotechar, mode) = csv_parameters(outpath)
      writer = csv.writer(outfile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
      writer.writerow(head)
      for row in reader:
        if accepted(row, tid_column, sid_column, aid_column):
          # Accepted record
          if row[tid_column] in all:
            row[aid_column] = ''  # Clobber distracting accepted id
            writer.writerow(row)
        else:
          # Synonym record?
          if row[aid_column] in all:
            row[pid_column] = ''    # Clobber distracting parent id
            if sid_column and row[sid_column] == "accepted":
              print("** Corrupt taxonomic status for %s" % row[tid_column])
            writer.writerow(row)

# Transitive closure of accepted records

def closure(topo, root_id):
  all = {}
  empty = []
  def descend(id):
    if id in all:
      print("** Node has multiple parents: %s" % (id))
    else:
      all[id] = True
      for child in topo.get(id, empty):
        descend(child)
  descend(root_id)
  return all

def read_topology(tax_path):
  children = {}
  (delimiter, quotechar, mode) = csv_parameters(tax_path)
  with open(tax_path, "r") as infile:
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
    head = next(reader)
    print("Header row:", head)
    tid_column = head.index("taxonID") 
    pid_column = head.index("parentNameUsageID")
    aid_column = head.index("acceptedNameUsageID")
    sid_column = head.index("taxonomicStatus")

    if tid_column == None:      # usually 0
      print("** No taxonID column found")
      return -1
    if aid_column == None:      # usually 0
      print("** No acceptedNameUsageID column found")
      return -1

    for row in reader:
      if accepted(row, tid_column, sid_column, aid_column):
        tid = row[tid_column]
        parent_id = row[pid_column]
        if parent_id != '':
          childs = children.get(parent_id, None)
          if childs:
            childs.append(tid)
          else:
            children[parent_id] = [tid]
  return children

def accepted(row, tid_column, sid_column, aid_column):
  tid = row[tid_column]
  if tid == '':
    return False
  if sid_column:
    if row[sid_column] == 'synonym':
      return False
  if aid_column:
    aid = row[aid_column]
    if aid == '' or aid == tid:
      return True
  return False

def csv_parameters(path):
  if path.endswith(".csv"):
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    return ("\t", "\a", csv.QUOTE_NONE)

# main(checklist, taxonomy, root_id, outfile)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--taxonomy', help="""taxonomy from which to extract
                            hierarchy; defaults to source""")
  parser.add_argument('source', help='taxonomy or checklist from which to extract subset')
  parser.add_argument('id', help="taxon id of subset's root")
  parser.add_argument('--out', help='where to store the subset')
  args = parser.parse_args()
  main(args.source, args.taxonomy or args.source, args.id, args.out)
