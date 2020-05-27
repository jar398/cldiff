"""
 Makes a subset of a checklist based on one subtree of a taxonmoy.

 python3 subset_dwc.py [--taxonomy tax_dwc] source_dwc out_dwc

"""

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

    tid_column = None
    if "taxonID" in head:
      tid_column = head.index("taxonID") 
    aid_column = None
    if "acceptedNameUsageID" in head:
      aid_column = head.index("acceptedNameUsageID")

    if tid_column == None:      # usually 0
      print("No taxonID column found")

    with open(outpath, "w") as outfile:
      (delimiter, quotechar, mode) = csv_parameters(outpath)
      writer = csv.writer(outfile, delimiter=delimiter, quotechar=quotechar, quoting=mode)
      writer.writerow(head)
      for row in reader:
        if tid_column != None:
          tid = row[tid_column]
          if tid in all: 
            writer.writerow(row)
        if aid_column != None:
          aid = row[aid_column]
          if aid in all:
            writer.writerow(row)

def closure(topo, root_id):
  all = {}
  empty = []
  def descend(id):
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
    for row in reader:
      child_id = row[tid_column]
      parent_id = row[pid_column]
      childs = children.get(parent_id, None)
      if childs:
        childs.append(child_id)
      else:
        children[parent_id] = [child_id]
  return children

def csv_parameters(path):
  if path.endswith(".csv"):
    print("CSV")
    return (",", '"', csv.QUOTE_MINIMAL)
  else:
    print("TSV")
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
