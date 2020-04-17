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
  (delimiter, quotechar) = choose_csv_parameters(checklist)
  with open(checklist, "r") as infile:
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar)
    head = next(reader)
    tid_column = head.index("taxonID") 
    aid_column = None
    if "acceptedNameUsageID" in head:
      aid_column = head.index("acceptedNameUsageID")

    with open(outpath, "w") as outfile:
      (delimiter, quotechar) = choose_csv_parameters(outpath)
      writer = csv.writer(outfile, delimiter=delimiter, quotechar=quotechar, quoting=csv.QUOTE_NONE)
      writer.writerow(head)
      for row in reader:
        tid = row[tid_column]
        if tid in all: 
          writer.writerow(row)
        elif aid_column:
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
  (delimiter, quotechar) = choose_csv_parameters(tax_path)
  with open(tax_path, "r") as infile:
    reader = csv.reader(infile, delimiter=delimiter, quotechar=quotechar)
    head = next(reader)
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

def choose_csv_parameters(outpath):
  _, extension = os.path.splitext(outpath)
  return (',', '"') if extension == "csv" else ('\t', '\a')

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
