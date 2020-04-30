"""
 Converts an NCBI dump to a DwCA TNU (taxon) file.
 Fold the scientific name (-> canonicalName) and authority (if it
 extends the scientific name) into the taxon record.

 python3 ncbi_to_dwca.py from to
   'from' is a directory containing .dmp files from ncbi
   'to' is a directory that will contain files for DwCA contents

 E.g.
  python3 ncbi_to_dwca.py ~/ncbi/2015-01-01/dump ~/ncbi/2015-01-01/dwca
"""

import sys, os, csv, argparse

def main(indir, outpath):
  assert os.path.exists(indir)
  nodes = read_nodes(os.path.join(indir, "nodes.dmp"))
  names = read_names(os.path.join(indir, "names.dmp"))
  merged = read_merged(os.path.join(indir, "merged.dmp"))
  (names, scinames, authorities) = move_names_to_nodes(names, nodes)
  emit_dwc(nodes, names, scinames, authorities, merged, outpath)

def write_row(writer,
              taxonID, parentNameUsageID, taxonRank,
              acceptedNameUsageID, scientificName, canonicalName,
              taxonomicStatus, nomenclaturalStatus):
  writer.writerow([taxonID, parentNameUsageID, taxonRank,
                   acceptedNameUsageID, scientificName, canonicalName,
                   taxonomicStatus, nomenclaturalStatus])

def emit_dwc(nodes, names, scinames, authorities, merged, outpath):
  outdir = os.path.basename(outpath)
  if not os.path.isdir(outdir): os.mkdir(outdir)
  (delimiter, quotechar) = choose_csv_parameters(outpath)
  print ("Writing", outpath)
  with open(outpath, "w") as outfile:
    writer = csv.writer(outfile, delimiter=delimiter, quotechar=quotechar)
    write_row(writer,
              "taxonID", "parentNameUsageID", "taxonRank",
              "acceptedNameUsageID", "scientificName", "canonicalName",
              "taxonomicStatus", "nomenclaturalStatus")
    for (id, parent_id, rank) in nodes:
      write_row(writer,
                id, parent_id, rank,
                None, authorities.get(id, None), scinames.get(id, None),
                "accepted", None)
    for (id, text, kind, spin) in names:
      if "BOLD:" in text or "bold:" in text:
        z = text.split("BOLD:")
        if len(z) == 2:
          write_row(writer,
                    id + ".BOLD", None, None,
                    id, None, "BOLD:" + z[1],
                    "synonym", "apparent BOLD id")
        else:
          print("Malformed BOLD: name, %s" % z)
      if kind != "scientific name":
        # synonym is a taxonomic status, not a nomenclatural status
        if kind == "synonym": kind = None
        minted = id + "." + str(spin)
        write_row(writer,
                  minted, None, None,
                  id, None, text,
                  "synonym", kind)
    for (old_id, new_id) in merged:
      canonical = "%s merged into %s" % (old_id, new_id)
      write_row(writer,
                old_id, None, None,
                new_id, None, canonical,
                "synonym", "merged id")


def move_names_to_nodes(names, nodes):
  keep = []
  scinames = {}
  for row in names:
    (id, text, kind, spin) = row
    if kind == "scientific name":
      scinames[id] = text
    else:
      keep.append(row)
  names = keep
  print (len(scinames), "canonicalNames (NCBI scientific names)")
  keep = []
  authorities = {}
  for row in names:
    (id, text, kind, spin) = row
    if kind == "authority":
      probe = scinames.get(id, None)
      if probe and text.startswith(probe):
        authorities[id] = text
      else:
        keep.append(row)
    else:
      keep.append(row)
  names = keep
  print (len(authorities), "authorities")
  return (names, scinames, authorities)

def read_nodes(nodes_path):
  nodes = []
  # Read the nodes file
  with open(nodes_path, "r") as infile:
    for row in csv.reader(infile, delimiter="\t", quotechar=quotechar_for_tab):
      # tax_id, |, parent tax_id, |, rank, ... other stuff we don't use ...
      nodes.append((row[0], row[2], row[4]))
  print (len(nodes), "nodes")
  return nodes

def read_names(names_path):
  names = []
  # Read the names file
  with open(names_path, "r") as infile:
    # Depends on names being grouped by taxa
    previous_id = None
    spin = -1
    for row in csv.reader(infile, delimiter='\t', quotechar=quotechar_for_tab):
      # tax_id, |, text, |, <unused>, |, name class, ... other stuff we don't use ...
      id = row[0]
      if str(id) != previous_id:
        spin = 1
        previous_id = id
      else:
        spin += 1
      names.append((id, row[2], row[6], spin))
  print (len(names), "names")
  return names

def read_merged(merged_path):
  merged = []
  # Read the merged file
  with open(merged_path, "r") as infile:
    for row in csv.reader(infile, delimiter="\t", quotechar=quotechar_for_tab):
      # old_tax_id, |, new_tax_id
      merged.append((row[0], row[2]))
  print (len(merged), "merged")
  return merged

def choose_csv_parameters(outpath):
  _, extension = os.path.splitext(outpath)
  return (',', '"') if extension == "csv" else ('\t', quotechar_for_tab)
quotechar_for_tab = '\a'

# When invoked from command line:

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('dump', help='directory containing taxdump files')
  parser.add_argument('--out', help='where to store the DwC version')
  args = parser.parse_args()
  main(args.dump, args.out)
