# Thanks to Open Tree

# Depth/mutex counts the number of splittings, starting at 0 = fewest
# bifurcations (everything is in one group; 'forest').

# Approximately 50 currently named ranks

rank_configuration = [
  [[],
   ["atom"]],  # added by jar for fun

  [["no rank",       # NCBI
    "population"],     
   ["cluster",      # added by jar, for SILVA ?
    "subform",
    "forma",
    "form"]],           # 2016 GBIF

  [["subvariety",
    "varietas"],
   ["variety",
    "natio"]],          # worms

  [["subspecies",
    "infraspecificname"],
   ["species",
    "species subgroup",
    "species group"]],

  [["subgenus"],
   ["genus",
    "subtribe",
    "tribe",
    "supertribe"]],     # worms

  [["subfamily"],
   ["family",
    "superfamily",
    "subsection",     # worms
    "section"]],        # worms

  [["parvorder",
    "infraorder",
    "suborder"],
   ["order",
    "superorder",
    "subcohort",      # NCBI
    "cohort"]],         # NCBI Polyneoptera

  [["subterclass",    # worms Colobognatha
    "infraclass",
    "subclass"],
   ["class",
    "superclass"]],

  [["subdivision",     # worms
    "infraphylum",     # worms
    "subphylum"],
   ["phylum",
    "superphylum"]],

  [["infrakingdom",     # worms
    "division",         # h2007
    "subkingdom"],
   ["kingdom",
    "superkingdom",
    "domain",
    "root",
    "checklist"]],

  [[], ["forest"]]    # added by JAR for fun
  ]

linnean_increment = 100000

# The least deep rank is 'forest'
forest = 0

# The deepest rank is 'atom'  (identity for 'min')
atom = (len(rank_configuration)+1) * linnean_increment

def process_ranks(groups):
  name_to_mutex_table = {}
  mutex_to_name_table = {}
  def set_mutex(name, mutex):
    name_to_mutex_table[name] = mutex
    mutex_to_name_table[mutex] = name
  foo = 10
  sublinnean_increment = int(linnean_increment / foo)
  mutex = atom
  for (colinnea, linnea) in groups:
    d = mutex
    colinnea.reverse()
    for name in colinnea:
      d += sublinnean_increment
      set_mutex(name, d)
    d = mutex
    for name in linnea:
      set_mutex(name, d)
      d -= sublinnean_increment
    mutex -= linnean_increment
  return (name_to_mutex_table, mutex_to_name_table)

(name_to_mutex_table, mutex_to_name_table) = \
  process_ranks(rank_configuration)

def name_to_mutex(name):
  # If absent, return the identity for min
  return name_to_mutex_table.get(name)

def mutex_to_name(mutex):       # For display purposes
  # If absent, return the numeric mutex I guess.  Can do better.
  return mutex_to_name_table.get(mutex, mutex)

root = name_to_mutex("root")

# ---

def self_test():
  for rank in ["subspecies", "atom", "forest"]:
    print("mutex %s -> %s" % (rank, name_to_mutex(rank)))

if __name__ == '__main__': self_test()
