# Properties

import collections
import dribble

# Property: record -> value

Property = \
  collections.namedtuple('Property',
                         ['uid',
                          'uri',
                          'pet_name',     # local to this code base
                          'specificity']) # with regard to taxon identity

properties_by_uri = {}
properties_by_pet_name = {}
properties_by_specificity = []

def by_specificity(specificity):
  return properties_by_specificity[specificity]  

def by_name(name):
  assert not "\t" in name
  prop = properties_by_pet_name.get(name)
  if not prop:
    dribble.log("No such property: %s" % name)
  return prop

def uri_to_pet_name(uri):
  return uri.split('/')[-1]

def uri_to_property(uri):
  probe = properties_by_uri.get(uri)
  if probe: return probe
  pet_name = uri_to_pet_name(uri)
  assert pet_name
  uid = len(properties_by_specificity)
  spec = uid
  sel = Property(uid, uri, pet_name, spec)
  properties_by_uri[uri] = sel
  properties_by_pet_name[pet_name] = sel
  properties_by_specificity.append(sel)
  return sel

def alias(this, to_that):
  if this in properties_by_pet_name:
    raise "URI pet_name collision:  \n%s  \n%s" % (to_that, this)
  sel = by_name(to_that)
  properties_by_pet_name[this] = sel

# Initialize non-field properties

# We do this in order from least specific to most specific,
# for simplicity... low order bits first

# dwciri: (http://rs.tdwg.org/dwc/iri/)

# The entire record for the taxon in the checklist
record    = uri_to_property("data:,property/record")

# dwc:TaxonID - "This term is no longer recommended for use."
#   http://rs.tdwg.org/dwc/terms

# URIs (columns) I'm considering adding:
#   - one URI for each table I set up in the analysis
#     - something for what I'm currently calling record id
#       (unique record key across all inputs)...
#     - mutexes

"""
GBIF properties...
  number_of_synonyms
  datasetID
  originalNameUsageID
  genericName
  nameAccordingTo
  namePublishedIn
  taxonRemarks
  kingdom
  phylum
  class
  order
  family
  genus
"""

uris = [
        "http://purl.org/dc/terms/source",     # EOL
        "https://github.com/jar398/biodiversity/wiki/subtree",  # subtree identity
        "https://github.com/jar398/biodiversity/wiki/tips",     # same tips
        "https://github.com/jar398/biodiversity/wiki/direct",   # not via synonym
        "http://rs.tdwg.org/dwc/terms/taxonomicStatus",     # flush?
        "http://rs.tdwg.org/dwc/terms/nomenclaturalStatus",
        "https://github.com/jar398/biodiversity/wiki/number_of_children",
        "https://github.com/jar398/biodiversity/wiki/number_of_synonyms",
        "http://rs.tdwg.org/dwc/terms/verbatimTaxonRank",   # don't use this URI...
        "http://rs.tdwg.org/dwc/terms/taxonRank",
        "http://rs.tdwg.org/dwc/terms/scientificNameAuthorship",
        "http://rs.tdwg.org/dwc/terms/nameAccordingToID",
        "http://rs.tdwg.org/dwc/terms/taxonID",    # exclude from diff?
        "https://github.com/jar398/biodiversity/wiki/gbif_id",
        "https://github.com/jar398/biodiversity/wiki/ncbi_id",
        "http://eol.org/schema/EOLid",  # page ID (via Jen)
        "http://rs.tdwg.org/dwc/terms/vernacularName",
        "http://rs.tdwg.org/dwc/terms/parentNameUsageID",
        "http://rs.tdwg.org/dwc/terms/namePublishedInYear",
        "http://rs.tdwg.org/dwc/terms/family",    # GBIF
        "http://rs.tdwg.org/dwc/terms/genus",    # GBIF
        "http://rs.tdwg.org/dwc/terms/specificEpithet",
        "http://rs.tdwg.org/dwc/terms/infraspecificEpithet",
        "https://github.com/jar398/biodiversity/wiki/quasi_epithet",
        "http://rs.tdwg.org/dwc/terms/acceptedNameUsageID",
        "http://rs.gbif.org/terms/1.0/canonicalName",
        "http://rs.tdwg.org/dwc/terms/scientificName",
        "https://github.com/jar398/biodiversity/wiki/particles",
]
for uri in uris: uri_to_property(uri)
  
alias("NCBI Taxonomy ID", "ncbi_id")

number_of_properties = len(properties_by_specificity)

if __name__ == '__main__':
  import sys
  import yaml
  with open(sys.argv[1], 'r') as stream:
    print(yaml.load(stream))
