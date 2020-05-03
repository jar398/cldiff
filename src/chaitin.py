# This is for testing the alignment logic.

import sys, csv

# Parses head, returns rest of string
#   isinstance(x, list)

def record(id, name, children, parent_id):
  id = str(id)
  if parent_id != None:
    parent_id = str(parent_id)
  return [[str(id), name, parent_id]] + children

def get_name(row): return row[1]

# Returns a list of records

def chaitin_head(s, i, parent):
  if s[i] == '(':
    (lst, j) = chaitin_list(s, i+1, i)
    head = lst[0]
    tail = lst[1:]
    return (record(i, get_name(head), tail, parent), j)
  else:
    return (record(i, s[i], [], parent), i+1)

# Returns a list of records

def chaitin_list(s, i, parent):
  if i == len(s):
    return ([], i)
  elif s[i] == ')':
    return ([], i+1)
  else:
    (head, j) = chaitin_head(s, i, parent)
    (tail, k) = chaitin_list(s, j, parent)
    return (head + tail, k)

def parse(s):
  yield ["taxonID", "canonicalName", "parentNameUsageID"]
  (lst, i) = chaitin_list(s, 0, None)
  for rec in lst:
    yield rec

def dump(s, path):
  rows = parse(s)
  if path == "-":
    dump_chaitin(rows, sys.stdout)
  else:
    with open(path, "w") as outfile:
      dump_chaitin(rows, outfile)

def dump_chaitin(rows, io):
  writer = csv.writer(io)
  for row in rows:
    writer.writerow(row)

def self_test():

  def test(x):
    print ("%s:" % x)
    for rec in parse(x):
      print (rec)
    print

  test("abc")
  test("(a)bc")
  test("(az)bc")
  test("q(az)bc")
  test("qr(az)bc")
  test("(((a)))")
  test("(a(b(c)))")

if __name__ == '__main__':
  self_test()
  dump("(a(bc))((ab)c)", '-')
