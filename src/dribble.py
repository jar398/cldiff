
import sys
import checklist as cl

dribble_file = None

def log(message):
  print(message)
  if dribble_file:
    print(message, file=dribble_file)

#confusing = 'Loris tardigradus'
confusing = 'Galagoides'

def watch(node):
  if node != None and node != cl.forest_tnu:
    return cl.get_name(node).startswith(confusing)
  return False
