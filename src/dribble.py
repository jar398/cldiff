
import sys

dribble_file = None

def log(message):
  print(message)
  if dribble_file:
    print(message, file=dribble_file)
