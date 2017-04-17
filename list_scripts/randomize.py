import sys
import random

def main():
  script, filename = sys.argv
  thelist = list()
  with open(filename, 'r') as f:
    for line in f:
      thelist.append(line)
    random.shuffle(thelist)

  with open("random" + filename, 'w') as f:
    for name in thelist:
      f.write(name)

if __name__ == "__main__": main()
