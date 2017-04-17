import sys

def main():
  script, filename = sys.argv
  thelist = list()
  with open(filename, 'r') as f:
    for line in f:
      thelist.append(line)
    thelist = sorted(thelist, key=str.lower)

  with open("organized" + filename, 'w') as f:
    for name in thelist:
      f.write(name)

if __name__ == "__main__": main()
