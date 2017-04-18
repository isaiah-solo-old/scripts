import sys

def main():

  # If incorrect number of arguments
  if len(sys.argv) != 2:
    print "Usage: python sort.py <file>"
    sys.exit(0)

  # Read arguments
  script, filename = sys.argv

  # Create new list to store file contents
  thelist = []

  # Open file and check if valid
  try:
    # Copy file contents to list
    with open(filename, 'r') as f:
      for line in f:
        thelist.append(line)
      thelist = sorted(thelist, key=str.lower)

  # If file is not valid
  except IOError:
    print filename + " not a valid file"
    sys.exit(0)

  # Create new file to store edited list
  newfilename = "sorted_" + filename
  with open(newfilename, 'w') as f:
    for name in thelist:
      f.write(name)

    # If all passes, give success message
    print newfilename + " was successfully created"

if __name__ == "__main__": main()
