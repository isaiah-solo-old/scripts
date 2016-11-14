#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>

using namespace std;

int main(int argc, char** argv) {
  if (argc != 2) {
    cout << "Usage: organize <filename>" << endl;
    return 1;
  }

  string line;
  vector<string> lines;
  ifstream infile (argv[1]);

  if (! infile.is_open()) {
    cout << "Error opening \'" << argv[1] << "\'" << endl;
    return 1;
  }

  while (getline (infile, line)) {
    lines.push_back(line);
  }
  infile.close();
  remove(argv[1]);

  sort(lines.begin(), lines.end());

  ofstream outfile;
  outfile.open(argv[1]);
  for (vector<string>::iterator itor = lines.begin(); itor != lines.end(); ++itor) {
    outfile << *itor << endl;
  }
  outfile.close();

  lines.clear();
}
