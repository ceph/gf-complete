#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
using namespace std;

#define VIT(i, v) for (i = 0; i < v.size(); i++) 
#define IT(it, ds) for (it = ds.begin(); it != ds.end(); it++)
#define FUP(i, n) for (i = 0; i < n; i++)

typedef map<int, string> ISmap;
typedef map<int, int> IImap;
typedef map<string, double> SDmap;

typedef ISmap::iterator ISmit;
typedef IImap::iterator IImit;
typedef SDmap::iterator SDmit;

typedef vector <string> SVec;

void StoSVec(string &s, SVec &sv)
{
  istringstream ss;
  string s2;

  ss.clear();
  ss.str(s);
  while (ss >> s2) sv.push_back(s2);
}

main()
{
  string s, k;
  double d, b;
  int i;
  SVec sv;
  SDmap bmap;
  SDmit bmit;

  while (getline(cin, s)) {
    sv.clear();
    StoSVec(s, sv);
    
    if (sv[0] == "Seed:") {
      b = 0;
      for (i = 0; i < 2; i++) {
        getline(cin, s);
        sv.clear();
        StoSVec(s, sv);
        sscanf(sv[3].c_str(), "%lf", &d);
        if (d > b) b = d;
      }
      getline(cin, s);
      sv.clear();
      StoSVec(s, sv);
      k = sv[2];
      k += " ";
      k += sv[3];
      for (i = 4; i < sv.size(); i++) {
        if (sv[i] != "-") {
          k += " ";
          k += sv[i];
        }
      }
      if (bmap[k] < b) bmap[k] = b;
    }
  }

  IT(bmit, bmap) {
    printf("%10.4lf %s\n", bmit->second, bmit->first.c_str());
  }
}
