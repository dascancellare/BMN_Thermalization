#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <fstream>

using namespace std;

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

void internal_crash(int line,const char *file,const char *temp,...);

string combine(const char *format,...);

//! read an element from input file
template <class T> void read(T &out,ifstream &in,string is)
{
  string s;
  if(!(in>>s)) CRASH("impossible to read expect string \"%s\"",is.c_str());
  if(s!=is) CRASH("obtained %s while reading %s",s.c_str(),is.c_str());
  if(!(in>>out)) CRASH("reading data");
}

int file_exists(string path);

#endif
