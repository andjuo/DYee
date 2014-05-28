#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


//------------------------------------------------------------------------
#ifndef MyTools_HH

inline bool PosOk(size_t pos) { return (pos==std::string::npos) ? 0 : 1; }

//------------------------------------------------------------------------

template<class T>
inline bool PosOk(const std::string &s, const T& substr) {
  return (s.find(substr)==std::string::npos) ? 0 : 1;
}

#endif
//----------------------------------------------------------------------
//----------------------------------------------------------------------

void prepareSigmaZfile() {
  TString inpFileName="dir-CSInfo/inpCS_sigmaZ.txt";
  TString outFileName="dir-CSInfo/inpCS_sigmaZ.dat";
  std::ifstream fin(inpFileName);
  std::ofstream fout(outFileName);
  std::string s;
  while (!fin.eof() && getline(fin,s)) {
    TString name, errStr;
    if (PosOk(s,"acc.")) {
      std::cout << "work on <" << s << ">\n";
      size_t pos2=s.find_first_of('&');
      std::string s2=s.substr(0,pos2-1);
      size_t pos=s2.find_last_not_of(' ');
      name=s2.substr(0,pos);
      name.ReplaceAll(" ","_");
      std::cout << "name=<" << name << ">\n";
      double centralVal=0., totErr=0., lumiErr=0.;
      const char *ptr=s.c_str();
      int digit=0, first_digit=1;
      for (size_t p=pos2+1; p<s.size(); ++p) {
	if (!digit) {
	  if ((ptr[p]>='0') && (ptr[p]<='9')) {
	    digit=1;
	    double term=atof(ptr+p);
	    if (!first_digit) {
	      lumiErr= term*term;
	      totErr += term*term;
	      std::cout <<" adding err=" << term << "\n";
	    }
	    else centralVal=term;
	    first_digit=0;
	  }
	}
	if (digit) {
	  if (((ptr[p]<'0') || (ptr[p]>'9')) && (ptr[p]!='.')) {
	    digit=0;
	    errStr.Append(' ');
	  }
	  else {
	    errStr.Append(ptr[p]);
	  }
	}
      }
      std::stringstream ss;
      std::cout << name.Length() << "\n";
      ss << name << std::string(23-name.Length(),' ');
      ss << " " << centralVal << " " << sqrt(totErr-lumiErr)
	 << "    " << errStr << "\n";
      std::cout << ss.str();
      fout << ss.str();
    }
  }
  fin.close();
  fout.close();
  std::cout << "file <" << outFileName << "> created\n";
  return;
}

//----------------------------------------------------------------------
