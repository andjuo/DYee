#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <iostream>

void printField(TString fname="", TString fieldType="M", TString fieldName="")
{
  TFile f(fname,"read");
  if (!f.IsOpen()) {
    std::cout << "failed to open file <" << fname << ">\n";
    return;
  }
  std::cout << fname << " file contents:\n";
  f.ls();

  std::cout << "\n\n";
  if (fieldName.Length()) {
    if (fieldType==TString("M")) {
      TMatrixD* m=(TMatrixD*)f.Get(fieldName);
      if (!m) {
	std::cout << "failed to get the matrix field <" << fieldName << ">\n";
      }
      else {
	std::cout << "matrix field " << fieldName << "\n";
	m->Print();
	delete m;
      }
    }
    else if (fieldType==TString("V")) {
      TVectorD* m=(TVectorD*)f.Get(fieldName);
      if (!m) {
	std::cout << "failed to get the vector field <" << fieldName << ">\n";
      }
      else {
	std::cout << "vector field " << fieldName << "\n";
	m->Print();
	delete m;
      }
    }
  }
  f.Close();
}
