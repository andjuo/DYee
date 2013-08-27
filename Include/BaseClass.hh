#ifndef BaseClass_HH
#define BaseClass_HH

#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <string>

// ------------------------------------------------------------
// ------------------------------------------------------------

class BaseClass_t {
  std::string className;
public:
  BaseClass_t(const std::string &theClassName) : className(theClassName) {
    if (className.size()==0) className="Base_t";
  }

protected:
  void printWarning(const char *msg) const {
    printf("warning in %s:",className.c_str());
    if (msg) std::cout << msg;
    std::cout << std::endl;
  }

  void printError(const char *msg) const {
    printf("error in %s:",className.c_str());
    if (msg) std::cout << msg;
    std::cout << std::endl;
  }

  template<class T>
  void printError(const char *format, const T &var) const {
    printf("error in %s:",className.c_str());
    if (format) printf(format,var);
    std::cout << std::endl;
  }

  template<class T1, class T2>
  void printError(const char *format, const T1 &var1, const T2 &var2) const {
    printf("error in %s:",className.c_str());
    if (format) printf(format,var1,var2);
    std::cout << std::endl;
  }

  template<class T1, class T2, class T3>
  void printError(const char *format, const T1 &var1, const T2 &var2, const T3 &var3) const {
    printf("error in %s:",className.c_str());
    if (format) printf(format,var1,var2,var3);
    std::cout << std::endl;
  }

  void printError(const char *format, const TString &var) const {
    TString tmp=(var.Length()) ? var : " ";
    printError(format,tmp.Data());
  }

  void printError(const char *format, const TString &var1, const TString &var2) const {
    TString tmp1=(var1.Length()) ? var1 : " ";
    TString tmp2=(var2.Length()) ? var2 : " ";
    printError(format,tmp1.Data(),tmp2.Data());
  }

  template<class T>
  void printError(const char *format, const TString &var1, const T &var2) const {
    TString tmp1=(var1.Length()) ? var1 : " ";
    printError(format,tmp1.Data(),var2);
  }

  template<class T>
  void printError(const char *format, const T &var1, const TString &var2) const {
    TString tmp2=(var2.Length()) ? var2 : " ";
    printError(format,var1,tmp2.Data());
  }

  int reportError(const char *msg) const {
    printError(msg);
    return 0;
  }

  template<class T>
  int reportError(const char *format, const T &var) const {
    printError(format,var);
    return 0;
  }

  template<class T1, class T2>
  int reportError(const char *format, const T1 &var1, const T2 &var2) const {
    printError(format,var1,var2);
    return 0;
  }

  template<class T1, class T2, class T3>
  int reportError(const char *format, const T1 &var1, const T2 &var2, const T3 &var3) const {
    printError(format,var1,var2,var3);
    return 0;
  }
  
  /*
  int reportError(const char *format, const TString &var) const {
    TString tmp=(var.Length()) ? var : " ";
    printError(format,tmp.Data());
    return 0;
  }
  */

};

// ------------------------------------------------------------
// ------------------------------------------------------------

#ifndef MYTOOLS_HH
//------------------------------------------------------------------------------------------------------------------------

inline
void HERE(const char *msg) {
  std::cout << ((msg) ? msg : "HERE") << std::endl;
}

//--------------------------------------------------------------

template<class Type_t>
inline
void HERE(const char *format, const Type_t &a) {
  std::cout << Form(format,a) << std::endl;
}

//--------------------------------------------------------------

template<class Type1_t, class Type2_t>
inline
void HERE(const char *format, const Type1_t &a, const Type2_t &b) {
  std::cout << Form(format,a,b) << std::endl;
}

//--------------------------------------------------------------

inline
void HERE(const char *format, const TString &a) {
  const char *aStr=(a.Length()) ? a.Data() : " ";
  std::cout << Form(format,aStr) << std::endl;
}
//--------------------------------------------------------------

inline
void HERE(const char *format, const TString &a, const TString &b) {
  const char *aStr=(a.Length()) ? a.Data() : " ";
  const char *bStr=(b.Length()) ? b.Data() : " ";
  std::cout << Form(format,aStr,bStr) << std::endl;
}

//--------------------------------------------------------------

template<class Type_t>
inline
void HERE(const char *format, const TString &a, const TString &b, const Type_t &x ) {
  const char *aStr=(a.Length()) ? a.Data() : " ";
  const char *bStr=(b.Length()) ? b.Data() : " ";
  std::cout << Form(format,aStr,bStr,x) << std::endl;
}

#endif

// ------------------------------------------------------------

#endif
