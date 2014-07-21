#ifndef TESTHELPERROOUNFOLDBASISSPLINES_H
#define TESTHELPERROOUNFOLDBASISSPLINES_H

#include "RooUnfoldBasisSplines.h"

//class to access the private elements of RooUnfoldBasisSplines for testing

class testHelperRooUnfoldBasisSplines
{
  const RooUnfoldBasisSplines* ourFriend;

 public:
 testHelperRooUnfoldBasisSplines(const RooUnfoldBasisSplines* fr)
   : ourFriend(fr) {return;};

  ~testHelperRooUnfoldBasisSplines() {return;};

  Double_t GetTau() const {return ourFriend->_tau;};
  Int_t GetM0() const {return ourFriend->_m0;};
  Int_t GetIauto() const {return ourFriend->_iauto;};

};

#endif// TESTHELPERROOUNFOLDBASISSPLINES_H
