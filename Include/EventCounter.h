#ifndef EventCounter_H
#define EventCounter_H

#include <TROOT.h>
#include <TString.h>
#include <iostream>


struct EventCounter_t {
public:
  TString name;
  ULong_t numEvents;
  ULong_t numEventsInGenAcceptance;
  ULong_t numEventsPassedEvtTrigger;
  ULong_t numDielectronsUnweighted;
  //ULong_t numPassedGoodPV;
  // weighted counts
  double numDielectrons;
  double numDielectronsGenMatched;
  double numDielectronsGoodEta;
  double numDielectronsGoodEt;
  double numDielectronsGoodEtEta;
  double numDielectronsHLTmatched;
  double numDielectronsIDpassed;
  double numDielectronsGoodMass;
  double numDielectronsPass[2];
  double numMultiDielectronsOk;
  double scale;
  int ignoreScale; // whether scale change should be allowed
 public:
  EventCounter_t(TString set_name) :
    name(set_name),
    numEvents(0), 
    numEventsInGenAcceptance(0),
    numEventsPassedEvtTrigger(0),
    numDielectronsUnweighted(0),
    //numPassedGoodPV(0),
    numDielectrons(0), 
    numDielectronsGenMatched(0),
    numDielectronsGoodEta(0), numDielectronsGoodEt(0),
    numDielectronsGoodEtEta(0),
    numDielectronsHLTmatched(0),
    numDielectronsIDpassed(0),
    numDielectronsGoodMass(0),
    numMultiDielectronsOk(0),
    scale(1.),
    ignoreScale(0)
  {
    numDielectronsPass[0]=0.;
    numDielectronsPass[1]=0.;
  }

  EventCounter_t(const EventCounter_t &e) {
    this->assign(e);
  }

  void clear(TString set_name="") {
    name=set_name;
    numEvents=0; 
    numEventsInGenAcceptance=0;
    numEventsPassedEvtTrigger=0;
    numDielectronsUnweighted=0;
    numDielectrons=0;
    numDielectronsGenMatched=0;
    numDielectronsGoodEta=0; numDielectronsGoodEt=0;
    numDielectronsGoodEtEta=0;
    numDielectronsHLTmatched=0;
    numDielectronsIDpassed=0;
    numDielectronsGoodMass=0;
    numDielectronsPass[0]=0.;
    numDielectronsPass[1]=0.;
    numMultiDielectronsOk=0,
    scale=1.;
  }

  TString &editName() { return name; }

  void numEvents_inc() { numEvents++; }
  void numEventsPassedAcceptance_inc() { numEventsInGenAcceptance++; }
  void numEventsPassedEvtTrigger_inc() { numEventsPassedEvtTrigger++; }
  void numDielectrons_inc() { numDielectrons+=scale; numDielectronsUnweighted++; }
  void numDielectronsGenMatched_inc() { numDielectronsGenMatched+=scale; }
  void numDielectronsGoodEta_inc() { numDielectronsGoodEta+=scale; }
  void numDielectronsGoodEt_inc() { numDielectronsGoodEt+=scale; }
  void numDielectronsGoodEtEta_inc() { numDielectronsGoodEtEta+=scale; }
  void numDielectronsHLTmatched_inc() { numDielectronsHLTmatched+=scale; }
  void numDielectronsIDpassed_inc() { numDielectronsIDpassed+=scale; }
  void numDielectronsGoodMass_inc() { numDielectronsGoodMass+=scale; }
  void numMultiDielectronsOk_inc() { numMultiDielectronsOk+=scale; }

  void numDielectronsPass_inc() { 
    numDielectronsPass[0] += scale;
    numDielectronsPass[1] += scale*scale;
  }

  void setIgnoreScale(int ignore) { ignoreScale=ignore; if (ignore) scale=1.; }
  void setScale(double sc) { if (!ignoreScale) scale=sc; }

  void assign(const EventCounter_t &e, TString new_name="") {
    if (new_name.Length()==0) name=e.name+TString("_clone");
    else name=new_name;
    numEvents=e.numEvents;
    numEventsInGenAcceptance=e.numEventsInGenAcceptance;
    numEventsPassedEvtTrigger=e.numEventsPassedEvtTrigger;
    numDielectronsUnweighted=e.numDielectronsUnweighted;
    numDielectrons=e.numDielectrons;
    numDielectronsGenMatched=e.numDielectronsGenMatched;
    numDielectronsGoodEta=e.numDielectronsGoodEta;
    numDielectronsGoodEt=e.numDielectronsGoodEt;
    numDielectronsGoodEtEta=e.numDielectronsGoodEtEta;
    numDielectronsHLTmatched=e.numDielectronsHLTmatched;
    numDielectronsIDpassed=e.numDielectronsIDpassed;
    numDielectronsGoodMass=e.numDielectronsGoodMass;
    numDielectronsPass[0]=e.numDielectronsPass[0];
    numDielectronsPass[1]=e.numDielectronsPass[1];
    numMultiDielectronsOk=e.numMultiDielectronsOk;
    scale=e.scale;
  }

  void add(const EventCounter_t &e) {
    numEvents+=e.numEvents;
    numEventsInGenAcceptance+=e.numEventsInGenAcceptance;
    numEventsPassedEvtTrigger+=e.numEventsPassedEvtTrigger;
    numDielectronsUnweighted+=e.numDielectronsUnweighted;
    numDielectrons+=e.numDielectrons;
    numDielectronsGenMatched+=e.numDielectronsGenMatched;
    numDielectronsGoodEta+=e.numDielectronsGoodEta;
    numDielectronsGoodEt+=e.numDielectronsGoodEt;
    numDielectronsGoodEtEta+=e.numDielectronsGoodEtEta;
    numDielectronsHLTmatched+=e.numDielectronsHLTmatched;
    numDielectronsIDpassed+=e.numDielectronsIDpassed;
    numDielectronsGoodMass+=e.numDielectronsGoodMass;
    numDielectronsPass[0]+=e.numDielectronsPass[0];
    numDielectronsPass[1]+=e.numDielectronsPass[1];
    numMultiDielectronsOk+=e.numMultiDielectronsOk;
    if (scale!=e.scale) scale=-1;
  }

  void print(int level=0) const { print(std::cout,level); }

  void print(std::ostream &out, int level=0) const  {
    const char *line="-----------------------------------------------\n";
    const char *format="%9.4lf";
    out << line;
    out << "eventCounter " << this->name << " info (level=" << level << "):\n";
    out << Form("   numEvents                = %lu\n",this->numEvents);
    out << Form("   numEventsInGenAcceptance = %lu\n",this->numEventsInGenAcceptance);
    if (level<=1) {
      out << Form("   numEventsPassedEvtTrigger= %lu\n",this->numEventsPassedEvtTrigger);
      out << Form("   numDielectronsUnweighted = %lu\n",this->numDielectronsUnweighted);
    }
    if (level==0) {
      out << "   numDielectrons           = " << Form(format,this->numDielectrons) << "\n";
      if (this->numDielectronsGenMatched) {
	out << "   numDielectronsGenMatched = " << Form(format,this->numDielectronsGenMatched) << "\n";
      }
      out << "   numDielectronsGoodEta    = " << Form(format,this->numDielectronsGoodEta) << "\n";
      out << "   numDielectronsGoodEt     = " << Form(format,this->numDielectronsGoodEt) << "\n";
      out << "   numDielectronsGoodEtEta  = " << Form(format,this->numDielectronsGoodEtEta) << "\n";
      out << "   numDielectronsHLTmatched = " << Form(format,this->numDielectronsHLTmatched) << "\n";
      out << "   numDielectronsIDpassed   = " << Form(format,this->numDielectronsIDpassed) << "\n";
      out << "   numDielectronsGoodMass   = " << Form(format,this->numDielectronsGoodMass) << "\n";
      out << "   numDielectronsPass       = " << Form(format,this->numDielectronsPass[0]) << " +/- " << Form(format,sqrt(this->numDielectronsPass[1])) << "\n";
    }
    out   << "  numMultiDielectronsOk     = " << Form(format,this->numMultiDielectronsOk) << "\n";
    if (this->ignoreScale) out << "   scale ignored\n";
    else out << Form("   last scale = %9.4e\n",this->scale);
    out << line;
  }

  friend std::ostream& operator<<(std::ostream& out, const EventCounter_t &e) {
    e.print(out);
    return out;
  }
};

// --------------------------------------------------

struct EventCounterExt_t : public EventCounter_t {
  double numDielectronsOkPosSS;
  double numDielectronsOkNegSS;
  //double numDielectronsZPeak;
public:
  EventCounterExt_t(TString set_name="") : EventCounter_t(set_name),
			numDielectronsOkPosSS(0.),
			numDielectronsOkNegSS(0.)
  {}

  EventCounterExt_t(const EventCounter_t &e) : EventCounter_t("") {
    std::cout << "cannot call constructor EventCounterExt_t(EventCounter_t)\n";
    std::cout << "  --\n" << e << "\n";
  }

  EventCounterExt_t(const EventCounterExt_t &e) :
    EventCounter_t(e),
    numDielectronsOkPosSS(e.numDielectronsOkPosSS), 
    numDielectronsOkNegSS(e.numDielectronsOkNegSS) {
  }

  void clear(TString new_name="") {
    EventCounter_t::clear(new_name);
    numDielectronsOkPosSS=0.;
    numDielectronsOkNegSS=0.;
  }

  void numDielectronsOkSameSign_inc(int sign) {
    if (sign>0) numDielectronsOkPosSS+=scale;
    else numDielectronsOkNegSS+=scale;
  }

  int numDielectronsOkSameSign_inc(int q1, int q2) {
    if (q1!=q2) return 0;
    if (q1>0) numDielectronsOkPosSS+=scale;
    else numDielectronsOkNegSS+=scale;
    return 1;
  }

  void assign(const EventCounter_t &e, TString new_name="") {
    if (0) { std::cout << new_name; e.print(); } // avoid compiler warning
    std::cout << "cannot call EventCounterExt_t::assign(EventCounter_t)\n";
  }

  void assign(const EventCounterExt_t &e, TString new_name="") {
    EventCounter_t::assign(e,new_name);
    numDielectronsOkPosSS=e.numDielectronsOkPosSS;
    numDielectronsOkNegSS=e.numDielectronsOkNegSS;
  }

  void add(const EventCounter_t &e) {
    if (0) e.print(); // avoid compiler warning
    std::cout << "cannot call eventCounter_t::add(eventCounter_t)\n";
  }

  void add(const EventCounterExt_t &e) {
    EventCounter_t::add(e);
    numDielectronsOkPosSS+=e.numDielectronsOkPosSS;
    numDielectronsOkNegSS+=e.numDielectronsOkNegSS;
  }

  void print(int level=0) const { print(std::cout,level); }

  void print(std::ostream &out, int level=0) const { 
    const char *line="-----------------------------------------------\n";
    const char *format="%9.4lf";
    out << line;
    out << "eventCounter " << this->name << " info (level=" << level << "):\n";
    out << Form("   numEvents                = %lu\n",this->numEvents);
    out << Form("   numEventsInGenAcceptance = %lu\n",this->numEventsInGenAcceptance);
    if (level<=1) {
    out << Form("   numEventsPassedEvtTrigger= %lu\n",this->numEventsPassedEvtTrigger);
    }
    if (level==0) {
      out << Form("   numDielectronsUnweighted = %lu\n",this->numDielectronsUnweighted);
      out << "   numDielectrons           = " << Form(format,this->numDielectrons) << "\n";
      if (this->numDielectronsGenMatched) {
	out << "   numDielectronsGenMatched = " << Form(format,this->numDielectronsGenMatched) << "\n";
      }
      out << "   numDielectronsGoodEta    = " << Form(format,this->numDielectronsGoodEta) << "\n";
      out << "   numDielectronsGoodEt     = " << Form(format,this->numDielectronsGoodEt) << "\n";
      out << "   numDielectronsGoodEtEta  = " << Form(format,this->numDielectronsGoodEtEta) << "\n";
      out << "   numDielectronsHLTmatched = " << Form(format,this->numDielectronsHLTmatched) << "\n";
      out << "   numDielectronsIDpassed   = " << Form(format,this->numDielectronsIDpassed) << "\n";
      out << "   numDielectronsGoodMass   = " << Form(format,this->numDielectronsGoodMass) << "\n";
      out << "   numDielectronsPass       = " << Form(format,this->numDielectronsPass[0]) << " +/- " << Form(format,sqrt(this->numDielectronsPass[1])) << "\n";
      out << "   numDielectronsOkPosSS    = " << Form(format,this->numDielectronsOkPosSS) << "\n";
      out << "   numDielectronsOkNegSS    = " << Form(format,this->numDielectronsOkNegSS) << "\n";
    }
    out   << "  numMultiDielectronsOk     = " << Form(format,this->numMultiDielectronsOk) << "\n";
    if (this->ignoreScale) out << "   scale ignored\n";
    else out << Form("   last scale = %9.4e\n",this->scale);
    out << line;
  }

  friend std::ostream& operator<<(std::ostream& out, const EventCounterExt_t &e) {
    e.print(out);
    return out;
  }

};


// --------------------------------------------------

#endif
