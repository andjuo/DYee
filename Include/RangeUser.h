#ifndef RangeUser_H
#define RangeUser_H

#include "../Include/ComparisonPlot.hh"
#include "../Include/CPlot.hh"


// --------------------------------------

struct RangeUser_t {
  double xmin,xmax;
  double ymin,ymax;
  int logX,logY;
public:
  RangeUser_t(double set_xmin=0., double set_xmax=-1.,
	      double set_ymin=0., double set_ymax=-1.,
	      int set_logX=-1, int set_logY=-1) :
    xmin(set_xmin), xmax(set_xmax),
    ymin(set_ymin), ymax(set_ymax),
    logX(set_logX), logY(set_logY)
  {}

  void apply(CPlot &cp) {
    if (xmax>xmin) cp.SetXRange(xmin,xmax);
    if (ymax>ymin) cp.SetYRange(ymin,ymax);
    if (logX!=-1) cp.SetLogx(logX);
    if (logY!=-1) cp.SetLogy(logY);
  }

  void apply(ComparisonPlot_t &cp) {
    if (xmax>xmin) cp.SetXRange(xmin,xmax);
    if (ymax>ymin) cp.SetYRange(ymin,ymax);
    if (logX!=-1) cp.SetLogx(logX);
    if (logY!=-1) cp.SetLogy(logY);
  }

  friend std::ostream& operator<<(std::ostream& out, const RangeUser_t &r) {
    out << "RangeUser(x: " << r.xmin << "," << r.xmax 
	<< "; y: " << r.ymin << "," << r.ymax 
	<< "; log: " << r.logX << "," << r.logY << ")";
    return out;
  }

};

// -------------------------------------

#endif
