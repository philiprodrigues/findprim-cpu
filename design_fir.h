#ifndef DESIGN_FIR_H
#define DESIGN_FIR_H

#include <vector>
#include <cmath>

// Functions to design a lowpass filter with the given number of taps
// and cutoff, as a fraction of the Nyquist frequency. Copied from the
// corresponding functions in scipy and converted to C++

#define PI 3.1416

std::vector<double> hamming(int M)
{
  std::vector<double> ret;
  for(int n=0; n<M; ++n){
    ret.push_back(0.54-0.46*cos(2.0*PI*n/(M-1)));
  }
  return ret;
}

double sinc(double x)
{
  if(x==0) return 1;
  return sin(PI*x)/(PI*x);
}

std::vector<double> firwin(int N, double cutoff)
{
  int alpha = N/2;
  std::vector<double> window=hamming(N);
  std::vector<double> ret;
  double sum=0;
  for(int m=0; m<N; ++m){
    double val=window[m]*sinc(cutoff*(m-alpha));
    ret.push_back(val);
    sum+=val;
  }
  for(int m=0; m<N; ++m){ ret[m]/=sum; }

  return ret;
}

#endif
