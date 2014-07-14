/* compile without generating an executable: g++ -std=c++11 -I/opt/local/include -c tasep.cpp */
#include<iostream>
#include<fstream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include "tasep_ode.hpp"
#include "math_utils.hpp"

using boost::numeric::odeint::integrate;
void tasep::operator() (const state_type &p, state_type &dpdt, double t)
{
  // length of rate vector is one element more than p vector
  size_t n(rate.size()-1);
  dpdt[0] = rate[0]*(1-p[0]) - rate[1]*p[0]*(1-p[1]);
  for (size_t i=1; i!=n-1; ++i) {
    dpdt[i] = rate[i]*p[i-1]*(1-p[i]) - rate[i+1]*p[i]*(1-p[i+1]);
  }
  dpdt[n-1] = rate[n-1]*p[n-2]*(1-p[n-1]) - rate[n]*p[n-1];
}

double tasep::numeric_steady_state(result_buffer& observer, state_type& p, double t0, double t1, double steady_threshold)
{
  double residue_norm(1000);
  size_t steps=0;
  while (residue_norm >= steady_threshold) {
    steps += integrate(*this, p, t0, t1, 0.1, observer);
    state_type residue(p.size(),0);
    (*this)(p,residue, 0);
    residue_norm = l2norm(residue);
    if (t1>10000) return residue_norm;
    t0 = t1;
    t1 *= 2;
  }
  double psum=0;
  for (auto pi: p) psum += pi;
  //time, steps, residue, sum of p, peptide length
  std::cout<<t1<<"\t"<<steps<<"\t"<<residue_norm<<"\t"<<psum<<"\t"<<p.size()<<std::endl;
  return residue_norm;
}

double tasep::numeric_steady_state(state_type& p, double t0, double t1, double steady_threshold)
{
  double residue_norm(1000);
  size_t steps(0);
  while (residue_norm >= steady_threshold) {
    steps += integrate(*this, p, t0, t1, 0.1);
    state_type residue(p.size(),0);
    (*this)(p,residue, 0);
    residue_norm = l2norm(residue)/residue.size();
    if (t0>1e4) {
      //std::cout<<"time: "<<t0<<" steps: "<<steps<<" flow diff per position: "<<residue_norm<<" plen: "<<p.size()<<std::endl;
      return residue_norm;
      //std::cout<<residue_norm<<" "<<std::flush;
    }
    t0 = t1;
    t1 *= 2;
  }
  //time, steps, residue, sum of p, peptide length
  //std::cout<<"time: "<<t0<<" steps: "<<steps<<" flow diff: "<<residue_norm<<" plen: "<<p.size()<<std::endl;
  return residue_norm;
}

void result_buffer::operator() (const state_type &p, double t)
{
  plist.push_back(p);
  tlist.push_back(t);
}

void result_buffer::write_to_file(const char* fname)
{
  std::ofstream ofile(fname);
  for (size_t i = 0; i!=tlist.size(); ++i){
    ofile<<tlist[i]<<"\t";
    for (auto p: plist[i])
      ofile<<p<<"\t";
    ofile<<std::endl;
  }
  ofile.close();
}
