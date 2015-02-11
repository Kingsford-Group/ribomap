/**
	This file is part of the Ribomap suite
	-- an automatic pipeline for quantifying isoform-level
	ribosome profiles from ribosome profiling data.


	Copyright 2015 Hao Wang

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

	    http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
**/



#include <vector>
#include <array>
#include <string>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <thread>
#include <iostream>

#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "transcript_model.hpp"
#include "tasep_ode.hpp"

/***
 * class rate_model:
 * store model parameters, model computation for each transcript
 */
ribo_model::ribo_model(const char* tfa_fname, const ribo_profile& profiler, const transcript_info& tinfo): erate(vector<double>(61,0)), _refID2pID(profiler.refID2pID)
{
  vector<rid_t> refID(profiler.get_expressed_transcript_ids());
  build_codon_rateID_map();
  build_transcript_rate_id_vec(tfa_fname, refID, tinfo);
}

ribo_model::ribo_model(const char* tfa_fname, const vector<rid_t>& refIDs, const transcript_info& tinfo): erate(vector<double>(61,0))
{
  build_codon_rateID_map();
  build_transcript_rate_id_vec(tfa_fname, refIDs, tinfo);
}

void ribo_model::random_initiation_rate(double imin, double imax, long seed)
{
  mt19937_64 gen(seed);
  uniform_real_distribution<double> irate_sample(imin,imax);
  for (size_t i=0; i!=model.size(); ++i)
    set_irate(i,irate_sample(gen));
}

void ribo_model::uniform_initiation_rate(double irate)
{
  for (size_t i=0; i!=model.size(); ++i)
    set_irate(i,irate);
}

void ribo_model::initiation_rate_from_file(const char* irate_fname)
{
  ifstream ifile(irate_fname);
  rid_t refID, t;
  double rate;
  while (ifile >> refID >> rate) {
    t = _refID2pID.at(refID);
    set_irate(t, rate);
  }
  ifile.close();
}

void ribo_model::initiation_rate_to_file(const char* irate_fname)
{
  ofstream ofile(irate_fname);
  rid_t refID, t;
  double rate;
  for (auto p: _refID2pID) {
    refID = p.first;
    t = p.second;
    rate = get_irate(t);
    ofile<< refID << " " << rate<<endl;
  }
  ofile.close();
}

void ribo_model::elongation_rate_from_file(const char* erate_fname)
{
  ifstream ifile(erate_fname);
  string codon;
  double rate;
  while (ifile >> codon >> rate) {
    int eid = _codon2id[codon];
    set_erate(eid, rate);
  }
  ifile.close();
}

void ribo_model::update_profile(const ribo_profile& rprofile, bool ignore_no_map, int nproc)
{
  omp_set_num_threads(nproc);
  #pragma omp parallel for schedule(static)
  for (size_t t=0; t<rprofile.number_of_transcripts(); ++t) {
    if (ignore_no_map and rprofile.get_tot_count(t) == 0) continue;
    tmodel& m = get_model_ref(t);
    if (get_obj(t)<0) continue;
    vector<double> rate_vec(build_rate_vec(t));
    double trate = compute_profile(rate_vec, m.profile);
    if (trate>0)
      set_trate(t, trate);
    else
      compute_translation_rate(t, rate_vec);
    compute_transcript_abundance(t, rprofile.get_read_assignments(t));
  }
}

void ribo_model::compute_translation_rate(rid_t t, const vector<double>& rate)
{
  assert(rate.size() == codon_vec_len(t)+1);
  size_t n(rate.size());
  double trate_avg(rate[0]*(1-get_prob(t,0)));
  for (size_t i=1; i!=n-1; ++i)
    trate_avg += rate[i]*get_prob(t,i-1)*(1-get_prob(t,i));
  trate_avg += rate[n-1]*get_prob(t,n-2);
  trate_avg /= n;
  set_trate(t,trate_avg);
}

void ribo_model::compute_transcript_abundance(rid_t t, const vector<double>& count)
{
  vector<double> c(count);
  const vector<double>& p = get_profile(t);
  double a(0), b(0); //, c(0);
  for (size_t i=0; i!=p.size(); ++i) {
    a += std::pow(p[i], 2);
    b -= 2*p[i]*c[i];
    //c += std::pow(count[i],2);
  }
  double tcnt = -b/(2*a);
  assert(tcnt>=0);
  set_tcnt(t,tcnt);
}

// double ribo_model::fit_single_irate(rid_t t, const vector<double>& count_vec, bool normalized)
// {
//   //cout<<"fit_single_irate: "<<t<<" "<<std::this_thread::get_id()<<endl;
//   //if (t>10) exit(1);
//   tmodel& m = get_model_ref(t);
//   double& start_point = m.irate;
//   vector<double> r(build_rate_vec(t));
//   if (r.size() != m.profile.size() + 1) {
//     cout<<"in ribo_model::fit_irate, rate_vec.size != profile.size + 1 !"<<endl;
//     return 0;
//   }
//   vector<double> c(count_vec);
//   if (normalized) normalize_vec(c);
//   irate_finder optimizer(c, r, m);
//   double obj_end(0);
//   try {
//     obj_end = dlib::find_min_single_variable(optimizer, start_point, 1e-200, 1e200, 1e-3, 1000);
//   }
//   catch (dlib::optimize_single_variable_failure& e) {
//     if (e.info == "timeout")
//       set_obj(t,-2);
//     else
//       set_obj(t,-1);
//     cout<<" b"<<t<<"b "<<flush;
//     //cout<<"optimization for transcript "<<t<<" not converged!"<<endl;
//     return -1;
//   }
//   set_obj(t,obj_end);
//   //cout<<t<<" "<<flush;
//   //cout<<"initiation rate: "<<start_point<<" objective function value: "<<obj_end<<endl;
//   return obj_end;
// }

// void ribo_model::fit_irate_block(const ribo_profile& rprofile, rid_t t_offset, int num_rate)
// {
//   // cout<<t_offset<<endl;
//   //cout<<"fit_irate_block: "<<t_offset<<" "<<std::this_thread::get_id()<<endl;
//   for (size_t t=t_offset; t!=t_offset+num_rate; ++t) {
//     // // if read assignments not abundant enough, ignore
//     // if (rprofile.get_tot_count(t) < rprofile.len(t) ) {
//     //   set_obj(t,0);
//     //   continue;
//     // }
//     fit_single_irate(t,rprofile.get_read_assignments(t),true);
//   }
// }

// void ribo_model::fit_irates_in_parallel(const ribo_profile& rprofile, int num_threads)
// {
//   int block_size = max(int(model.size()/double(num_threads)),1);
//   int offset(0), remain(model.size());
//   vector<thread> threads;
//   for (int i=0; i!=num_threads; ++i) {
//     threads.emplace_back(thread(&ribo_model::fit_irate_block, this, std::ref(rprofile), offset, block_size));
//     offset += block_size;
//     if (model.size()-offset <= block_size) break;
//   }
//   if (model.size()-1 >= offset)
//     threads.emplace_back(thread(&ribo_model::fit_irate_block, this, std::ref(rprofile), offset, model.size()-offset));
//   for (auto& t: threads) { t.join(); }
// }

ostream& ribo_model::print_erates(ostream& out)
{
  for (auto it=_codon2id.begin(); it!=_codon2id.end(); ++it)
    out<<it->first<<":"<<get_erate(it->second)<<" ";
  return out;
}

void ribo_model::build_codon_rateID_map()
{
  array<char, 4> base{'A','C','G','T'};
  array<char, 4> codon;
  rid_t i(0);
  codon[3] = '\0';
  for (auto b0: base) {
    codon[0] = b0;
    for (auto b1: base) {
      codon[1] = b1;
      for (auto b2: base) {
        codon[2] = b2;
        string c(codon.data());
        // ignore stop codons
        if (c=="TAA" or c=="TAG" or c=="TGA") continue;
        _codon2id[c] = i;
        //cout<<c<<" "<<i<<endl;
	++i;
      }
    }
  }
}

void ribo_model::build_transcript_rate_id_vec(const char* tfa_fname, const vector<rid_t>& refID, const transcript_info& tinfo)
{
  fasta_reader transcript_fa(tfa_fname);
  for (size_t i=0; i!=refID.size(); ++i) {
    vector<int> rate_id_vec(tseq2rid_vec(refID[i], transcript_fa, tinfo));
    state_type profile(rate_id_vec.size(), 0);
    //assert(rate_id_vec.size()>=3);
    model.emplace_back(tmodel{0,0,0,profile,0,rate_id_vec});
  }
}

vector<int> ribo_model::tseq2rid_vec(unsigned refid, const fasta_reader& transcript_fa, const transcript_info& tinfo)
{
  int start(tinfo.cds_start(refid)+3);
  int stop(tinfo.cds_stop(refid)-3);
  string tseq(transcript_fa.read_region(refid, start, stop));
  vector<int> rate_id_vec;
  for (size_t i=0; i<tseq.size(); i+=3) {
    string codon(tseq.substr(i,3));
    // stop codon in the middle fo the seq are mapped to 'W' according to
    // gencode translation fasta file
    if (codon=="TGA") codon="TGG";
    try {
      int rate_id = _codon2id.at(codon);
      rate_id_vec.push_back(rate_id);
    }
    catch (const std::out_of_range& e) {
      cerr<<"codon "<<codon<<" does not exist in codon map!\n";
      //JUST FOR TESTING PURPOSE HAVE TO DELETE LATER
      rate_id_vec.push_back(0);
    }
  }
  assert(rate_id_vec.size()==tinfo.cds_pep_len(refid)-2);
  return rate_id_vec;
}

vector<double> ribo_model::build_rate_vec(rid_t t)
{
  vector<double> rate_vec(codon_vec_len(t),0);
  // build elongation rate vector
  for (size_t i=0; i!=rate_vec.size(); ++i)
    rate_vec[i] = get_erate(model[t].codon_id_map[i]);
  // insert initiation rate at the beginning of rate vector
  rate_vec.insert(rate_vec.begin(), get_irate(t));
  return rate_vec;
}

// double compute_profile_simulator(const vector<double>& rate, state_type& profile)
// {
//   polysome tasep_tmp(&rate);
//   if (profile.size() != rate.size()-1) {
//     cout<<profile.size()<<" "<<rate.size()<<endl;
//     cout<<"profile length not equal to rate vector length-1!"<<endl;
//     return -1;
//   }
//   tasep_tmp.run();
//   profile = tasep_tmp.get_Aprob();
//   return tasep_tmp.compute_translation_rate();
// }

double compute_profile(const vector<double>& rate, state_type& profile)
{
  tasep tasep_tmp(rate);
  state_type& p = profile;
  if (rate.size() < 3) return -1;
  if (p.size() != rate.size()-1) {
    cout<<p.size()<<" "<<rate.size()<<endl;
    cout<<"profile length not equal to rate vector length-1!"<<endl;
    return -1;
  }
  double residule_norm(tasep_tmp.numeric_steady_state(p));
  return 0;
}
