#ifndef FOOTPRINT_ASSIGNER_HPP
#define FOOTPRINT_ASSIGNER_HPP

#include <unordered_map>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <chrono>

#include <boost/dynamic_bitset.hpp>
#include <dlib/matrix.h>

using namespace std;

//------class forward declarations------//
class transcript_info;
class fasta_reader;

class ribo_profile;
class ribo_model;
class irate_finder;

//------copy aliases from other headers to here------//
struct fp_record;
using fp_list_t = vector<fp_record>;
using state_type = vector<double>;
using time_point = chrono::time_point<std::chrono::system_clock>;
using time_period = chrono::duration<double>;


//------aliases in this header------//
using column_vec = dlib::matrix<double,0,1>;
using rid_t = unsigned;

//------const------//
const double EPSILON = 1e-7;
const unsigned SEED = 619048235;

//------classes------//
//------profile info per transcript------//
struct tprofile{
  // expected read count on each position 
  double tot_count;
  vector<double> count;
  // transcript abundance from sailfish
  double tot_abundance;
};
//------expected count info per read per transcript-location pair------//
struct read_loc_count {
  rid_t tid;
  rid_t loc;
  double count;
};
//------count info per read------//
struct read_count { 
  double tot_count;
  vector<read_loc_count> loc_count_list;
};

/***
 * class ribo_profile:
 * store the profile from the model, the expected read assignment, and perform the EM
 */
class ribo_profile{
public:
  vector<read_count> read_count_list; 
  unordered_map<rid_t, rid_t> refID2pID; //ref index from fasta --> profile index
  ribo_profile(const transcript_info& tinfo, const char* sf_fname, double abundance_cutoff = 1);
  bool initialize_read_count(const fp_list_t& fp_codon_list, bool normalize = true);
  vector<rid_t> get_expressed_transcript_ids() const;
  bool is_expressed(rid_t refID) const { return nonzero_abundance_vec[refID]; }
  const vector<double>& get_read_assignments(rid_t t) const { return profile[t].count;}
  double get_count(rid_t t, rid_t i) const { return profile[t].count[i]; }
  double get_tot_abundance(rid_t t) const { return profile[t].tot_abundance; }
  double get_tot_count(rid_t t) const { return profile[t].tot_count; }
  rid_t get_transcript_index(rid_t refID) const { return refID2pID.at(refID); }
  size_t len(rid_t t) const { return profile[t].count.size(); }
  size_t number_of_transcripts() const { return profile.size(); }
  bool single_map_read_count(const fp_list_t& fp_codon_list);
  bool assign_reads();
  bool assign_reads_with_model(const ribo_model& rmodel);
  void update_count_profile();
private:
  vector<tprofile> profile;
  // bit vector: 1--transcript has non-zero abundance; 0--zero abundance
  boost::dynamic_bitset<> nonzero_abundance_vec; 
  double total_read_count;
  void reset_read_count();
  void include_abundant_transcript(rid_t refID) { nonzero_abundance_vec.set(refID); }
  void add_read_count(rid_t t, rid_t i, double count) { profile.at(t).count.at(i) += count; }
  void add_tot_count(rid_t t, double count) { profile[t].tot_count += count; }
};

//------tasep model parameters and info per transcript------//
struct tmodel{
  double irate; //initiation rate
  double trate; //translation rate
  state_type profile; //steady-state probability vector of tasep model
  double obj;
  vector<int> codon_id_map; //codon id vector to convert to elongation rate seq
};
/***
 * class rate_model:
 * store model parameters, model computation for each transcript
 */
class ribo_model {
public:
  vector<tmodel> model;
  vector<double> erate;
  ribo_model(const char* tfa_fname, const ribo_profile& profiler, const transcript_info& tinfo);
  // model parameter initialization
  void set_irate(rid_t t, double rate) { model[t].irate = rate; }
  void random_initiation_rate(double imin, double imax, long seed=SEED);
  void uniform_initiation_rate(double irate);
  void initiation_rate_from_file(const char* irate_fname);
  void elongation_rate_from_file(const char* erate_fname);
  // initiation rate computation
  double fit_single_irate(rid_t t, const vector<double>& count_vec, bool normalized=false);
  void fit_irate_block(const ribo_profile& rprofile, rid_t t_offset, int num_rate=10);
  void fit_irates_in_parallel(const ribo_profile& rprofile, int num_threads);
  // model profile computation
  void update_profile(const ribo_profile& rprofile, bool ignore_no_map=true);
  void compute_translation_rate(rid_t t, const vector<double>& rate);
  // access class element
  tmodel& get_model_ref(rid_t t) {return model[t]; }
  double get_irate(rid_t t) const { return model[t].irate; }
  double get_trate(rid_t t) const { return model[t].trate; }
  double get_obj(rid_t t) const { return model[t].obj; }
  const state_type& get_profile(rid_t t) const { return model[t].profile; }
  double get_prob(rid_t t, rid_t i) const { return model[t].profile[i]; }
  vector<double> build_rate_vec(rid_t t);
  // output results
  ostream& print_erates(ostream& out=cout);
private:
  unordered_map<string,int> codon2id;
  const unordered_map<rid_t, rid_t>& refID2pID;
  void build_codon_rateID_map();
  void build_transcript_rate_id_vec(const char* tfa_fname, const vector<rid_t>& refID, const transcript_info& tinfo);
  vector<int> tseq2rid_vec(unsigned refid, const fasta_reader& transcript_fa, const transcript_info& tinfo);
  void set_trate(rid_t t, double rate) { model[t].trate = rate; }
  void set_obj(rid_t t, double obj) { model[t].obj = obj; }
  void set_erate(int eid, double rate) { erate[eid] = rate; }
  double get_erate(int eid) const { return erate[eid]; }
  size_t codon_vec_len(rid_t t) const { return model[t].codon_id_map.size(); }
};

/***
 * class irate_finder:
 * functor of the optimization package for finding a local optimal initiation rate
 * that makes the model profile and the assignment profile the most similar
 */
class irate_finder {
public:
  irate_finder(const vector<double>& count_vec, vector<double>& rate_vec, tmodel& model): count_vec_(count_vec), rate_vec_(rate_vec), model_(model), tstart_(chrono::system_clock::now()), tnow_(chrono::system_clock::now()), tinterval_(0), counter_(0) {}
  double operator() (double irate);
private:
  const vector<double>& count_vec_;
  vector<double>& rate_vec_;
  tmodel& model_;
  time_point tstart_, tnow_;
  time_period tperiod_;
  double tinterval_;
  int counter_;
  double compute_objective(const state_type& p) const;
  double multinomial_likelihood(const state_type& p) const;
  double kl_divergence(const state_type& p) const;
};

//------function declarations------//
bool compute_profile(const vector<double>& rate, state_type& profile);
template<class vector_class>
ostream& print_vec(vector_class v, ostream& out=cout)
{ 
  for(auto vi:v) out<<vi<<" "; 
  return out;
}

template<class vector_class>
bool any_zeros(vector_class vec)
{
  for (auto v: vec)
    if (v==0) return true;
  return false;
}

void eps_to_zeros(vector<double>& vec, double eps=1e-50);
double vec_dist(const column_vec& vec_a, const column_vec& vec_b);
#endif
