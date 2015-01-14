#ifndef TRANSCRIPT_MODEL_HPP
#define TRANSCRIPT_MODEL_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include <random>
#include <iostream>

#include "utils.hpp"

using namespace std;

//------class forward declarations------//
class transcript_info;
class ribo_profile;
class fasta_reader;

//------tasep model parameters and info per transcript------//
struct tmodel{
  double irate; //initiation rate
  double trate; //translation rate
  double tcnt;
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
  ribo_model(const char* tfa_fname, const vector<rid_t>& refIDs, const transcript_info& tinfo);
  // model parameter initialization
  void set_irate(rid_t t, double rate) { model[t].irate = rate; }
  void random_initiation_rate(double imin, double imax, long seed=SEED);
  void uniform_initiation_rate(double irate);
  void initiation_rate_to_file(const char* irate_fname);
  void initiation_rate_from_file(const char* irate_fname);
  void elongation_rate_from_file(const char* erate_fname);
  // initiation rate computation
  double fit_single_irate(rid_t t, const vector<double>& count_vec, bool normalized=false);
  void fit_irate_block(const ribo_profile& rprofile, rid_t t_offset, int num_rate=10);
  void fit_irates_in_parallel(const ribo_profile& rprofile, int num_threads);
  // model profile computation
  void update_profile(const ribo_profile& rprofile, bool ignore_no_map=true, int nproc=1);
  void compute_translation_rate(rid_t t, const vector<double>& rate);
  void compute_transcript_abundance(rid_t, const vector<double>& count);
  // access class element
  tmodel& get_model_ref(rid_t t) {return model[t]; }
  size_t number_of_transcripts() const { return model.size(); }
  size_t len(rid_t t) const { return model[t].profile.size(); }
  double get_irate(rid_t t) const { return model[t].irate; }
  double get_trate(rid_t t) const { return model[t].trate; }
  double get_tcnt(rid_t t) const { return model[t].tcnt; }
  double get_obj(rid_t t) const { return model[t].obj; }
  const state_type& get_profile(rid_t t) const { return model[t].profile; }
  double get_prob(rid_t t, rid_t i) const { return model[t].profile[i]; }
  vector<double> build_rate_vec(rid_t t);
  // output results
  ostream& print_erates(ostream& out=cout);
private:
  unordered_map<string,int> _codon2id;
  const unordered_map<rid_t, rid_t> _refID2pID;
  void build_codon_rateID_map();
  void build_transcript_rate_id_vec(const char* tfa_fname, const vector<rid_t>& refID, const transcript_info& tinfo);
  vector<int> tseq2rid_vec(unsigned refid, const fasta_reader& transcript_fa, const transcript_info& tinfo);
  void set_trate(rid_t t, double rate) { model[t].trate = rate; }
  void set_tcnt(rid_t t, double tcnt) { model[t].tcnt = tcnt; }
  void set_obj(rid_t t, double obj) { model[t].obj = obj; }
  void set_erate(int eid, double rate) { erate[eid] = rate; }
  double get_erate(int eid) const { return erate[eid]; }
  size_t codon_vec_len(rid_t t) const { return model[t].codon_id_map.size(); }
};


double compute_profile(const vector<double>& rate, state_type& profile);
#endif
