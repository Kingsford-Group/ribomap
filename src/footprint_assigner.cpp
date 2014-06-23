#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <functional>
#include <iostream>

#include <dlib/optimization.h>

#include "gencode_parser.hpp"
#include "footprint_assigner.hpp"
#include "bam_parser.hpp"
#include "tasep_ode.hpp"
#include "math_utils.hpp"

/***
 * class ribo_profile: 
 * store the profile from the model, the expected read assignment, and perform the EM
 */
ribo_profile::ribo_profile(const transcript_info& tinfo, const char* sf_fname, double abundance_cutoff): nonzero_abundance_vec(boost::dynamic_bitset<>(tinfo.total_count())), total_read_count(0)
{
  rid_t pid(0);
  ifstream ifile(sf_fname);
  while(ifile.peek() == '#'){
    string line;
    getline(ifile,line);
  }
  string header;
  int len;
  double tpk, rpkm, kpkm, enk, enr;
  double total_abundance(0);
  while(ifile >> header >> len >> tpk >> rpkm >> kpkm >> enk >> enr){
    size_t i = header.find("|");
    if (i==string::npos){
      cerr<<"invalid transcript header!\n";
      exit(1);
    }
    if (tpk < abundance_cutoff) continue;
    string tid(header.substr(0,i));
    rid_t rid(tinfo.get_refID(tid));
    int plen(tinfo.cds_pep_len(rid));
    // total number of (t,i) fragments: #transcript x plen
    total_abundance += tpk * plen;

    // initialize profile list
    vector<double> count(plen,0);
    profile.emplace_back(tprofile{0, count, tpk});
    refID2pID[rid] = pid++;
    include_abundant_transcript(rid);
    //if (profile.size()==40) break;
  }
  // normalize abundance
  for (size_t t = 0; t!=profile.size(); ++t)
    profile[t].tot_abundance /= total_abundance;
}

bool ribo_profile::initialize_read_count(const fp_list_t&  fp_codon_list, bool normalize)
{
  for (auto r: fp_codon_list) {
    vector<read_loc_count> loc_count_list;
    double tot_count(r.count);
    for (auto pos: r.al_loci) {
      rid_t refID(pos.refID); //transcript index
      rid_t t(get_transcript_index(refID));
      int start(pos.start), stop(pos.stop);
      if (start<0 or stop>len(t)) {
	cout<<"profile index out of bound! readID:"<<*r.seqs.begin()<<" ";
	cout<<r.al_loci.size()<<" "<<refID<<" "<<start<<"-"<<stop<<" "<<len(t)<<endl;
	continue;
      }
      for (rid_t i=start; i!=stop; ++i)
	loc_count_list.emplace_back(read_loc_count{t,i,tot_count});
    }// for codon range list
    read_count_list.emplace_back(read_count{tot_count,loc_count_list});
    total_read_count += tot_count;
  }// for r
  if (normalize)
    for (auto& r: read_count_list)
      r.tot_count /= total_read_count;
  return false;
}

vector<rid_t> ribo_profile::get_expressed_transcript_ids() const
{
  vector<rid_t> refID_vec(refID2pID.size());
  for (auto p: refID2pID)
    refID_vec[p.second] = p.first;
  return refID_vec;
}

// bowtie-best: reads only map to one place, no need to build read_count_list
bool ribo_profile::single_map_read_count(const fp_list_t& fp_codon_list)
{
  for (auto r: fp_codon_list) {
    auto& pos(r.al_loci[0]);
    rid_t refID(pos.refID); //transcript index in ribo_profile.profile
    rid_t t(get_transcript_index(refID));
    // sanity check for transcript index
    if (t>number_of_transcripts()-1){
      cout<<"transcript index out of bound! readID:"<<*r.seqs.begin()<<" ";
      cout<<t<<" "<<number_of_transcripts()<<endl;
    }
    int start(pos.start), stop(pos.stop);
    // sanity check for profile index
    if (start<0 or stop>len(t)) {
      cout<<"profile index out of bound! readID:"<<*r.seqs.begin()<<" ";
      cout<<r.al_loci.size()<<" "<<refID<<" "<<start<<"-"<<stop<<" "<<len(t)<<endl;
    }
    for (rid_t i=start; i!=stop; ++i){
      add_read_count(t,i,r.count);
      add_tot_count(t,r.count);
    }
  }// for r
  return false;
}

// assign reads based on transcript abundance only 
bool ribo_profile::assign_reads()
{
  for (auto& r: read_count_list) {
    double current_read_abd(0);
    // round 1: get expected values
    // numerators stored in read_local_count.count
    // denominators stored in current_read_abd
    for (auto& p: r.loc_count_list) {
      p.count = get_tot_abundance(p.tid);
      current_read_abd += p.count;
    }//round 1
    if (current_read_abd == 0) continue;
    // round 2: update counts
    for (auto& p: r.loc_count_list)
      p.count *= r.tot_count/current_read_abd;
  }// for r
  return false;
}

// assign reads based on transcript abundance x location abundance
bool ribo_profile::assign_reads_with_model(const ribo_model& rmodel)
{
  for (auto& r: read_count_list) {
    double current_read_abd(0);
    // round 1: get expected values
    // numerators stored in read_local_count.count
    // denominators stored in current_read_abd
    for (auto& p: r.loc_count_list) {
      p.count = get_tot_abundance(p.tid) * rmodel.get_prob(p.tid,p.loc);
      current_read_abd += p.count;
    }//round 1
    if (current_read_abd == 0) continue;
    // round 2: update counts
    for (auto& p: r.loc_count_list)
      p.count *= r.tot_count/current_read_abd;
  }// for r
  return false;
}

void ribo_profile::reset_read_count()
{
  for (auto& t: profile) {
    t.tot_count = 0;
    fill_n(t.count.begin(), t.count.size(), 0);
  }
}

void ribo_profile::update_count_profile()
{
  reset_read_count();
  for (auto r: read_count_list) {
    for (auto p: r.loc_count_list) {
      add_read_count(p.tid,p.loc,p.count);
      add_tot_count(p.tid,p.count);
    }
  }
}

/***
 * class rate_model:
 * store model parameters, model computation for each transcript
 */
ribo_model::ribo_model(const char* tfa_fname, const ribo_profile& profiler, const transcript_info& tinfo): erate(vector<double>(61,0)), refID2pID(profiler.refID2pID)
{
  vector<rid_t> refID(profiler.get_expressed_transcript_ids());
  build_codon_rateID_map();
  build_transcript_rate_id_vec(tfa_fname, refID, tinfo);
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
    t = refID2pID.at(refID);
    set_irate(t, rate);
  }
  ifile.close();
}

void ribo_model::elongation_rate_from_file(const char* erate_fname)
{
  ifstream ifile(erate_fname);
  string codon;
  double rate;
  while (ifile >> codon >> rate) {
    int eid = codon2id[codon];
    set_erate(eid, rate);
  }
  ifile.close();
}

void ribo_model::update_profile(const ribo_profile& rprofile, bool ignore_no_map)
{
  for (size_t t=0; t!=rprofile.number_of_transcripts(); ++t) {
    if (ignore_no_map and rprofile.get_tot_count(t) == 0) continue;
    tmodel& m = get_model_ref(t);
    if (get_obj(t)<0) continue;
    vector<double> rate_vec(build_rate_vec(t));
    compute_profile(rate_vec, m.profile);
    compute_translation_rate(t, rate_vec);
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

double ribo_model::fit_single_irate(rid_t t, const vector<double>& count_vec, bool normalized)
{
  //cout<<"fit_single_irate: "<<t<<" "<<std::this_thread::get_id()<<endl;
  tmodel& m = get_model_ref(t);
  double& start_point = m.irate;
  vector<double> r(build_rate_vec(t));
  if (r.size() != m.profile.size() + 1) {
    cout<<"in ribo_model::fit_irate, rate_vec.size != profile.size + 1 !"<<endl;
    return 0;
  }
  vector<double> c(count_vec);
  if (normalized) normalize_vec(c);
  irate_finder optimizer(c, r, m);
  double obj_end(0);
  try {
    obj_end = dlib::find_min_single_variable(optimizer, start_point, 1e-200, 1e200, 1e-3, 1000); 
  }
  catch (dlib::optimize_single_variable_failure& e) {
    if (e.info == "timeout")
      set_obj(t,-2);
    else 
      set_obj(t,-1);
    cout<<" b"<<t<<"b "<<flush;
    //cout<<"optimization for transcript "<<t<<" not converged!"<<endl;
    return -1;
  }
  set_obj(t,obj_end);
  //cout<<t<<" "<<flush;
  //cout<<"initiation rate: "<<start_point<<" objective function value: "<<obj_end<<endl;
  return obj_end;
}

void ribo_model::fit_irate_block(const ribo_profile& rprofile, rid_t t_offset, int num_rate)
{
  // cout<<t_offset<<endl;
  //cout<<"fit_irate_block: "<<t_offset<<" "<<std::this_thread::get_id()<<endl;
  for (size_t t=t_offset; t!=t_offset+num_rate; ++t) {
    // if read assignments not abundant enough, ignore
    if (rprofile.get_tot_count(t) < rprofile.len(t) ) {
      set_obj(t,0);
      continue;
    }
    fit_single_irate(t,rprofile.get_read_assignments(t),true);
  }
}

void ribo_model::fit_irates_in_parallel(const ribo_profile& rprofile, int num_threads)
{
  int block_size = max(int(model.size()/double(num_threads)),1);
  int offset(0), remain(model.size());
  vector<thread> threads;
  for (int i=0; i!=num_threads; ++i) {
    threads.emplace_back(thread(&ribo_model::fit_irate_block, this, std::ref(rprofile), offset, block_size));
    offset += block_size;
    if (model.size()-offset <= block_size) break;
  }
  if (model.size()-1 >= offset)
    threads.emplace_back(thread(&ribo_model::fit_irate_block, this, std::ref(rprofile), offset, model.size()-offset));
  for (auto& t: threads) { t.join(); }
}

ostream& ribo_model::print_erates(ostream& out)
{
  for (auto it=codon2id.begin(); it!=codon2id.end(); ++it)
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
        codon2id[c] = i;
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
    assert(rate_id_vec.size()>=3);
    model.emplace_back(tmodel{0,0,profile,0,rate_id_vec});
  }
}

vector<int> ribo_model::tseq2rid_vec(unsigned refid, const fasta_reader& transcript_fa, const transcript_info& tinfo)
{
  int start(tinfo.cds_start(refid));
  int stop(tinfo.cds_stop(refid));
  string tseq(transcript_fa.read_region(refid, start, stop));
  vector<int> rate_id_vec;
  for (size_t i=0; i<tseq.size(); i+=3) {
    string codon(tseq.substr(i,3));
    // stop codon in the middle fo the seq are mapped to 'W' according to
    // gencode translation fasta file
    if (codon=="TGA") codon="TGG";
    try {
      int rate_id = codon2id.at(codon);
      rate_id_vec.push_back(rate_id);
    }
    catch (const std::out_of_range& e) {
      cerr<<"codon "<<codon<<" does not exist in codon map!\n";
      //JUST FOR TESTING PURPOSE HAVE TO DELETE LATER
      rate_id_vec.push_back(0);
    }
  }
  assert(rate_id_vec.size()==tinfo.cds_pep_len(refid));
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

/***
 * class irate_finder:
 * functor of the optimization package for finding a local optimal initiation rate
 * that makes the model profile and the assignment profile the most similar
 */
double irate_finder::operator() (double irate) 
{
  vector<double>& rate_vec(rate_vec_);
  rate_vec[0] = irate;
  state_type& p(model_.profile);
  compute_profile(rate_vec,p);
  state_type p_tmp(p);
  eps_to_zeros(p_tmp);
  if (any_zeros(p_tmp)) return 1e200;
  normalize_vec(p_tmp);
  // sample time regularly
  // quit if sovling time too long
  //if (counter_%10 == 0) {
    tnow_ = chrono::system_clock::now();
    tperiod_ = tnow_-tstart_;
    tinterval_ = tperiod_.count();
    if ((counter_ == 0 and tinterval_ > 4) or tinterval_ > 400)
      throw dlib::optimize_single_variable_failure("timeout");
    //cout<<std::this_thread::get_id()<<" "<<tinterval_<<" "<<counter_<<" | "<<flush;
  //}
  ++counter_;
  return compute_objective(p_tmp);
}

double irate_finder::compute_objective(const state_type& p) const
{
  //return euclidean_dist(p,count_vec_);
  return kl_divergence(p);
  //return -multinomial_likelihood(p);
}

double irate_finder::multinomial_likelihood(const state_type& p) const
{
  double loglikelihood(0);
  assert( p.size() == count_vec_.size());
  for (size_t i=0; i!=p.size(); ++i) {
    loglikelihood += count_vec_[i] * log(p[i]);
    if (std::isnan(loglikelihood)) {
      cout<<"log likelihood became nan! count: "<<count_vec_[i]<<" profile: "<<p[i]<<endl;
      break;
    }
  }
  return loglikelihood;
}

double irate_finder::kl_divergence(const state_type& p) const
{
  double divergence(0);
  assert( p.size() == count_vec_.size());
  for (size_t i=0; i!=p.size(); ++i) {
    if (count_vec_[i] == 0) continue;
    divergence += count_vec_[i]*log(count_vec_[i]/p[i]);
    if (std::isnan(divergence)) {
      cout<<"kl divergence became nan! count: "<<count_vec_[i]<<" profile: "<<p[i]<<endl;
      break;
    }
  }
  return divergence;
}

bool compute_profile(const vector<double>& rate, state_type& profile)
{
  tasep tasep_tmp(rate);
  state_type& p = profile;
  if (p.size() != rate.size()-1) {
    cout<<p.size()<<" "<<rate.size()<<endl;
    cout<<"profile length not equal to rate vector length-1!"<<endl;
    return true;
  }
  double residule_norm(tasep_tmp.numeric_steady_state(p));
  return false;
}

/* differences between rate vectors */
double vec_dist(const column_vec& vec_a, const column_vec& vec_b)
{
  double sqr_sum = 0;
  assert(vec_a.size()==vec_b.size());
  for (size_t i=0; i!=vec_a.size(); ++i)
    sqr_sum += pow(vec_a(i)-vec_b(i),2);
  return sqrt(sqr_sum);
}

void eps_to_zeros(vector<double>& vec, double eps)
{
  for (auto& v: vec) 
    if (abs(v)<eps) 
      v = 0;
}

