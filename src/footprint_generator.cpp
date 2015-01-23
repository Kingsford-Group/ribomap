#include <string>
#include <ios>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <numeric>
#include <cstdint>

#include "utils.hpp"
#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "transcript_model.hpp" 

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 14) {
    cout<< "Usage: ./footprint_generator transcript_fa cds_range sailfish_result erate_fn ilow ihigh read_cnt read_len offset log_fname fq_fname nproc cutoff"<<endl;
    exit(1);
  }

  // cout.setf(ios::scientific);
  // cout.precision(3);
  
  // command parser
  char* ref_fa=argv[1];
  char* cds_range=argv[2]; 
  char* sf_fn=argv[3];
  char* erate_fn=argv[4];
  char* log_fn=argv[10];
  char* fq_fn=argv[11];
  double ilow(std::stof(argv[5])), ihigh(std::stof(argv[6])), cutoff(std::stof(argv[13]));
  int read_len(std::stoi(argv[8])), offset(std::stoi(argv[9])), nproc(std::stoi(argv[12]));
  uint_fast64_t read_count(std::stoi(argv[7]));
  //cds range
  cout<<"getting transcript info...\n";
  transcript_info tinfo(ref_fa, cds_range);
  //profile
  cout<<"constructing profile class...\n";
  ribo_profile rprofile(tinfo, sf_fn, "sailfish", cutoff);
  cout<<"number of transcripts in profile class: "<<rprofile.number_of_transcripts()<<endl;
  //model
  cout<<"initializing tasep model parameters..."<<endl;
  ribo_model ribo_rate(ref_fa, rprofile, tinfo);
  cout<<"tRNA abundance as elongation rate"<<endl;
  ribo_rate.elongation_rate_from_file(erate_fn);
  cout<<"initiation rate range: "<<ilow<<" "<<ihigh<<endl;
  ribo_rate.random_initiation_rate(ilow,ihigh,SEED);
  cout<<"compute model profile..."<<endl;
  time_point start, end;
  start = chrono::system_clock::now();
  ribo_rate.update_profile(rprofile,false,nproc);
  end = chrono::system_clock::now();
  time_period t = end-start;
  cout<<"solving the entire system uses: "<<t.count()<<" secs"<<endl;
  cout<<"writing synthetic log..."<<endl;
  //log
  ofstream logfile(log_fn);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  logfile<<"number of reads: "<<read_count<<endl;
  logfile<<"read length: "<<read_len<<endl;
  logfile<<"erates: ";
  ribo_rate.print_erates(logfile)<<endl;
  cout<<"getting total abundance of the transcriptome..."<<endl;
  double tot_abd(0);
  for (size_t t=0; t!= rprofile.number_of_transcripts(); ++t) {
    auto& p = ribo_rate.get_profile(t);
    double pabd(0);
    pabd = accumulate(p.begin(),p.end(), pabd);
    double tabd(rprofile.get_tot_abundance(t));
    tot_abd += tabd * pabd;
  }
  cout<<"total abundance is "<<tot_abd<<endl;
  cout<<"in order to generate "<<read_count<<" footprints ";
  read_count /= tot_abd;
  cout<<read_count<<" mRNA fragments are needed. "<<endl;
  int refID(0), readID(0);
  int64_t tot_count(0);
  double abundance(0);  
  fasta_reader transcript_fa(ref_fa);
  vector<string> read_buffer;
  string quality(read_len,'a');
  // pipeline of generating synthetic footprints
  // for each transcript:
  for (size_t t = 0; t!= rprofile.number_of_transcripts(); ++t) {
    // step 1: get model profile
    abundance = rprofile.get_tot_abundance(t);
    if (abundance==0) continue;
    refID = refID_vec[t];
    tot_count = round(read_count * abundance);
    auto& p = ribo_rate.get_profile(t);
    if (accumulate(p.begin(), p.end(), double(0))==0) continue;
    // step 2: compute count profile
    vector<uint_fast64_t> profile_count(p.size(),0);
    for (size_t i=0; i!=p.size(); ++i)
      profile_count[i] = round(p[i]* tot_count);
    if (accumulate(profile_count.begin(), profile_count.end(), double(0))==0) continue;
    // step 3: profile write to log
    logfile<<"refID: "<<refID<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID)<<endl;
    logfile<<"irate: "<<ribo_rate.get_irate(t)<<endl;
    logfile<<"trate: "<<ribo_rate.get_trate(t)<<endl;
    logfile<<"plen: "<<p.size()<<endl;
    logfile<<"tcnt: "<<tot_count<<endl;
    logfile<<"tabd: "<<rprofile.get_tot_abundance(t)*rprofile.len(t)<<endl;
    logfile<<"mprofile: "<<p<<endl;
    logfile<<"rprofile: "<<profile_count<<endl;
    // step 4: generate reads
    for (size_t i=0; i!=p.size(); ++i) {
      if (profile_count[i] == 0) continue;
      // first and last codons are skipped
      unsigned mid_pos = (i+1)*3 + tinfo.cds_start(refID);
      // sanity check whether I got this
      if (mid_pos == 0) {
    	cout<<"mid_pos cannot be zero"<<endl;
    	continue;
      }
      unsigned start = max((unsigned)0,mid_pos-offset);
      uint64_t stop = min((uint64_t)mid_pos+(read_len-offset), transcript_fa.length(refID));
      if (stop-start < read_len) continue;
      // fasta reader is zero-based, tinfo.cds_start has been adjusted for that
      string tseq(transcript_fa.read_region(refID, start, stop));
      string comment(tinfo.get_tid(refID)+" "+to_string(refID)+" codon: "+to_string(i)+" base: "+to_string(start)+"-"+to_string(stop));
      // generate reads
      for (size_t dup_cnt=0; dup_cnt!=profile_count[i]; ++dup_cnt) {
    	string read("@"+to_string(readID)+" "+comment+"\n"+
    		    tseq+"\n"+
    		    "+"+to_string(readID)+" "+comment+"\n"+
    		    quality+"\n");
    	read_buffer.push_back(read);
    	++readID;
      }// for dup_cnt
    }// for i
  }// for t
  logfile.close();
  cout<<"generated "<<read_buffer.size()<<" reads."<<endl;
  cout<<"reshuffle reads..."<<endl;
  random_shuffle(read_buffer.begin(), read_buffer.end());
  cout<<"write to fastq file..."<<endl;
  ofstream ofile(fq_fn);
  for (auto read: read_buffer)
    ofile<<read;
  ofile.close();
  return 0;
}
