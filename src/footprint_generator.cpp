/*
g++ -std=c++11 -I/home/hw1/.local/include -lz -lbz2 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -Wno-deprecated-declarations -Wno-virtual-move-assign -I/opt/local/include -L/opt/local/lib -lboost_system -lboost_filesystem -O3 -o footprint_generator footprint_generator.cpp footprint_assigner.cpp tasep_ode.cpp bam_parser.cpp gencode_parser.cpp utils.cpp

// TEST CASE TINY:
// reads from a real gene with 3 isoforms with synthetic transcript abundance
./footprint_generator ../../data/human_mouse_rp_Guo/fasta/ENSG00000145592.9.fq /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_transcripts_filter.fa /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.annotation.gtf gene_ENSG00000145592.9.txt /home/hw1/riboseq/results/outputs/single_gene/truth_ENSG00000145592.9.log

// TEST CASE FINAL:
// all expressed transcripts with tpm > 1 (35,410)
./footprint_generator ../data/human_mouse_rp_Guo/fasta/all_expressed.fq /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_transcripts_filter.fa /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.annotation.gtf /home/hw1/riboseq/data/human_mouse_rp_Guo/sailfish_quant/quant.sf /home/hw1/riboseq/results/outputs/footprint_generator_logs/footprint_all_expressed.log
*/

#include <ios>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <numeric>

#include "math_utils.hpp"
#include "gencode_parser.hpp"
#include "footprint_assigner.hpp" 

using namespace std;
using time_point = chrono::time_point<std::chrono::system_clock>;
using time_period = chrono::duration<double>;

int main(int argc, char ** argv)
{
  if (argc != 6) {
    cout<< "Usage: ./footprint_generator fastq_fname transcript_fasta gtf_fn sailfish_result log_fname"<<endl;
    exit(1);
  }

  // cout.setf(ios::scientific);
  // cout.precision(3);
  
  cout<<"getting transcript info...\n";
  transcript_info tinfo(argv[2], argv[3]);
  //profile
  cout<<"constructing profile class...\n";
  ribo_profile rprofile(tinfo, argv[4]);
  cout<<"number of transcripts in profile class: "<<rprofile.number_of_transcripts()<<endl;
  //rate
  ribo_model ribo_rate(argv[2], rprofile, tinfo);
  cout<<"tRNA abundance as elongation rate"<<endl;
  ribo_rate.elongation_rate_from_file("elongation_rate_human.txt");
  cout<<"initiation rate is randomly set to be 100 times smaller than elongation rate range"<<endl;
  ribo_rate.random_initiation_rate(0.03,0.3,SEED);
  // cout<<"uniform rate as true rate"<<endl;
  // ribo_rate.uniform_initializer(0.5,0.5);
  // log file to store both model profile and actual profile
  cout<<"writing synthetic log..."<<endl;
  ofstream logfile(argv[5]);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  int64_t read_count = 100000000;
  int read_len = 36;
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  logfile<<"number of reads: "<<read_count<<endl;
  logfile<<"read length: "<<read_len<<endl;
  logfile<<"erates: ";
  ribo_rate.print_erates(logfile)<<endl;
  cout<<"compute model profile..."<<endl;
  ribo_rate.update_profile(rprofile,false);
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

  cout<<"generating synthetic footprints..."<<endl;
  fasta_reader transcript_fa(argv[2]);
  vector<string> read_buffer;
  int64_t tot_count(0);
  int plen(0), refID(0), readID(0);
  double abundance(0);
  string quality(read_len,'a');
  // pipeline of generating synthetic footprints
  // for each transcript:
  for (size_t t = 0; t!= rprofile.number_of_transcripts(); ++t) {
    // step 1: get model profile
    refID = refID_vec[t];
    plen = rprofile.len(t);
    tot_count = round(read_count * rprofile.get_tot_abundance(t));
    auto& p = ribo_rate.get_profile(t);
    // step 2: compute count profile
    vector<int64_t> profile_count(plen,0);
    for (size_t i=0; i!=plen; ++i)
      profile_count[i] = round(p[i]* tot_count);
    // step 3: profile write to log
    logfile<<"refID: "<<refID<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID)<<endl;
    logfile<<"irate: "<<ribo_rate.get_irate(t)<<endl;
    logfile<<"trate: "<<ribo_rate.get_trate(t)<<endl;
    logfile<<"len: "<<plen<<endl;
    logfile<<"tcnt: "<<tot_count<<endl;
    logfile<<"tabd: "<<rprofile.get_tot_abundance(t)<<endl;
    logfile<<"mprofile: ";
    print_vec(p,logfile)<<endl;
    logfile<<"rprofile: ";
    print_vec(profile_count,logfile)<<endl;
    // step 4: generate reads
    for (size_t i=0; i!=plen; ++i) {
      if (profile_count[i] == 0) continue;
      // no need to skipping the first codon now
      unsigned mid_pos = i*3 + tinfo.cds_start(refID);
      // sanity check whether I got this
      if (mid_pos == 0) {
	cout<<"mid_pos cannot be zero"<<endl;
	continue;
      }
      unsigned start = max((unsigned)0,mid_pos-read_len/2);
      uint64_t stop = min((uint64_t)mid_pos+read_len/2, transcript_fa.length(refID));
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
  ofstream ofile(argv[1]);
  for (auto read: read_buffer)
    ofile<<read;
  ofile.close();
  return 0;
}
