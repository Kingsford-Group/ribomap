#include <ios>
#include <fstream>
#include <cstdlib>

#include "gencode_parser.hpp"
#include "footprint_assigner.hpp"
#include "bam_parser.hpp"
#include "math_utils.hpp"

//------function forward declarations-------//
bool translation_pipeline(const transcript_info& tinfo, const char* bam_fname, const char* tfasta_fname, const char* sf_fname, const char* log_fname, int offset);
bool bowtie_best_log(const transcript_info& tinfo, const char* bam_fname, const char* tfasta_fname, const char* sf_fname, const char* log_fname, int offset);

int main(int argc, char ** argv)
{
  if (argc != 8) {
    cout<< "Usage: ./ribomap input_dir bam_core transcript_fasta gtf_fn sailfish_result output_dir footprint_offset"<<endl;
    exit(1);
  }

  string idir(argv[1]);
  string file_core(argv[2]);
  string odir(argv[6]);
  string nodup_bam = idir + "/" + file_core + "_nodup.bam";
  string bb_bam = idir + "/" + file_core + "_best.bam";
  string ribomap_log = odir +"/" + file_core + "_ribomap.log";
  string bb_log = odir + "/" + file_core + "_bb.log";
  int offset = atoi(argv[7]);

  cout<<"getting transcript info...\n";
  transcript_info tinfo(argv[3], argv[4]);
  //bowtie_best_log(tinfo, bb_bam.c_str(), argv[3], argv[5], bb_log.c_str(),offset);
  translation_pipeline(tinfo, nodup_bam.c_str(), argv[3], argv[5], ribomap_log.c_str(),offset);
  return 0;
}

bool translation_pipeline(const transcript_info& tinfo, const char* bam_fname, const char* tfasta_fname, const char* sf_fname, const char* log_fname, int offset)
{
  //profile
  cout<<"constructing profile class...\n";
  ribo_profile rprofile(tinfo, sf_fname);
  cout<<"number of transcripts in profile class: "<<rprofile.number_of_transcripts()<<endl;
  //reads
  cout<<"getting read info...\n";
  fp_list_t fp_codon_list;
  expressed_read_codon_ranges_from_bam(fp_codon_list, bam_fname, tinfo, rprofile,offset);
  cout<<"number of reads within cds ranges: "<<fp_codon_list.size()<<endl;
  cout<<"initializing read counts..."<<endl;
  rprofile.initialize_read_count(fp_codon_list, false);
  cout<<"assigning reads..."<<endl;
  rprofile.assign_reads();
  cout<<"update count profile..."<<endl;
  rprofile.update_count_profile();
  cout<<"write results to log..."<<endl;
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  ofstream logfile(log_fname);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  for (size_t t=0; t!=rprofile.number_of_transcripts(); ++t) {
    logfile<<"refID: "<<refID_vec[t]<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID_vec[t])<<endl;
    logfile<<"tabd: "<<rprofile.get_tot_abundance(t)<<endl;
    logfile<<"rprofile: ";
    print_vec(rprofile.get_read_assignments(t),logfile)<<endl;
  }
  logfile.close();
  return false;
}

bool bowtie_best_log(const transcript_info& tinfo, const char* bam_fname, const char* tfasta_fname, const char* sf_fname, const char* log_fname, int offset)
{
  //profile
  cout<<"constructing profile class...\n";
  ribo_profile rprofile(tinfo, sf_fname);
  cout<<"number of transcripts in profile class: "<<rprofile.number_of_transcripts()<<endl;
  cout<<"getting read info...\n"; 
  fp_list_t fp_codon_list;
  expressed_read_codon_ranges_from_bam(fp_codon_list, bam_fname, tinfo, rprofile, offset);
  cout<<"number of read within cds ranges: "<<fp_codon_list.size()<<endl;
  cout<<"getting ribosome profile from bowtie best"<<endl;
  rprofile.single_map_read_count(fp_codon_list);
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  ofstream logfile(log_fname);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  for (size_t t=0; t!=rprofile.number_of_transcripts(); ++t) {
    logfile<<"refID: "<<refID_vec[t]<<"\ntid: "<<tinfo.get_tid(refID_vec[t])<<endl;
    print_vec(rprofile.get_read_assignments(t), logfile)<<"\n";
  }
  logfile.close();
  return false;
}
