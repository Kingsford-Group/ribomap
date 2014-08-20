#include <ios>
#include <fstream>
#include <cstdlib>
#include <numeric>

#include "ezOptionParser.hpp"
#include "gencode_parser.hpp"
#include "ribomap_profiler.hpp"
#include "bam_parser.hpp"
#include "abundance_rank.hpp"

//------function forward declarations-------//
void usage(ez::ezOptionParser& opt);
bool translation_pipeline(const transcript_info& tinfo, const char* bam_fname, const char* tfasta_fname, const char* fname, const char* log_fname, const string& filetype, int offset);

int main(int argc, const char ** argv)
{
  // command line parser
  ez::ezOptionParser opt;
  opt.overview = "assign reads to transcript locations based on the estimated transcript abundance, output a per-transcript profile file.";
  opt.syntax = "mapper [options]";
  opt.example = "mapper --bam input_bam --fasta transcript_fasta --gtf gtf_fn --sf sailfish_result --out output_profile --offset footprint_offset\n";
  opt.footer = "";

  opt.add( "", // Default.
           0, // Required?
           1, // Number of args expected.
           0, // Delimiter if expecting multiple args.
           "Print help/usage message", // Help description.
           "-h",     // Flag token.
           "-help", // Flag token.
           "--help", // Flag token.
           "--usage" // Flag token.
           );

  opt.add("", 1, 1, 0, "input bam", "-b", "--bam");
  opt.add("", 1, 1, 0, "transcript fasta", "-f", "--fasta");
  opt.add("", 1, 1, 0, "transcript gtf", "-g", "--gtf");
  opt.add("", 0, 1, 0, "sailfish result", "-s", "--sf");
  opt.add("", 0, 1, 0, "cufflinks result", "-c", "--cl");
  opt.add("", 0, 1, 0, "express result", "-e", "--ep");
  opt.add("", 1, 1, 0, "output profile", "-o", "--out");
  opt.add("", 1, 1, 0, "footprint P-site offset", "-p", "--offset");

  opt.parse(argc, argv);

  if (opt.isSet("-h")) {
    usage(opt);
    return 1;
  }

  vector<string> badOptions;
  uint32_t i;
  if(!opt.gotRequired(badOptions)) {
    for(i=0; i < badOptions.size(); ++i) {
      cerr << "ERROR: Missing required option " << badOptions[i] << ".\n";
    }
    usage(opt);
    return 1;
  }

  string bam_fname, transcript_fa, transcript_gtf, ifname, ofname, filetype;
  opt.get("-b")->getString(bam_fname);
  opt.get("-f")->getString(transcript_fa);
  opt.get("-g")->getString(transcript_gtf);

  // get at least one transcript abundance file name
  if (opt.isSet("-s")) {
    opt.get("-s")->getString(ifname);
    filetype = "sailfish";
  }
  else if (opt.isSet("-c")) {
    opt.get("-c")->getString(ifname);
    filetype = "cufflinks";
  }
  else if (opt.isSet("-e")) {
    opt.get("-e")->getString(ifname);
    filetype = "express";
  }
  else {
    cerr<<"ERROR: no transcript abundance file provided!\n";
    return 1;
  }
    
  opt.get("-o")->getString(ofname);
  int offset;
  opt.get("-p")->getInt(offset);

  cout<<"getting transcript info...\n";
  transcript_info tinfo(transcript_fa.c_str(), transcript_gtf.c_str());
  translation_pipeline(tinfo, bam_fname.c_str(), transcript_fa.c_str(), ifname.c_str(), ofname.c_str(), filetype, offset);
  return 0;
}

void usage(ez::ezOptionParser& opt)
{
  std::string usage;
  opt.getUsage(usage);
  std::cout << usage;
}

bool translation_pipeline(const transcript_info& tinfo, const char* bam_fname, const char* tfasta_fname, const char* ifname, const char* log_fname, const string& filetype, int offset)
{
  //profile
  cout<<"constructing profile class...\n";
  ribo_profile rprofile(tinfo, ifname, filetype, 1);
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
  cout<<"rank transcripts..."<<endl;
  abundance_rank rank(rprofile, tinfo);
  rank.get_rank(100);
  string fn(log_fname);
  size_t i(fn.find_last_of("."));
  string fn_core(fn.substr(0,i));
  rank.write_diff_list(fn_core, 10);
  cout<<"write profile results..."<<endl;
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  ofstream logfile(log_fname);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  for (size_t t=0; t!=rprofile.number_of_transcripts(); ++t) {
    auto& p = rprofile.get_read_assignments(t);
    double rabd = accumulate(p.begin(),p.end(),double(0));
    double tabd = rprofile.get_tot_abundance(t);
    logfile<<"refID: "<<refID_vec[t]<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID_vec[t])<<endl;
    logfile<<"rabd: "<<rabd<<endl;
    logfile<<"tabd: "<<tabd<<endl;
    logfile<<"te: "<<rabd/tabd<<endl;
    logfile<<"rprofile: ";
    print_vec(p,logfile)<<endl;
  }
  logfile.close();
  return false;
}
