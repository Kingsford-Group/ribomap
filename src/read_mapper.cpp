#include <ios>
#include <fstream>
#include <cstdlib>
#include <numeric>

#include "ezOptionParser.hpp"
#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "bam_parser.hpp"
#include "utils.hpp"
#include "abundance_rank.hpp"

//------function forward declarations-------//
void usage(ez::ezOptionParser& opt);
bool readmapper_pipeline(const transcript_info& tinfo, const char* mRNA_bam, const char* ribo_bam, const char* ifname, const char* log_fname, const string& filetype, int offset);

int main(int argc, const char ** argv)
{
  //command line parser
  ez::ezOptionParser opt;
  opt.overview = "assign reads to transcript locations based on the estimated transcript abundance, output a per-transcript profile file.";
  opt.syntax = "riboprof [options]";
  opt.example = "riboprof --mrnabam mRNA_bam --ribobam ribo_bam --fasta transcript_fasta --cds_range cds_range_file --sf salmon_result --offset footprint_offset --out output_profile\n";
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

  opt.add("", 1, 1, 0, "RNAseq alignment bam", "-m", "--mrnabam");
  opt.add("", 1, 1, 0, "riboseq alignment bam", "-r", "--ribobam");
  opt.add("", 1, 1, 0, "transcript fasta", "-f", "--fasta");
  opt.add("", 0, 1, 0, "cds range file", "-cds", "--cds_range");
  opt.add("", 0, 1, 0, "sailfish result", "-s", "--sf");
  opt.add("", 0, 1, 0, "cufflinks result", "-c", "--cl");
  opt.add("", 0, 1, 0, "express result", "-e", "--ep");
  opt.add("", 1, 1, 0, "footprint P-site offset", "-p", "--offset");
  opt.add("", 1, 1, 0, "output profile", "-o", "--out");

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

  string mRNA_bam, ribo_bam, transcript_fa, cds_range(""), ifname, ofname, filetype;
  opt.get("-m")->getString(mRNA_bam);
  opt.get("-r")->getString(ribo_bam);
  opt.get("-f")->getString(transcript_fa);
  if (opt.isSet("-cds"))
    opt.get("-cds")->getString(cds_range);
  else
    cerr<<"cds_range file not provided, assume transcript fasta only contains cds sequences.\n";

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
  transcript_info tinfo(transcript_fa.c_str(), cds_range.c_str());
  cout<<"total number of transcripts in transcriptome: "<<tinfo.total_count()<<endl;
  readmapper_pipeline(tinfo, mRNA_bam.c_str(), ribo_bam.c_str(), ifname.c_str(), ofname.c_str(), filetype, offset);
  return 0;
}

void usage(ez::ezOptionParser& opt)
{
  std::string usage;
  opt.getUsage(usage);
  std::cout << usage;
}


ribo_profile assign_reads(const transcript_info& tinfo, const char* bam_fname, const char* ifname, const string& filetype, int offset)
{
  cout<<"constructing profile class...\n";
  ribo_profile profiler(tinfo, ifname, filetype, 0.01);
  cout<<"number of transcripts in profile class: "<<profiler.number_of_transcripts()<<endl;
  //reads
  cout<<"getting read info...\n";
  fp_list_t fp_codon_list;
  expressed_read_codon_ranges_from_bam(fp_codon_list, bam_fname, tinfo, profiler, offset, "-");
  cout<<"initializing read counts..."<<endl;
  profiler.initialize_read_count(fp_codon_list, false);
  cout<<"assigning reads..."<<endl;
  profiler.assign_reads();
  cout<<"update count profile..."<<endl;
  profiler.update_count_profile();
  return profiler;
}

bool readmapper_pipeline(const transcript_info& tinfo, const char* mRNA_bam, const char* ribo_bam, const char* ifname, const char* log_fname, const string& filetype, int offset)
{
  //profile
  ribo_profile mprofile = assign_reads(tinfo, mRNA_bam, ifname, filetype, offset);
  ribo_profile rprofile = assign_reads(tinfo, ribo_bam, ifname, filetype, offset);

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
    vector<double> rp=rprofile.get_read_assignments(t);
    vector<double> mp=mprofile.get_read_assignments(t);
    // bias correction by normalizing ribo profile with mRNA profile
    vector<double> np(rp);
    for (size_t i=0; i!= np.size(); ++i) {
      if ( mp[i] != 0 ) 
	np[i] /= mp[i];
      else
	np[i] = 0;
    }
    // ribosome loads: sum of ribosome profile
    double rabd = accumulate(rp.begin(), rp.end(), double(0));
    double tabd = rprofile.get_tot_abundance(t);
    logfile<<"refID: "<<refID_vec[t]<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID_vec[t])<<endl;
    logfile<<"rabd: "<<rabd<<endl;
    logfile<<"tabd: "<<tabd<<endl;
    logfile<<"te: "<<rabd/tabd<<endl;
    logfile<<"ribo profile: "<<rp<<endl;
    logfile<<"normalized ribo profile: "<<np<<endl;
    logfile<<"mRNA profile: "<<mp<<endl;
  }
  logfile.close();
  return false;
}
