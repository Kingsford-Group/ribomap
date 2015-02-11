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



#include <ios>
#include <fstream>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <string>

#include "ezOptionParser.hpp"
#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "bam_parser.hpp"
#include "utils.hpp"
#include "abundance_rank.hpp"
#include "profile_writer.hpp"

const int MAX_READ_LEN = 200;
//------function forward declarations-------//
void usage(ez::ezOptionParser& opt);
bool readmapper_pipeline(const transcript_info& tinfo, const char* mRNA_bam, const char* ribo_bam, int lmin, int lmax,bool useSecondary, bool useRC, const string& offset, const char* tabd_fname, const string& tabd_type, double tabd_cutoff, const string& log_prefix);

int main(int argc, const char ** argv)
{
  //command line parser
  ez::ezOptionParser opt;
  opt.overview = "assign reads to transcript locations based on the estimated transcript abundance, output a per-transcript profile file.";
  opt.syntax = "riboprof [options]";
  opt.example = "riboprof --fasta transcript_fasta --cds_range cds_range_file --mrnabam mRNA_bam --ribobam ribo_bam --min_fplen 25 --max_fplen 33 --offset footprint_offset --sf salmon_result --tabd_cutoff 0.01  --out output_prefix\n";
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
  opt.add("", 1, 1, 0, "minimum read length for keeping a read", "-lmin", "--min_fplen");
  opt.add("", 1, 1, 0, "maximum read length for keeping a read", "-lmax", "--max_fplen");
  opt.add("", 1, 1, 0, "footprint P-site offset", "-p", "--offset");
  opt.add("", 1, 1, 0, "transcript fasta", "-f", "--fasta");
  opt.add("", 0, 1, 0, "cds range file", "-cds", "--cds_range");
  opt.add("", 0, 1, 0, "sailfish result", "-s", "--sf");
  opt.add("", 0, 1, 0, "cufflinks result", "-c", "--cl");
  opt.add("", 0, 1, 0, "express result", "-e", "--ep");
  opt.add("", 1, 1, 0, "transcript abundance cutoff", "-threshold", "--tabd_cutoff");
  opt.add("", 1, 1, 0, "output profile", "-o", "--out");
  opt.add("", 0, 0, 0, "use secondary alignments", "-sec", "--useSecondary");
  opt.add("", 0, 0, 0, "use alignments with RC flag", "-rc", "--useRC");
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

  string transcript_fa, cds_range(""), mRNA_bam, ribo_bam, offset, tabd_fname, tabd_type, oprefix;
  opt.get("-f")->getString(transcript_fa);
  if (opt.isSet("-cds"))
    opt.get("-cds")->getString(cds_range);
  else
    cerr<<"cds_range file not provided, assume transcript fasta only contains cds sequences.\n";
  opt.get("-m")->getString(mRNA_bam);
  opt.get("-r")->getString(ribo_bam);
  opt.get("-p")->getString(offset);
  // get at least one transcript abundance file name
  if (opt.isSet("-s")) {
    opt.get("-s")->getString(tabd_fname);
    tabd_type = "sailfish";
  }
  else if (opt.isSet("-c")) {
    opt.get("-c")->getString(tabd_fname);
    tabd_type = "cufflinks";
  }
  else if (opt.isSet("-e")) {
    opt.get("-e")->getString(tabd_fname);
    tabd_type = "express";
  }
  else {
    cerr<<"ERROR: no transcript abundance file provided!\n";
    return 1;
  }
  opt.get("-o")->getString(oprefix);
  int lmin(0), lmax(0);
  opt.get("-lmin")->getInt(lmin);
  opt.get("-lmax")->getInt(lmax);
  double tabd_cutoff(0);
  opt.get("-threshold")->getDouble(tabd_cutoff);
  bool useSecondary = opt.isSet("-sec");
  bool useRC = opt.isSet("-rc");
  cout<<"getting transcript info...\n";
  transcript_info tinfo(transcript_fa.c_str(), cds_range.c_str());
  cout<<"total number of transcripts in transcriptome: "<<tinfo.total_count()<<endl;
  readmapper_pipeline(tinfo, mRNA_bam.c_str(), ribo_bam.c_str(), lmin, lmax, useSecondary, useRC, offset, tabd_fname.c_str(), tabd_type, tabd_cutoff, oprefix);
  return 0;
}

void usage(ez::ezOptionParser& opt)
{
  std::string usage;
  opt.getUsage(usage);
  std::cout << usage;
}

bool readmapper_pipeline(const transcript_info& tinfo, const char* mRNA_bam, const char* ribo_bam, int lmin, int lmax, bool useSecondary, bool useRC, const string& offset, const char* tabd_fname, const string& tabd_type, double tabd_cutoff, const string& log_prefix)
{
  cout<<"assigning ribo-seq reads..."<<endl;
  cout<<"constructing profile class..."<<endl;
  ribo_profile rprofile(tinfo, tabd_fname, tabd_type, tabd_cutoff);
  cout<<"number of transcripts in profile class: "<<rprofile.number_of_transcripts()<<endl;
  cout<<"loading reads from bam..."<<endl;
  fp_list_t fp_rec;
  try {
    int psite_offset = stoi(offset);
    expressed_read_bases_from_bam(fp_rec, ribo_bam, tinfo, rprofile, "-", lmin, lmax, useSecondary, false, psite_offset);
  }
  catch (const std::invalid_argument& ia) {
    expressed_read_bases_from_bam(fp_rec, ribo_bam, tinfo, rprofile, "-", lmin, lmax, useSecondary, false, offset.c_str());
  }
  cout.precision(10);
  cout<<"assigning reads to frame 0 loci..."<<endl;
  if (not rprofile.assign_reads(fp_rec, unordered_set<int>{FRAME0}))
    cout<<"total input count not matching total assigning count!"<<endl;
  cout<<"assigning reads to frame 1 and 2 loci..."<<endl;
  if (not rprofile.assign_reads(fp_rec, unordered_set<int>{FRAME1, FRAME2}))
    cout<<"total input count not matching total assigning count!"<<endl;
  cout<<"assigning reads to UTR loci..."<<endl;
  if (not rprofile.assign_reads(fp_rec, unordered_set<int>{UTR5, UTR3}))
    cout<<"total input count not matching total assigning count!"<<endl;
  cout<<"assigning RNA-seq reads..."<<endl;
  ribo_profile mprofile(tinfo, tabd_fname, tabd_type, tabd_cutoff);
  cout<<"number of transcripts in profile class: "<<mprofile.number_of_transcripts()<<endl;
  cout<<"loading reads from bam..."<<endl;
  fp_rec.clear();
  expressed_read_bases_from_bam(fp_rec, mRNA_bam, tinfo, mprofile, "-", lmin, MAX_READ_LEN, true, useRC, -1);
  cout<<"assigning reads..."<<endl;
  if (not mprofile.assign_reads(fp_rec,unordered_set<int>{UTR5, FRAME0, FRAME1, FRAME2, UTR3}))
    cout<<"total input count not matching total assigning count!"<<endl;
  cout<<"write profile results..."<<endl;
  generate_profile_report(log_prefix, rprofile, mprofile, tinfo);
  cout<<"rank transcripts..."<<endl;
  abundance_rank rank(rprofile, tinfo);
  cout<<"transcripts in rank list: "<<rank.size()<<endl;
  rank.get_rank(100);
  rank.write_diff_list(log_prefix, 10);
  return false;
}
