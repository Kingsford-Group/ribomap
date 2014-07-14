#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include "ezOptionParser.hpp"

using namespace std;
using seq_id_map_t = map<string, vector<string> >;

void usage(ez::ezOptionParser& opt);
void merge_dup_fq_to_fa(const string& ifname, const string& ofname, int trimlen);

int main(int argc, const char ** argv)
{
  // command line parser
  ez::ezOptionParser opt;
  opt.overview = "trim fastq seqs, merge duplicate seqs, output merged/unmerged fasta files";
  opt.syntax = "preprocess_fq [options]";
  opt.example = "preprocess_fq -i input.fq -od output_directory -l trimlen\n";
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

  opt.add("", 1, 1, 0, "input fastq", "-i", "--in");
  opt.add("", 1, 1, 0, "output fasta", "-o", "--out");
  opt.add("25", 0, 1, 0, "trim length", "-l", "--trimlen");

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

  string ifname, ofname;
  opt.get("-i")->getString(ifname);
  opt.get("-o")->getString(ofname);
  int trimlen;
  opt.get("-l")->getInt(trimlen);

  merge_dup_fq_to_fa(ifname, ofname, trimlen);
  cerr<<"\nDone."<<endl;
  return 0;
}

void usage(ez::ezOptionParser& opt)
{
  std::string usage;
  opt.getUsage(usage);
  std::cout << usage;
}

void merge_dup_fq_to_fa(const string& ifname, const string& ofname, int trimlen)
{
  cerr<<"start merging duplicate sequences..."<<endl;
  seq_id_map_t seq2ids;
  string line, header, seq;
  ifstream ifile(ifname);
  int ln(0);
  getline(ifile, line);
  while(ifile) {
    if (line.front() == '@') {
      header = line;
      getline(ifile, line);
      if (line.size()>1) {
	header[0] = '>';
	seq = line.substr(0,trimlen);
	seq2ids[seq].push_back(header);
      }
      getline(ifile, line);
      getline(ifile, line);
    }
    if (++ln % 1000 == 0) {
      std::cerr << "\r\rParsed " << ln << " reads";
    }
    getline(ifile, line);
  }
  ifile.close();
  cerr<<"writing seqs to fasta..."<<endl;
  ln = 0;
  ofstream ofile(ofname);
  for (auto p: seq2ids) {
    seq = p.first;
    auto header_orig = p.second[0];
    auto idx = header_orig.find(" ");
    header = header_orig.substr(0,idx)+string("_")+to_string(p.second.size());
    ofile<<header<<endl;
    ofile<<seq<<endl;
    if (++ln % 1000 == 0) {
      std::cerr << "\r\rwrote " << ln << " reads";
    }
  }
  ofile.close();
}
