/*
g++ -std=c++11 -I/home/hw1/.local/include/ -lz -lbz2 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -Wno-deprecated-declarations -I/opt/local/include -L/opt/local/lib -lboost_system -lboost_filesystem -O3 -o gencode_parser gencode_parser.cpp utils.cpp

./gencode_parser /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_transcripts_filter.fa /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.pc_translations_filter.fa /data/iGenomes/Homo_sapiens/GenCode/gencode.v18.annotation.gtf
*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>
#include "gencode_parser.hpp"
#include "utils.hpp"


int test_case_gencode_parser(int argc, char ** argv)
{
  if (argc != 4) {
    cout << "Usage: gencode_parser transcript_fasta peptide_fasta gtf_fname" << endl;
    return 1;
  }
  transcript_info tinfo(argv[1], argv[3]);
  return 0;
  for (size_t i=0; i!=tinfo.total_count(); ++i) {
    cout<<tinfo.get_tid(i)<<" "<<tinfo.cds_start(i)<<"-"<<tinfo.cds_stop(i)<<endl;
    if (i>20) break;
  }
  
  // check whether encoding is right
  fasta_reader transcript_fa(argv[1]), peptide_fa(argv[2]);
  cout<<"check whether encoding is right\n";
  cout<<tinfo.total_count()<<endl;
  int diff_cnt = 0;
  for (size_t i=0; i!=tinfo.total_count(); ++i){
    int start = tinfo.cds_start(i), stop = tinfo.cds_stop(i);
    int plen = tinfo.cds_pep_len(i);
    if (plen < 3) continue;
    string pconvert, pseq(peptide_fa.read_seq(i)), tseq(transcript_fa.read_region(i,start,stop));
    encode_peptide(tseq,pconvert);
    if (tinfo.phase(i) != 0)
      pconvert = "X"+pconvert;
    if (pseq.compare(pconvert) != 0){
      cout<<pseq<<endl;
      cout<<pconvert<<endl;
      ++diff_cnt;
    }
    if (i>100) return 0;
  }
  cout<<diff_cnt<<endl;
  return 0;
}

transcript_info::transcript_info(const char* tfname, const char* gtf_fname, const char* cereal_name)
{
  // file exists, load in values
  if (fileExists(string(cereal_name))) {
    cout<<"read transcripts from object cache\n";
    ifstream is(cereal_name);
    cereal::BinaryInputArchive iarchive(is);
    iarchive(*this);
  }
  // file does not exist, save value to file
  else {
    cout<<"read in transcriptome\n";
    get_info_from_fasta(tfname);
    build_tid_idx_map();
    cout<<"read in phase\n";
    get_info_from_gtf(gtf_fname);
    cout<<"adjust cds ranges\n";
    adjust_cds_ranges();
    cout<<"cache does not exist, creat cereal object\n";
    ofstream os(cereal_name);
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(*this);
  }
}

/*
transcript.fa header:
0 transcript-id|
1 gene-id|
2 Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
3 Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
4 transcript-name|
5 gene-name|
6 sequence-length|
7 5'-UTR (3'-UTR if reverse strand) location in the transcript|
8 CDS location in the transcript|
9 3'-UTR (5'-UTR if reverse strand) location in the transcript
 */
bool transcript_info::get_info_from_fasta(const char* tfname)
{
  ifstream tfile(tfname);
  if (!tfile.good())
    return 1;
  seqan::RecordReader<ifstream, seqan::SinglePass<> > treader(tfile);
  seqan::CharString theader,tseq;
  // int i(0);
  while(!atEnd(treader)){
    if (readRecord(theader, tseq, treader, seqan::Fasta()) != 0){
      std::cerr<<"ERROR reading FASTA "<<tfname<<std::endl;
      return 1;
    }
    string theader_str(toCString(theader));
    vector<string> twords;
    string_split(theader_str, "|", twords);

    //fill in transcript record:
    string tid(twords[0]), gid(twords[1]), tname(twords[4]), gname(twords[5]);
    int start = 0, stop = 0;
    for (auto it = twords.begin()+7; it!=twords.end(); ++it){
      string tmp_wrd = *it;
      if (tmp_wrd.find("CDS:")==0){
	string_lstrip("CDS:",*it);
	vector<string> cds_range;
	string_split(*it, "-", cds_range);
	start = atoi(cds_range[0].c_str());
	stop = atoi(cds_range[1].c_str());
      }
    }//for twords
    tlist.push_back(tprop(tid,gid,tname,gname,start,stop));
    // cout<<tid<<" "<<gid<<" "<<tname<<" "<<gname<<" "<<start<<"-"<<stop<<endl;
    // if (++i>10) break;
  }
  return false;
}

void transcript_info::build_tid_idx_map()
{
  for(int i=0; i!=tlist.size(); ++i)
    tid2refid[tlist[i].tid] = i;
}

/*
gtf fields:
0: seqname
1: source
2: feature
3: start
4: end
5: score
6: strand
7: frame
8: attribute
      gene_id ENSGXXXXXXXXXXX *
      transcript_id ENSTXXXXXXXXXXX *
      gene_type list of biotypes
gene_status {KNOWN, NOVEL, PUTATIVE}
      gene_name string
      transcript_type list of biotypes
      transcript_status {KNOWN, NOVEL, PUTATIVE}
      transcript_name string
      level 1 (verified loci), 2 (manually annotated loci), 3 (automatically annotated loci)
*/
bool transcript_info::get_info_from_gtf(const char* gtf_fname)
{
  // Open input stream
  seqan::GffStream gffIn(gtf_fname, seqan::GffStream::READ, seqan::GffStream::GTF);
  if (!isGood(gffIn)){
    std::cerr << "ERROR: Could not open "<<gtf_fname<<"\n";
    return 1;
  }
  // Read the file record by record.
  seqan::GffRecord record;
  int i = 0;
  while (!atEnd(gffIn)){
    if (readRecord(record, gffIn) != 0){
      std::cerr << "ERROR: Problem reading from gtf file\n";
      return 1;
    }
    // only look at CDS entries
    if (record.type == "CDS"){
      // read in related fields
      string tid;
      int exon_num, phase;
      for (size_t i = 0; i<length(record.tagName); ++i){
	if (record.tagName[i] == "transcript_id")
	  tid = string(toCString(record.tagValue[i]));
	if (record.tagName[i] == "exon_number")
	  exon_num = atoi(toCString(record.tagValue[i]));
      }
      // check whether tid in transcript list
      auto it = tid2refid.find(tid);
      if (it != tid2refid.end()){
	int id = it->second;
	if (tlist[id].exon_num == MAX_EXON_NUM){
	  tlist[id].chrm = string(toCString(record.ref));
	  tlist[id].strand = (record.strand=='+') ? true : false;
	}
	// update phase
	if (tlist[id].exon_num > exon_num){
	  tlist[id].exon_num = exon_num;
	  tlist[id].phase = record.phase - '0';
	}// if exon_num update
	//if (id==10) break;
      }//if find key
    }// if CDS
  }// while not eof
  //for (auto t: tlist)
  //  cout<<t.tid<<" "<<t.phase<<endl;
  return false;
}

/*
peptide.fa header
0 transcript-id|
1 gene-id|
2 Havana-gene-id (if the gene contains manually annotated transcri\
pts, '-' otherwise)|
3 Havana-transcript-id (if this transcript was manually annotated,\
 '-' otherwise)|
4 transcript-name|
5 gene-name|
6 sequence-length
*/
void transcript_info::adjust_cds_ranges_pepfile(const char* pfname)
{
  fasta_reader peptide_fa(pfname);
  for (size_t i = 0; i!= tlist.size(); ++i){
    //cout<<tlist[i].tid<<" "<<tlist[i].start<<"-"<<tlist[i].stop;
    int plen = sequenceLength(peptide_fa.faiIndex, i);
    if (tlist[i].phase!=0)
      plen -= 1;
    // adjust cds range based on peptide fasta sequence length
    tlist[i].plen = plen;
    tlist[i].start += tlist[i].phase;
    tlist[i].stop = tlist[i].start+plen*3-1;
  }
}

void transcript_info::adjust_cds_ranges()
{
  // fasta_reader trans_fa(tfname);
  for (size_t i = 0; i!= tlist.size(); ++i){
    //cout<<"before adjust: "<<tlist[i].tid<<" "<<tlist[i].start<<"-"<<tlist[i].stop;
    tlist[i].start += tlist[i].phase;
    int tlen = tlist[i].stop - tlist[i].start + 1;
    if (tlen%3 != 0)
      tlist[i].stop -= tlen%3;
    // skipping the first and the last codon
    // since they are likely to be start and stop codons
    // the current code do not check the actual sequence to verify this statement
    tlist[i].start += 3;
    tlist[i].stop -= 3;
    tlist[i].plen = (tlist[i].stop - tlist[i].start + 1)/3;
  }
}

fasta_reader::fasta_reader(const char* fname)
{
  // Try to read the FAI index.
  bool readSuccess = (read(faiIndex, fname) == 0);
  if (!readSuccess) {
    cerr << "Could not read the FAI index.  Not fatal, we can just build it.\n";
    // Try to build the FAI index (in memory) if reading was unsuccessful.  If
    if (build(faiIndex, fname) != 0) {
      cerr << "FATAL: Could not build FAI index.\n";
      exit(1);
    }
    // building into memory succeeded, we try to write it out.
    if (write(faiIndex) != 0) {
      cerr << "FATAL: Could not write out FAI index after building.\n";
      exit(1);
    }
    read(faiIndex, fname);
  }
}

string fasta_reader::read_region(unsigned refID, unsigned start, unsigned stop) const
{
  seqan::CharString seq_seqan;
  if (readRegion(seq_seqan, faiIndex, refID, start, stop)) return "";
  return string(toCString(seq_seqan));
}

string fasta_reader::read_seq(unsigned refID) const
{
  seqan::CharString seq_seqan;
  if (readSequence(seq_seqan,faiIndex, refID)) return "";
  return string(toCString(seq_seqan));
}

uint64_t fasta_reader::length(unsigned refID) const
{
  return sequenceLength(faiIndex, refID);
}
void encode_peptide(const string& tseq, string& pseq)
{
  for (size_t i=0; i<tseq.size(); i+=3){
    string c(tseq.substr(i,3));
    string aa(codon2aa.at(c));
    pseq += aa;
  }
}