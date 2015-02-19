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
#include "reference_info_builder.hpp"
#include "utils.hpp"

int test_case()
{
  // const char* tfa("/home/hw1/scratch/ribojamdetector/transcriptome/protein_coding_33_filtered.fasta");
  // const char* pfa("/home/hw1/scratch/ribojamdetector/transcriptome/protein_coding_trans_33_filtered.fasta");
  // const char* cds("/home/hw1/scratch/ribojamdetector/transcriptome/cds_range.txt");
  // const char* tfa("/home/hw1/scratch/ribojamdetector/transcriptome/orf_coding.fasta");
  // const char* pfa("/home/hw1/scratch/ribojamdetector/transcriptome/orf_trans.fasta");
  // const char* cds("");
  const char* tfa("/home/hw1/scratch/ribomap/ref/human/gencode.v21.pc_transcripts_filter.fa");
  const char* pfa("/home/hw1/scratch/ribomap/ref/human/gencode.v21.pc_translations_filter.fa");
  const char* cds("/home/hw1/scratch/ribomap/ref/human/gencode.v21.pc_transcripts_cds.txt");

  transcript_info tinfo(tfa, cds);

  for (size_t i=0; i!=tinfo.total_count(); ++i) {
    cout<<tinfo.get_tid(i)<<" "<<tinfo.cds_start(i)<<"-"<<tinfo.cds_stop(i)<<endl;
    if (i>20) break;
  }
  
  // check whether encoding is right
  fasta_reader transcript_fa(tfa), peptide_fa(pfa);
  cout<<"transcriptome size: "<<tinfo.total_count()<<endl;
  cout<<"check whether encoding is right\n";
  int diff_cnt = 0;
  for (size_t i=0; i!=tinfo.total_count(); ++i){
    int start = tinfo.cds_start(i), stop = tinfo.cds_stop(i);
    int plen = tinfo.cds_pep_len(i);
    if (plen < 3) continue;
    string pconvert, pseq(peptide_fa.read_seq(i)), tseq(transcript_fa.read_region(i,start,stop));
    if (encode_peptide(tseq,pconvert)) {
      cout<<i<<" "<<transcript_fa.transcript_name(i)<<" "<<peptide_fa.transcript_name(i)<<endl;
      cout<<pseq<<endl;
      cout<<pconvert<<endl;
    }
    if (pseq.front() == 'X')
      pseq = pseq.substr(1);
    int l = std::min(pconvert.size(), pseq.size());
    pconvert = pconvert.substr(0, l);
    pseq = pseq.substr(0,l);
    if (pseq.compare(pconvert) != 0){
      cout<<i<<" "<<transcript_fa.transcript_name(i)<<" "<<peptide_fa.transcript_name(i)<<endl;
      cout<<pseq<<endl;
      cout<<pconvert<<endl;
      ++diff_cnt;
    }
  }
  cout<<"# transcripts that are encoded differently from the ground truth: "<<diff_cnt<<endl;
  return 0;

}

rid_t transcript_info::get_refID(const string& tid) const
{
  auto it = tid2refid.find(tid);
  if (it == tid2refid.end())
    return tid2refid.size();
  else
    return it->second; 
}

transcript_info::transcript_info(const char* tfa, const char* cds_range)
{
  get_info_from_fasta(tfa);
  build_tid_idx_map();
  if (cds_range)
    get_cds_range(cds_range);
}

bool transcript_info::get_info_from_fasta(const char* tfname)
{
  ifstream tfile(tfname);
  if (!tfile.good())
    return 1;
  seqan::RecordReader<ifstream, seqan::SinglePass<> > treader(tfile);
  seqan::CharString theader,tseq;
  while(!atEnd(treader)){
    if (readRecord(theader, tseq, treader, seqan::Fasta()) != 0){
      std::cerr<<"ERROR reading FASTA "<<tfname<<std::endl;
      return 1;
    }
    string tid(toCString(theader));
    size_t i(tid.find_first_of(" \t"));
    tid = tid.substr(0,i);
    int tlen = length(tseq);
    int start(0), end(tlen); // end position have to pass over the last position
    end -= (end-start)%3; 
    int plen = (end-start+1)/3;
    tlist.emplace_back(tprop(tid, start, end, plen, tlen));
  }
  build_tid_idx_map();
  return false;
}

void transcript_info::build_tid_idx_map()
{
  for(size_t i=0; i!=tlist.size(); ++i)
    tid2refid[tlist[i].tid] = i;
}

bool transcript_info::get_cds_range(const char* cds_fn)
{
  ifstream ifile(cds_fn);
  string tid;
  int start, stop;
  while(ifile >> tid >> start >> stop) {
    rid_t rid(get_refID(tid));
    if (rid==total_count()) {
      cerr<<tid<<" not found in transcript fasta!\n";
      continue;
    }
    set_cds_start(rid, start);
    set_cds_stop(rid, stop);
    set_cds_pep_len(rid, (stop-start+1)/3);
  }
  ifile.close();
  return false;
}

void transcript_info::get_gencode_info(const char* tfname, const char* gtf_fname)
{
  cout<<"read in transcriptome\n";
  get_gencode_info_from_fasta(tfname);
  build_tid_idx_map();
  cout<<"read in frame\n";
  get_info_from_gtf(gtf_fname);
  cout<<"adjust cds ranges\n";
  adjust_cds_ranges();
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
bool transcript_info::get_gencode_info_from_fasta(const char* tfname)
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
      if (startswith("CDS:", tmp_wrd)){
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
	  tlist[id].frame = record.phase - '0';
	}// if exon_num update
	//if (id==10) break;
      }//if find key
    }// if CDS
  }// while not eof
  //for (auto t: tlist)
  //  cout<<t.tid<<" "<<t.frame<<endl;
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
    if (tlist[i].frame!=0)
      plen -= 1;
    // adjust cds range based on peptide fasta sequence length
    tlist[i].plen = plen;
    tlist[i].start += tlist[i].frame;
    tlist[i].stop = tlist[i].start+plen*3-1;
    tlist[i].start -= 1; // gencode header is 1-based, fasta reader is 0-based
  }
}

void transcript_info::adjust_cds_ranges()
{
  // fasta_reader trans_fa(tfname);
  for (size_t i = 0; i!= tlist.size(); ++i){
    //cout<<"before adjust: "<<tlist[i].tid<<" "<<tlist[i].start<<"-"<<tlist[i].stop;
    tlist[i].start += tlist[i].frame;
    int tlen = tlist[i].stop - tlist[i].start + 1;
    if (tlen%3 != 0)
      tlist[i].stop -= tlen%3;
    // skipping the first and the last codon
    // since they are likely to be start and stop codons
    // the current code do not check the actual sequence to verify this statement
    tlist[i].start += 3;
    tlist[i].stop -= 3;
    tlist[i].plen = (tlist[i].stop - tlist[i].start + 1)/3;
    tlist[i].start -= 1; // gencode header is 1-based, fasta reader is 0-based
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

string fasta_reader::read_region(rid_t refID, unsigned start, unsigned stop) const
{
  seqan::CharString seq_seqan;
  if (readRegion(seq_seqan, faiIndex, refID, start, stop)) return "";
  return string(toCString(seq_seqan));
}

string fasta_reader::read_seq(rid_t refID) const
{
  seqan::CharString seq_seqan;
  if (readSequence(seq_seqan,faiIndex, refID)) return "";
  return string(toCString(seq_seqan));
}

uint64_t fasta_reader::length(rid_t refID) const
{
  return sequenceLength(faiIndex, refID);
}

string fasta_reader::transcript_name(rid_t refID) const
{
  return string(toCString(sequenceName(faiIndex, refID)));
}

bool encode_peptide(const string& tseq, string& pseq)
{
  for (size_t i=0; i<tseq.size(); i+=3){
    string c(tseq.substr(i,3));
    if (codon2aa.find(c) == codon2aa.end()) {
      cout<<c<<endl;
      return true;
    }
    string aa(codon2aa.at(c));
    pseq += aa;
  }
  return false;
}

