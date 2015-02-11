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



#include <map>
#include <set>
#include <vector>
#include <string>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <seqan/bam_io.h>

#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "bam_parser.hpp"
#include "utils.hpp"

void alignment_regions_to_codon_ranges(const rd_rec_map_t& fp_rec, const transcript_info& tinfo, fp_list_t& fp_codon_list, int offset)
{ 
  for (auto fp: fp_rec){
    vector<position> codon_list;
    for (auto bp: fp.second.al_loci){
      unsigned refID = bp.refID;
      int cds_begin = tinfo.cds_start(refID);
      int cds_end = tinfo.cds_stop(refID);
      int phase = tinfo.frame(refID);
      position cp;
      cp.type = FRAME0;
      read_t read_type; 
      if (offset==-1)
	read_type = base_range_to_codon_range(bp, cds_begin, cds_end, phase, cp);
      else
	read_type = base_range_to_middle_codon(bp, cds_begin, cds_end, offset, cp);
      if (read_type == FRAME0)
	codon_list.push_back(cp);
    }
    if (codon_list.size()!=0)
      fp_codon_list.emplace_back(fp_record{fp.second.count, codon_list, fp.second.seqs, fp.second.used});
  }
  // sanity check whether fp_rec is the same as bam summary
  double total=0, multi_mapped=0;
  for (auto rec: fp_codon_list){
    total += rec.count;
    if (rec.al_loci.size()>1)
      multi_mapped += rec.count;
  }
  cout<<"reads mapped to cds: "<<fp_codon_list.size()<<endl;
  printf("total: %.0f\tmulti_mapped: %.0f (%.2f %%)\n",total, multi_mapped, multi_mapped*100/total);
}

void assign_P_site(const rd_rec_map_t& fp_rec_in, const transcript_info& tinfo, fp_list_t& fp_rec_out, int offset)
{ 
  for (auto& ifp: fp_rec_in){
    vector<position> obp_list;
    for (auto& ibp: ifp.second.al_loci){
      unsigned refID = ibp.refID;
      int cds_begin = tinfo.cds_start(refID);
      int cds_end = tinfo.cds_stop(refID);
      position obp(ibp);
      if (offset==-1)
	read_type_from_range(ibp, cds_begin, cds_end, obp);
      else
	read_type_of_psite(ibp, cds_begin, cds_end, offset, obp);
      obp_list.push_back(obp);
    }
    fp_rec_out.emplace_back(fp_record{ifp.second.count, obp_list, ifp.second.seqs, ifp.second.used});
  }
  // sanity check whether fp_rec is the same as bam summary
  double total=0, multi_mapped=0;
  for (auto rec: fp_rec_out){
    total += rec.count;
    if (rec.al_loci.size()>1)
      multi_mapped += rec.count;
  }
  cout<<"total output footprint: "<<fp_rec_out.size()<<endl;
  printf("total: %.0f\tmulti_mapped: %.0f (%.2f %%)\n",total, multi_mapped, multi_mapped*100/total);
}

void assign_P_site(const rd_rec_map_t& fp_rec_in, const transcript_info& tinfo, fp_list_t& fp_rec_out, rlen2psite_t rl2p)
{
  for (auto& ifp: fp_rec_in){
    vector<position> obp_list;
    for (auto& ibp: ifp.second.al_loci){
      unsigned refID = ibp.refID;
      int cds_begin = tinfo.cds_start(refID);
      int cds_end = tinfo.cds_stop(refID);
      position obp(ibp);
      auto it = rl2p.find(ibp.stop-ibp.start);
      if (it==rl2p.end()) continue;
      read_type_of_psite(ibp, cds_begin, cds_end, it->second, obp);
      obp_list.push_back(obp);
    }
    fp_rec_out.emplace_back(fp_record{ifp.second.count, obp_list, ifp.second.seqs, ifp.second.used});
  }
  // sanity check whether fp_rec is the same as bam summary
  double total=0, multi_mapped=0;
  for (auto rec: fp_rec_out){
    total += rec.count;
    if (rec.al_loci.size()>1)
      multi_mapped += rec.count;
  }
  cout<<"total output footprint: "<<fp_rec_out.size()<<endl;
  printf("total: %.0f\tmulti_mapped: %.0f (%.2f %%)\n",total, multi_mapped, multi_mapped*100/total);
}

/* assign read type based on read range to output position and return read type */
read_t read_type_from_range(const position& ibp, int cds_begin, int cds_end, position& obp)
{
  obp.type = FRAME0;
  if (ibp.start < cds_begin)
    obp.type =  UTR5;
  else if (ibp.stop > cds_end)
    obp.type = UTR3;
  return obp.type;
}

/* get P-site from reads, asign read type based on P-site, and return read type */
read_t read_type_of_psite(const position& ibp, int cds_begin, int cds_end, int offset, position& obp)
{
  int psite = ibp.start + offset;
  // transcripts on the reverse complement strand has been flipped
  // therefore no need to consider the transcript strand
  if (psite < cds_begin)
    obp.type =  UTR5;
  // last codon position --> last base - 1 (zero index)
  else if (psite > cds_end-1)
    obp.type = UTR3;
  else {
    switch ((psite-cds_begin)%3) {
    case 0:
      obp.type = FRAME0;
      break;
    case 1:
      obp.type = FRAME1;
      break;
    case 2:
      obp.type = FRAME2;
      break;
    }
  }
  obp.start = psite;
  obp.stop = obp.start+1;
  return obp.type;
}

/*convert read mapping of transcript DNA base ranges to codon ranges */
read_t base_range_to_codon_range(const position& bp, int cds_begin, int cds_end, int phase, position& cp)
{
  //both alignment and codon_position are zero-based
  int pos_begin = bp.start;
  int pos_end = bp.stop;
  // transcripts on the reverse complement strand has been flipped
  // therefore no need to consider the transcript strand
  if (pos_end < cds_begin)
    return UTR5;
  // last codon position --> last base - 1 (zero index)
  else if (pos_begin > cds_end-1)
    return UTR3;
  else {
    // only consider read range within the cds
    pos_begin = max(pos_begin, cds_begin);
    pos_end = min(pos_end, cds_end);
    // // only consider codons that are entirely within the read range
    // pos_begin += (3-(pos_begin-cds_begin)%3)%3;
    // pos_end -= (pos_end-cds_begin)%3;
    int codon_begin = (pos_begin - cds_begin)/3;
    int codon_end = (pos_end - cds_begin)/3;
    if (phase!=0){
      codon_begin = max(codon_begin+1, cds_begin);
      codon_end = min(codon_end+1, cds_end);
    }
    cp.refID = bp.refID;
    cp.start = codon_begin;
    cp.stop = codon_end;
    return FRAME0;
  }
}

/*convert read mapping of transcript DNA base ranges to the middle codon */
read_t base_range_to_middle_codon(const position& bp, int cds_begin, int cds_end, int offset, position& cp)
{
  int middle = bp.start + offset;
  // transcripts on the reverse complement strand has been flipped
  // therefore no need to consider the transcript strand
  if (middle < cds_begin)
    return UTR5;
  // last codon position --> last base - 1 (zero index)
  if (middle > cds_end-1)
    return UTR3;
  // if middle codon not in frame, map to the closest codon
  int codon_begin;
  switch ((middle-cds_begin)%3) {
  case 2:
    codon_begin = (middle-cds_begin)/3 + 1;
    break;
  case 0:
  case 1:
  default:
    codon_begin = (middle-cds_begin)/3;
    break;
  }
  if (codon_begin >= (cds_end-cds_begin)/3 )
    return UTR3;
  else {
    int codon_end = codon_begin + 1;
    cp.refID = bp.refID;
    cp.start = codon_begin;
    cp.stop = codon_end;
    return FRAME0;
  }
}

bool expressed_read_codon_ranges_from_bam(fp_list_t& fp_codon_list, const char *fn, const transcript_info& tinfo, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC, int offset)
{
  cout<<"getting alignment records..."<<endl;
  rd_rec_map_t rd_rec;
  get_expressed_alignments_from_bam(rd_rec, fn, profiler, cnt_sep, lmin, lmax, useSecondary, useRC);
  cout<<"total number of reads: "<<rd_rec.size()<<endl;
  cout<<"convert read loci to codon ranges...\n";
  alignment_regions_to_codon_ranges(rd_rec, tinfo, fp_codon_list, offset);
  return false;
}

bool expressed_read_bases_from_bam(fp_list_t& fp_rec_out, const char* fn, const transcript_info& tinfo, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC, int offset)
{
  cout<<"getting alignment records..."<<endl;
  rd_rec_map_t rd_rec;
  get_expressed_alignments_from_bam(rd_rec, fn, profiler, cnt_sep, lmin, lmax, useSecondary, useRC);
  cout<<"total number of reads: "<<rd_rec.size()<<endl;
  cout<<"getting read type and p-sites...\n";
  assign_P_site(rd_rec, tinfo, fp_rec_out, offset);
  return false;
}

bool expressed_read_bases_from_bam(fp_list_t& fp_rec_out, const char* fn, const transcript_info& tinfo, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC, const char* offset_fn)
{
  cout<<"getting readlen mapping to P site offset..."<<endl;
  rlen2psite_t rl2p = get_readlen_psite_map(offset_fn);
  cout<<"getting alignment records..."<<endl;
  rd_rec_map_t rd_rec;
  get_expressed_alignments_from_bam(rd_rec, fn, profiler, cnt_sep, lmin, lmax, useSecondary, useRC);
  cout<<"total number of reads: "<<rd_rec.size()<<endl;
  cout<<"getting read type and p-sites...\n";
  assign_P_site(rd_rec, tinfo, fp_rec_out, rl2p);
  return false;
}

rlen2psite_t get_readlen_psite_map(const char* offset_fn)
{
  ifstream ifile(offset_fn);
  rlen2psite_t rl2p;
  int readlen, offset;
  while(ifile >> readlen >> offset)
    rl2p[readlen] = offset;
  ifile.close();
  return rl2p;
}

int get_count_from_fasta_header(const string header, const string& sep="_")
{
  vector<string> words;
  int count(1);
  string_split(header, sep, words);
  if (words.size()==1) return 1;
  else {
    try { count = std::stoi(words.back().c_str()); }
    catch (const std::invalid_argument& ia) {
      return 1;
    }
  }
  return count;
}

bool get_expressed_alignments_from_bam(rd_rec_map_t& rd_rec, const char *fn, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC)
{
  // open bam file
  seqan::BamStream bamIn(fn);
  if (!isGood(bamIn)){
    cerr << "ERROR: Could not open "<<fn<<"!"<<endl;
    exit(1);
  }
  seqan::BamAlignmentRecord bam_rec;
  while(!atEnd(bamIn)){
    if(readRecord(bam_rec,bamIn)!=0){
      cerr << "ERROR: Could not read bam record!\n";
      return true;
    }
    // get mapped reads
     if (!hasFlagUnmapped(bam_rec)){
       if (not useSecondary and hasFlagSecondary(bam_rec)) continue;
       // the reverse complement of a read map to the transcriptome
       // however ribo-seq is stranded, such case is therefore impossible
       // therefore such alignments should be discarded
       if (hasFlagRC(bam_rec) and not useRC) continue;
      // get alignment information
      unsigned refID = bam_rec.rID;
      if ( not profiler.is_expressed(refID)) continue;
      // get read name and sequence
      string name(toCString(bam_rec.qName));
      string seq(toCString(bam_rec.seq));
      //bool strand = !hasFlagRC(bam_rec);
      // build a list of mapped regions
      // length>1 if have spliced mappings
      // --> length(cigarElement) > 1
      int pos_begin(bam_rec.beginPos);
      // read length is the length on the transcript covered by the read
      int read_len(0);
      // parse cigar string to figure out read length
      bool skip(false);
      for (int i=0; i<length(bam_rec.cigar); ++i) {
	switch (bam_rec.cigar[i].operation) {
	case 'M':
	  // add match length to read length
	  read_len += bam_rec.cigar[i].count;
	  break;
	case 'D':
	  // Deletion is a gap to the read
	  // read covered one base longer on the transcript
	  read_len += 1;
	  break;
	case 'N':
	  // skipped regions indicate novel splices
	  // since reads are mapped to the transcriptome
	  // skip these alignments for now
	  skip = true;
	  break;
	case 'I':
	  // Insertion is a gap to the reference
	  // no change in read length
	case 'S':
	  // softclipping won't change read length
	default:
	  break;
	}
	if (skip) break;
      }
      if (skip or read_len<lmin or read_len>lmax) continue;
      int pos_end = pos_begin + read_len;
      position p{refID, pos_begin, pos_end, UNKNOWN};
      auto it = rd_rec.find(name);
      if (it == rd_rec.end()) {
	int count(get_count_from_fasta_header(name, cnt_sep));
	rd_rec.emplace(name,fp_record{count, vector<position>{p}, set<string>{seq}, false});
	//cout<<name<<" "<<count<<endl;
      }
      else {
	it->second.seqs.insert(seq);
	it->second.al_loci.push_back(p);
      }
    }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  return false;
}
