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



#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "profile_writer.hpp"
#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "utils.hpp"

void generate_profile_report(const string& fn_core, const ribo_profile& rprofile, const ribo_profile& mprofile, const transcript_info& tinfo)
{
  string fname = fn_core + ".base";
  generate_base_report(fname, rprofile, mprofile, tinfo);
  fname = fn_core + ".codon";
  generate_codon_report(fname, rprofile, mprofile, tinfo);
  fname = fn_core + ".stats";
  generate_statistics_report(fname, rprofile, tinfo);
}

void generate_base_report(const string& fname, const ribo_profile& rprofile, const ribo_profile& mprofile, const transcript_info& tinfo)
{
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  ofstream logfile(fname);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  for (size_t t=0; t!=rprofile.number_of_transcripts(); ++t) {
    double rabd = rprofile.get_tot_count(t);
    if (rabd==0) continue;
    rid_t refID(refID_vec[t]);
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
    logfile<<"refID: "<<refID<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID)<<endl;
    logfile<<"ribo profile: "<<rp<<endl;
    logfile<<"normalized ribo profile: "<<np<<endl;
    logfile<<"mRNA profile: "<<mp<<endl;
  }
  logfile.close();
}

vector<double> generate_codon_profile(rid_t refID, const ribo_profile& profiler, const transcript_info& tinfo)
{
  rid_t t = profiler.refID2pID.at(refID);
  int start(tinfo.cds_start(refID)), plen(tinfo.cds_pep_len(refID));
  vector<double> base_profile(profiler.get_read_assignments(t));
  vector<double> codon_profile(plen,0);
  for (size_t i=0; i!=codon_profile.size(); ++i)
    codon_profile[i] = base_profile[start + i*3];
  return codon_profile;
}

void generate_codon_report(const string& fname, const ribo_profile& rprofile, const ribo_profile& mprofile, const transcript_info& tinfo)
{
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  ofstream logfile(fname);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  for (size_t t=0; t!=rprofile.number_of_transcripts(); ++t) {
    double rabd = rprofile.get_tot_count(t);
    if (rabd==0) continue;
    rid_t refID(refID_vec[t]);
    vector<double> rp=generate_codon_profile(refID, rprofile, tinfo);
    vector<double> mp=generate_codon_profile(refID, mprofile, tinfo);
    // bias correction by normalizing ribo profile with mRNA profile
    vector<double> np(rp);
    for (size_t i=0; i!= np.size(); ++i) {
      if ( mp[i] != 0 )
        np[i] /= mp[i];
      else
        np[i] = 0;
    }
    logfile<<"refID: "<<refID<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID)<<endl;
    logfile<<"ribo profile: "<<rp<<endl;
    logfile<<"normalized ribo profile: "<<np<<endl;
    logfile<<"mRNA profile: "<<mp<<endl;
  }
  logfile.close();
}

void generate_statistics_report(const string& fname, const ribo_profile& profiler, const transcript_info& tinfo) 
{
  vector<rid_t> refID_vec(profiler.get_expressed_transcript_ids());
  ofstream logfile(fname);
  // logfile.setf(ios::scientific);
  // logfile.precision(3);
  for (size_t t=0; t!=profiler.number_of_transcripts(); ++t) {
    rid_t refID(refID_vec[t]);
    double rabd = profiler.get_tot_count(t);
    double tabd = profiler.get_tot_abundance(t) * profiler.len(t);
    logfile<<"refID: "<<refID<<endl;
    logfile<<"tid: "<<tinfo.get_tid(refID)<<endl;
    logfile<<"rabd: "<<rabd<<endl;
    logfile<<"tabd: "<<tabd<<endl;
    logfile<<"te: "<<rabd/tabd<<endl;
  }
  logfile.close();
}
