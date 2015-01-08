#ifndef PROFILE_WRITER_HPP
#define PROFILE_WRITER_HPP
#include <string>

class transcript_info;
class ribo_profile;

void generate_profile_report(const std::string& fn_core, const ribo_profile& rprofile, const ribo_profile& mprofile, const transcript_info& tinfo);
void generate_base_report(const std::string& fname, const ribo_profile& rprofile, const ribo_profile& mprofile, const transcript_info& tinfo);
void generate_codon_report(const std::string& fname, const ribo_profile& rprofile, const ribo_profile& mprofile, const transcript_info& tinfo);
void generate_statistics_report(const std::string& fname, const ribo_profile& profiler, const transcript_info& tinfo);

#endif
