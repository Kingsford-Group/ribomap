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
