/*
 *    sgRNAalignment: aligment-based search for sgRNAs
 *
 *    Copyright (C) 2016 Stanford University
 *             Timothy Daley
 *
 *    Authors: Timothy Daley
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>
#include <queue>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <tr1/unordered_map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>



#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>
#include <MappedRead.hpp>


using std::string;
using std::vector;
using std::min;
using std::endl;
using std::max;
using std::cerr;
using std::tr1::unordered_multimap;
using std::tr1::unordered_map;
using std::strcmp;
using std::make_pair;


inline bool
is_genomic_base(char n){
  switch(n){
    case 'A':
      return true;
    case 'C':
      return true;
    case 'G':
      return true;
    case 'T':
      return true;
    default:
      return false;
  }
}

inline char
complement(char n){
  switch(n){
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    case 'N':
      return 'N';
  }
  assert(false);
  return ' ';
}

inline char
to_upper(char n){
  switch(n){
    case 'a':
      return 'A';
    case 'c':
      return 'C';
    case 'g':
      return 'G';
    case 't':
      return 'T';
    case 'A':
      return 'A';
    case 'C':
      return 'C';
    case 'G':
      return 'G';
    case 'T':
      return 'T';
    default:
      return 'N';
  }
}

string
reverse_complement(const string seq){
  string rev_comp;
  rev_comp.resize(seq.size());
  for(size_t i = 0; i < seq.size(); i++){
    rev_comp[seq.size() - 1  - i] = complement(to_upper(seq[i]));
  }
  return rev_comp;
}

struct sgRNA {
  sgRNA () {}
  sgRNA(const string s, const string g)
  {seq = s; rev_comp = reverse_complement(s); target_name = g; count = 0;}
  
  string seq;
  string rev_comp;
  string target_name;
  size_t count;
  size_t index;
};


typedef unordered_multimap<string, sgRNA> SeedHash;

void
build_hashes(const size_t guide_length,
             vector<sgRNA> &guides,
             SeedHash &forward_1st_hash,
             SeedHash &forward_2nd_hash,
             SeedHash &revcomp_1st_hash,
             SeedHash &revcomp_2nd_hash){
  for(size_t i = 0; i < guides.size(); i++){
    string first_half = guides[i].seq.substr(0, guide_length/2);
    forward_1st_hash.insert(make_pair<string, sgRNA>(first_half, guides[i]));
    string first_half_rev_comp = reverse_complement(first_half);
    revcomp_1st_hash.insert(make_pair<string, sgRNA>(first_half_rev_comp, guides[i]));
  }
  for(size_t i = 0; i < guides.size(); i++){
    string second_half = guides[i].seq.substr(guide_length/2, guide_length/2);
    forward_2nd_hash.insert(make_pair<string, sgRNA>(second_half, guides[i]));
    string second_half_rev_comp =
      reverse_complement(second_half);
    revcomp_2nd_hash.insert(make_pair<string, sgRNA>(second_half_rev_comp, guides[i]));
  }
}

// wildcard = 'N'
int
LevenshteinMetric(const string &s1,
                  const string &s2){
  vector< vector<int> > dist(s1.size() + 1, vector<int>(s2.size() + 1));
  //initialization
  dist[0][0] = 0;
  for(size_t i = 1; i < dist.size(); i++){
    dist[i][0] = dist[i-1][0] + 1;
  }
  for(size_t j = 1; j < dist[0].size(); j++){
    dist[0][j] = dist[0][j-1] + 1;
  }
  // now work
  for(size_t i = 1; i < dist.size(); i++){
    for(size_t j = 1; j < dist[0].size(); j++){
      dist[i][j] = min(dist[i-1][j] + 1,
                       min(dist[i][j-1] + 1,
                           dist[i-1][j-1] +
                           ((to_upper(s1[i-1]) == to_upper(s2[j-1])) ? 0:1)));
    }
  }
  return dist[s1.size()][s2.size()];
}



size_t
read_fasta_batch(const string &input_file_name,
                 vector<string> &seqs,
                 vector<string> &names){
  std::ifstream in(input_file_name.c_str());
  if (!in) //if file doesn't open
    throw SMITHLABException("could not open input file: " + input_file_name);
  size_t n_regions = 0;
  size_t line_count = 0ul;
  string buffer;
  string current_seq;
  while(getline(in, buffer)){
    ++line_count;
    std::istringstream is(buffer);
    // > marks the name portions
    if(buffer[0] == '>'){
      names.push_back(buffer.substr(1, buffer.length() - 1));
      // first one do nothing, otherwise append the current sequence
      if(line_count > 1){
        seqs.push_back(current_seq);
        current_seq.clear();
        ++n_regions;
      }
    }
    else{
      current_seq.append(buffer);
    }
  }
  // need to append final sequence
  ++n_regions;
  seqs.push_back(current_seq);
  
  return n_regions;
}

// input file should be a text file with gene name then guide sequence
size_t
read_guides_from_text(const string &input_file_name,
                      vector<sgRNA> &guides){
  std::ifstream in(input_file_name.c_str());
  if (!in) //if file doesn't open
    throw SMITHLABException("could not open input file: " + input_file_name);
  size_t line_count = 0;
  string buffer;
  string seq;
  string gene;
  while(getline(in, buffer)){
    ++line_count;
    std::istringstream is(buffer);
    if(!(is >> gene >> seq))
      throw SMITHLABException("bad line at " + toa(line_count) +
                              " " + buffer);
    sgRNA guide = sgRNA(seq, gene);
    guides.push_back(guide);
  }
  return line_count;
}


void
match_guide(vector<sgRNA> &guides,
            const string &input_string,
            const size_t guide_length,
            const size_t max_dist,
            SeedHash &forward_1st_hash,
            SeedHash &forward_2nd_hash,
            SeedHash &revcomp_1st_hash,
            SeedHash &revcomp_2nd_hash){
  vector<size_t> distances;
  vector<size_t> indexes;
  for(size_t i = 0; i < input_string.size() - guide_length; i++){
    // test for match in hash tables, if so then do full distance
    // 1st hash table
    string test_string = input_string.substr(i, guide_length/2);
    std::pair<SeedHash::iterator, SeedHash::iterator>
      matches = forward_1st_hash.equal_range(test_string);
    for(SeedHash::iterator it = matches.first;
        it != matches.second;){
      // do full distance
      string test_guide = it->second.seq;
      size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
      if(dist <= max_dist){
        distances.push_back(dist);
        indexes.push_back(it->second.index);
      }
    }
    
    // 2nd hash table
    test_string = input_string.substr(guide_length/2, guide_length/2);
    matches = forward_2nd_hash.equal_range(test_string);
    for(SeedHash::iterator it = matches.first;
        it != matches.second;){
      // do full distance
      string test_guide = it->second.seq;
      size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
      if(dist <= max_dist){
        distances.push_back(dist);
        indexes.push_back(it->second.index);
      }
    }
    
    // 3rd hash table, reverse complement now
    test_string = reverse_complement(input_string.substr(0, guide_length/2));
    matches = revcomp_1st_hash.equal_range(test_string);
    for(SeedHash::iterator it = matches.first;
        it != matches.second;){
      // do full distance
      string test_guide = it->second.seq;
      size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
      if(dist <= max_dist){
        distances.push_back(dist);
        indexes.push_back(it->second.index);
      }
    }
    
    // 4th hash table
    test_string = reverse_complement(input_string.substr(guide_length/2, guide_length/2));
    matches = revcomp_1st_hash.equal_range(test_string);
    for(SeedHash::iterator it = matches.first;
        it != matches.second;){
        // do full distance
      string test_guide = it->second.seq;
      size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
      if(dist <= max_dist){
        distances.push_back(dist);
        indexes.push_back(it->second.index);
      }
    }
    // one match, update counts
    if(distances.size() == 1){
      guides[indexes[0]].count++;
    }
    else if(distances.size() > 1){
      size_t which_min = 0;
      bool IS_UNIQUE = true;
      for(size_t j = 1; j < distances.size(); j++){
        if(distances[j] < distances[which_min]){
          which_min = j;
          IS_UNIQUE = true;
        }
        else if(distances[j] == distances[which_min]){
          IS_UNIQUE = false;
        }
      }
      if(IS_UNIQUE)
        guides[indexes[which_min]].count++;
    }
  }
}



int
main(const int argc, const char **argv) {
  
  try {
    // option variables
    size_t max_dist = 1;
    string guide_file_name;
    string fastq_file_name;
    string counts_file_name;
    bool VERBOSE = false;
    
    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "-g guide file name"
                           " -f fastq file name");
    opt_parse.add_opt("max_dist", 'd', "maximum edit distance for matches"
                      " (default: " + toa(max_dist) + ")",
                      false, max_dist);
    opt_parse.add_opt("guides", 'g', "file name of the guides and associated genes, "
                      "tab-separated with gene first",
                      true, guide_file_name);
    opt_parse.add_opt("fastq", 'f', "sequenced fastq file name",
                      true, fastq_file_name);
    opt_parse.add_opt("output", 'o', "output file name",
                      false, counts_file_name);
    opt_parse.add_opt("VERBOSE", 'V', "verbose mode",
                      false, VERBOSE);
    
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || argc == 0 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    /**********************************************************************/
    

    
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}