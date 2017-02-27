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
#include <QualityScore.hpp>



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
to_upper(string s){
  for(size_t i = 0; i < s.size(); i++)
    s[i] = to_upper(s[i]);
  return s;
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
  sgRNA(const string s, const string g, const size_t i)
  {seq = s; rev_comp = reverse_complement(s); target_name = g; index = i; count = 0;}
  
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
    std::istringstream is(buffer);
    if(!(is >> gene >> seq))
      throw SMITHLABException("bad line at " + toa(line_count) +
                              " " + buffer);
    seq = to_upper(seq);
    sgRNA guide = sgRNA(seq, gene, line_count);
    guides.push_back(guide);
    ++line_count;
  }
  return line_count;
}


bool
match_guide(const bool VERBOSE,
            vector<sgRNA> &guides,
            const string &input_string,
            const size_t guide_length,
            const size_t max_dist,
            SeedHash &forward_1st_hash,
            SeedHash &forward_2nd_hash,
            SeedHash &revcomp_1st_hash,
            SeedHash &revcomp_2nd_hash){
  vector<size_t> distances;
  vector<size_t> indexes;
  if(VERBOSE)
    cerr << "testing string " << input_string << endl;
  
  for(size_t i = 0; i < input_string.size() - guide_length; i++){
    // test for match in hash tables, if so then do full distance
    // 1st hash table
    
    bool EXACT_MATCH = false;
    string test_string = input_string.substr(i, guide_length/2);
    if(VERBOSE)
      cerr << test_string << endl;
    if(VERBOSE)
      cerr << "forward 1st half hash" << endl;
    std::pair<SeedHash::iterator, SeedHash::iterator>
      matches = forward_1st_hash.equal_range(test_string);
    if(matches.first != matches.second){
      if(VERBOSE)
        cerr << "match found" << endl;
      for(SeedHash::iterator it = matches.first;
          it != matches.second;){
      // do full distance
        string test_guide = it->second.seq;
        if(VERBOSE)
          cerr << "guide: " << it->second.seq << '\t'
               << it->second.target_name << '\t'
               << it->second.index << endl;
        size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
        
        if(dist <= max_dist){
          distances.push_back(dist);
          indexes.push_back(it->second.index);
        }
        if(dist == 0){
          EXACT_MATCH = true;
          break;
        }
        it++;
      }
    }
    if(EXACT_MATCH)
      break;
    
    // 2nd hash table
    test_string = input_string.substr(guide_length/2, guide_length/2);
    if(VERBOSE)
      cerr << test_string << endl;
    if(VERBOSE)
      cerr << "forward 2nd half hash" << endl;
    matches = forward_2nd_hash.equal_range(test_string);
    if(matches.first != matches.second){
      if(VERBOSE)
        cerr << "match found" << endl;
      for(SeedHash::iterator it = matches.first;
          it != matches.second;){
      // do full distance
        string test_guide = it->second.seq;
        if(VERBOSE)
          cerr << "guide: " << it->second.seq << '\t'
          << it->second.target_name << '\t'
          << it->second.index << endl;
        size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
        if(dist <= max_dist){
          distances.push_back(dist);
          indexes.push_back(it->second.index);
        }
        if(dist == 0){
          EXACT_MATCH = true;
          break;
        }
        it++;
      }
    }
    if(EXACT_MATCH)
      break;

    
    // 3rd hash table, reverse complement now
    test_string = reverse_complement(input_string.substr(0, guide_length/2));
    if(VERBOSE)
      cerr << test_string << endl;
    if(VERBOSE)
      cerr << "reverse complement 1st half hash" << endl;
    matches = revcomp_1st_hash.equal_range(test_string);
    if(matches.first != matches.second){
      if(VERBOSE)
        cerr << "match found" << endl;
      for(SeedHash::iterator it = matches.first;
          it != matches.second;){
      // do full distance
        string test_guide = it->second.seq;
        if(VERBOSE)
          cerr << "guide: " << it->second.seq << '\t'
          << it->second.target_name << '\t'
          << it->second.index << endl;
        size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
        if(dist <= max_dist){
          distances.push_back(dist);
          indexes.push_back(it->second.index);
        }
        if(dist == 0){
          EXACT_MATCH = true;
          break;
        }
        it++;
      }
    }
    if(EXACT_MATCH)
      break;
    
    // 4th hash table
    test_string = reverse_complement(input_string.substr(guide_length/2, guide_length/2));
    if(VERBOSE)
      cerr << test_string << endl;
    if(VERBOSE)
      cerr << "reverse complement 2nd half hash" << endl;
    matches = revcomp_1st_hash.equal_range(test_string);
    if(matches.first != matches.second){
      if(VERBOSE)
        cerr << "match found" << endl;
      for(SeedHash::iterator it = matches.first;
          it != matches.second;){
        // do full distance
        if(VERBOSE)
          cerr << "guide: " << it->second.seq << '\t'
          << it->second.target_name << '\t'
          << it->second.index << endl;
        string test_guide = it->second.seq;
        size_t dist = LevenshteinMetric(input_string.substr(i, guide_length), test_guide);
        if(dist <= max_dist){
          distances.push_back(dist);
          indexes.push_back(it->second.index);
        }
        if(dist == 0){
          EXACT_MATCH = true;
          break;
        }
        it++;
      }
    }
    if(EXACT_MATCH)
      break;
  }
  
  // one match, update counts
  if(distances.size() == 1){
    guides[indexes[0]].count++;
    return true;
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
    if(IS_UNIQUE){
      guides[indexes[which_min]].count++;
      return true;
    }
  }
  return false;
}

// from smithlab_os
\

inline bool is_fastq_name_line(size_t line_count) {
  return ((line_count & 3ul) == 0ul);
}

inline bool is_fastq_sequence_line(size_t line_count) {
  return ((line_count & 3ul) == 1ul);
}

inline bool is_fastq_score_name_line(size_t line_count) {
  return ((line_count & 3ul) == 2ul);
}

inline bool is_fastq_score_line(size_t line_count) {
  return ((line_count & 3ul) == 3ul);
}

void
read_fastq_file(const string &input_file_name,
                vector<string> &names,
                vector<string> &sequences,
                vector<string> &scores) {
  
  static const size_t INPUT_BUFFER_SIZE = 1000000;
  
  std::ifstream in(input_file_name.c_str());
  if (!in) //if file doesn't open
    throw SMITHLABException("could not open input file: " + input_file_name);
  
  string s, name, scr;
  bool first_line = true;
  bool is_sequence_line = false, is_score_line = false;
  size_t line_count = 0;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException(
                              "Line in " + name + "\nexceeds max length: "
                              + toa(INPUT_BUFFER_SIZE));
    if (in.gcount() == 0)
      break;
    
    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    
    if (is_fastq_name_line(line_count)) {
      if (buffer[0] != '@')
        throw SMITHLABException("invalid FASTQ name line: " + string(buffer));
      if (first_line == false && s.length() > 0) {
        names.push_back(name);
        sequences.push_back(s);
        scores.push_back(scr);
      } else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("@ "));
      is_sequence_line = true;
    }
    if (is_fastq_sequence_line(line_count)) {
      assert(is_sequence_line);
      s = buffer;
      is_sequence_line = false;
    }
    if (is_fastq_score_name_line(line_count)) {
      if (buffer[0] != '+')
        throw SMITHLABException("invalid FASTQ score name line: " +
                                string(buffer));
      is_score_line = true;
    }
    if (is_fastq_score_line(line_count)) {
      assert(is_score_line);
      scr = buffer;
      is_score_line = false;
    }
    ++line_count;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
    scores.push_back(scr);
  }
}

void
write_guides(const string &counts_file_name,
             const vector<sgRNA> &guides){
  std::ofstream of;
  if (!counts_file_name.empty()) of.open(counts_file_name.c_str());
  std::ostream out(counts_file_name.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  for(size_t i = 0; i < guides.size(); i++){
    out << guides[i].target_name << '\t'
        << guides[i].seq << '\t' << guides[i].count << endl;
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
    
    vector<string> names;
    vector<string> sequences;
    vector<string> scores;
    if(VERBOSE)
      cerr << "reading in reads" << endl;
    read_fastq_file(fastq_file_name, names, sequences, scores);
    if(VERBOSE)
      cerr << "read in " << sequences.size() << " reads" << endl;
    if(VERBOSE)
      cerr << "reading in guides" << endl;
    vector<sgRNA> guides;
    read_guides_from_text(guide_file_name, guides);
    if(VERBOSE){
      cerr << "read in " << guides.size() << " guides" << endl;
      //for(size_t i = 0; i < guides.size(); i++){
      //  cerr << guides[i].target_name << '\t'
      //       << guides[i].seq << '\t' << guides[i].count << endl;
      //}
    }
    
    size_t guide_length = guides[0].seq.size();
    SeedHash forward_1st_hash, forward_2nd_hash, revcomp_1st_hash, revcomp_2nd_hash;
    if(VERBOSE)
      cerr << "building hash tables" << endl;
    build_hashes(guide_length, guides, forward_1st_hash, forward_2nd_hash,
                 revcomp_1st_hash, revcomp_2nd_hash);
    
    if(VERBOSE)
      cerr << "matching reads to guides" << endl;
    size_t n_matches = 0;
    for(size_t i = 0; i < sequences.size(); i++){
      bool MATCH_FOUND =
        match_guide(false, guides, sequences[i], guide_length,  max_dist,
                  forward_1st_hash, forward_2nd_hash,
                  revcomp_1st_hash, revcomp_2nd_hash);
      if(MATCH_FOUND)
        n_matches++;
      if(VERBOSE && (i % 1000000 == 0))
        cerr << "processed " << i << " reads" << endl;
    }
    
    if(VERBOSE){
      cerr << "# matches found = " << n_matches << endl;
      cerr << "percent mapped = " << static_cast<double>(n_matches)/sequences.size() << endl;
    }
    write_guides(counts_file_name, guides);
    
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