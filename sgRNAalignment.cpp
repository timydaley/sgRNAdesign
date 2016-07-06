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


#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>


using std::string;
using std::vector;
using std::min;
using std::endl;
using std::max;
using std::cerr;
using std::tr1::unordered_map;

size_t
read_fasta(const string &input_file_name,
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
  
char
complement(char n){
  switch(n)
  {
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

void
reverse_complement(const string seq,
                   string &rev_comp){
  rev_comp.clear();
  rev_comp.resize(seq.size());
  for(size_t i = 0; i < seq.size(); i++){
    rev_comp[seq.size() - 1  - i] = complement(seq[i]);
  }
}

bool
wildcard_seq_match(const char *first_seq,
                   const char *second_seq){
  if((*first_seq == '\0' && *second_seq != '\0')
     || (*first_seq != '\0' && second_seq == '\0')){
    throw SMITHLABException("length of strings don't match");
  }
  // reach end of string, done
  if(*first_seq == '\0' && *second_seq == '\0'){
    return true;
  }
  // match
  else if (*first_seq == *second_seq){
    return wildcard_seq_match(first_seq + 1, second_seq + 1);
  }
  // wildcard
  else if (*first_seq == 'N' || *second_seq == 'N'){
    return wildcard_seq_match(first_seq + 1, second_seq + 1);
  }
  // else no match
  return false;
}



void
propose_sgRNAs(const bool VERBOSE,
               const string PAM,
               const string region_seq,
               const size_t len_sgRNA,
               vector<string> &possible_sgRNAs){

  string region_seq_rev_comp;
  string rev_PAM;
  reverse_complement(region_seq, region_seq_rev_comp);
  reverse_complement(PAM, rev_PAM);
  const size_t PAM_len = PAM.size();
  
  // test for wildcards
  size_t first_wildcard_pos = PAM.find_first_not_of("ACGT");
  if(first_wildcard_pos == std::string::npos){
    if(VERBOSE){
      cerr << "no wildcard detected" << endl;
    }
    // no wildcards, just do simple matching
    for(size_t i = len_sgRNA; i < region_seq.size() - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(test_seq == PAM){
        possible_sgRNAs.push_back(region_seq.substr(i - len_sgRNA, len_sgRNA + PAM_len));
      }
      test_seq = region_seq_rev_comp.substr(i, PAM_len);
      if(test_seq == rev_PAM){
        possible_sgRNAs.push_back(region_seq_rev_comp.substr(i, len_sgRNA + PAM_len));
      }
    }
  }
  else{
    if(VERBOSE){
      cerr << "wildcard detected" << endl;
    }
    size_t last_wildcard_pos = PAM.find_last_not_of("ACGT");
    // test for more than one wildcard
    if(first_wildcard_pos != last_wildcard_pos){
      throw SMITHLABException("more than one wildcard found in PAM");
    }
    
    for(size_t i = len_sgRNA; i < region_seq.size() - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(wildcard_seq_match(test_seq.c_str(), PAM.c_str())){
        possible_sgRNAs.push_back(region_seq.substr(i - len_sgRNA, len_sgRNA + PAM_len));
      }
    }
    for(size_t i = 0; i < region_seq_rev_comp.size() - len_sgRNA - PAM_len; i++){
      string test_seq = region_seq_rev_comp.substr(i, PAM_len);
      if(wildcard_seq_match(test_seq.c_str(), rev_PAM.c_str())){
        possible_sgRNAs.push_back(region_seq_rev_comp.substr(i, len_sgRNA + PAM_len));
      }
    }
  }
}

int
main(const int argc, const char **argv) {
  
  try {
    // option variables
    size_t edit_dist = 2;
    size_t len_sgRNA = 20;
    string input_file_name;
    string genome_file_name;
    string output_file_name;
    string PAM_seq = "NGG";
    bool VERBOSE = false;
    
    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "-i input fasta file of proposed region"
                           " -g input genome file");
    opt_parse.add_opt("output", 'o', "output file",
                      false , output_file_name);
    opt_parse.add_opt("edit_distance", 'd', "maximum edit distance for matches"
                      " (default: " + toa(edit_dist) + ")",
                      false, edit_dist);
    opt_parse.add_opt("PAM", 'P', "PAM sequence (default: "
                      + PAM_seq + ")",
                      false, PAM_seq);
    opt_parse.add_opt("length", 'l', "length of sgRNAs (default: "
                      + toa(len_sgRNA) + ")",
                      false, len_sgRNA);
    opt_parse.add_opt("input", 'i', "input file of the DNA sequence of the "
                      "regions to search in fasta format",
                      true, input_file_name);
    opt_parse.add_opt("genome", 'g', "genome file in fasta format",
                      false, genome_file_name);
    opt_parse.add_opt("VERBOSE", 'V', "verbose mode",
                      false, VERBOSE);
    
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    
    if(VERBOSE){
      cerr << "PAM = " << PAM_seq << endl;
    }
    vector<string> seqs, names;
    size_t n_seqs = read_fasta(input_file_name, seqs, names);
    
    if(VERBOSE){
      cerr << "sequences read in: " << n_seqs << endl;
  
      for(size_t i = 0; i < min(seqs.size(), names.size()); i++){
        cerr << names[i] << endl;
        cerr << seqs[i] << endl;
      }
    }
    
    vector<string> possible_sgRNAs;
    propose_sgRNAs(VERBOSE, PAM_seq, seqs[0], len_sgRNA, possible_sgRNAs);
    
    if(VERBOSE){
      cerr << "proposed sgRNAs:" << endl;
      for(size_t i = 0; i < possible_sgRNAs.size(); i++){
        cerr << possible_sgRNAs[i] << endl;
      }
    }

    

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