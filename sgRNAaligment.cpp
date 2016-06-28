/*
 *    sgRNAalignment: aligment-based search for sgRNAs
 *
 *    Copyright (C) 2013-2015 Stanford University
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
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <GenomicRegion.hpp>
#include <RNG.hpp>
#include <smithlab_os.hpp>

char complement(char n)
{
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
  rev_comp.clear()
  rev_comp.resize(string_seq.size());
  for(size_t i = 0; i < seq.size(); i++){
    rev_comp[rev_comp.size() - i] = complement(seq[i]);
  }
}



size_t
propose_sgRNAs(const bool VERBOSE,
               const string PAM,
               const string region_seq,
               const size_t len_sgRNA;
               vector<string> &possible_sgRNAs){

  string region_seq_rev_comp;
  string rev_PAM;
  reverse_complement(region_seq, region_seq_rev_comp);
  reverse_complement(PAM, rev_PAM);
  
  size_t first_wildcard_pos = PAM.find_first_not_of("ACGT");
  if(first_wildcard_pos != string::npos){
    size_t last_wildcard_pos = PAM.find_last_not_of("ACGT");
    if(first_wildcard_pos != last_wildcard_pos){
      throw SMITHLABException("more than one wildcard found in PAM");
    }
    
   
    
  }
  else{
    // no wildcards
    for(size_t i = len_sgRNA; i < region_seq.size() - PAM.size(); i++){
      string test_seq = region_seq.substr(i, 3);
      if(test_seq == PAM){
        possible_sgRNAs.push_back(region_seq.substr(i - len_sgRNA, 23));
      }
      test_seq = region_seq_rev_comp.substr(i, 3);
      if(test_seq == rev_PAM){
        possible_sgRNAs.push_back(region_seq_rev_comp(i - len_sgRNA, 23));
      }
    }

  
}

int
main(const int argc, const char **argv) {
  
  try {
    // option variables
    size_t edit_dist = 2;
    size_t len_sgRNA = 20;
    string input_file_name;
    string output_file_name;
    string PAM_seq = "NGG";
    bool VERBOSE = false;
    
    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "-i input fasta file of proposed region",
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
                      + toa(length) + ")",
                      false, len_sgRNA);
    opt_parse.add_opt("input", 'i', "input file of the DNA sequence of the "
                      "regions to search in fasta format",
                      true, input_file_name);
    opt_parse.add_opt("genome", 'g', "genome file in fasta format",
                      true, genome_file_name);
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