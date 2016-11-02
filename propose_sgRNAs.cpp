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


struct sgRNA {
  sgRNA () {}
  sgRNA(const string s, const bool rc)
    {seq = s; rev_comp = rc; to_delete = false;}
  sgRNA(const string s, const string n, const bool rc)
    {seq = s; target_name = n; rev_comp = rc; to_delete = false;}
  sgRNA(const string s, const size_t k, const bool rc)
    {seq = s; key = k; rev_comp = rc; to_delete = false;}
  sgRNA(const string s, const string n, const size_t k, const bool rc)
    {seq = s; target_name = n; key = k; rev_comp = rc; to_delete = false;}

  
  string seq;
  string target_name;
  int key;
  bool rev_comp;
  bool to_delete;
  vector<MappedRead> matches;
};

bool
compare_sgRNA(const sgRNA &one, const sgRNA &two){
  return (one.seq < two.seq);
}

bool
equal_sgRNA_seq(const sgRNA &one, const sgRNA &two){
  return (one.seq == two.seq);
}

std::ostream&
operator<<(std::ostream& the_stream, const sgRNA &sgrna){

  the_stream << "sgRNA: " << sgrna.seq << endl;
  the_stream << "matches: " << endl;
  for(size_t i = 0; i < sgrna.matches.size(); i++){
    the_stream << sgrna.matches[i] << endl;
  }
  
  
  return the_stream;
}

typedef unordered_multimap<size_t, sgRNA> SeedHash;

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

/*
// from smithlab_utils.hpp
inline size_t
base2int(char c) {
  switch(c) {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    case 'a' : return 0;
    case 'c' : return 1;
    case 'g' : return 2;
    case 't' : return 3;
    default  : return 4;
  }
}
 */

// convert whole string to integer value for Rabin-Karp
inline size_t
string2int(const string seq){
  // sum_{i = 0}^{|seq| - 1} seq[i]*4^i
  size_t x = 0;
  for(size_t i = 0; i < seq.size(); i++)
    x += base2int(seq[i])*(1 << 2*i);
  return x;
}

inline char
int2base(const int i){
  switch(i) {
    case 0  : return 'A';
    case 1  : return 'C';
    case 2  : return 'G';
    case 3  : return 'T';
    default : return 'N';
  }
}

inline string
int2string(const int input_int,
           const size_t string_length){
  int remainder = input_int;
  string i2s;
  i2s.resize(string_length, 'A');
  for(size_t i = string_length; i > 0; i--){
    int divisor = 1 << 2*(i - 1);
    i2s[i - 1] = int2base(remainder/divisor);
    remainder = remainder % divisor;
  }
  return i2s;
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

// wildcard = 'N'
int
LevenshteinWildcardMetric(const string &s1,
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
                           dist[i-1][j-1] + ((to_upper(s1[i-1]) == to_upper(s2[j-1])
                                              || s1[i-1] == 'N'
                                              || s2[j-1] == 'N') ? 0:1)));
    }
  }
  return dist[s1.size()][s2.size()];
}

// wildcard = 'N'
int
MismatchWildcardMetric(const string &s1,
                       const string &s2){
  if(s1.size() != s2.size())
    cerr << "Cannote determine mismatch distance of string of two different lengths" << endl;
  assert(s1.size() == s2.size());
  
  int dist = 0;
  for(size_t i = 0; i < s1.size(); i++)
    dist += (s1[i] == s2[i] || s1[i] == 'N' || s2[i] == 'N') ? 0:1;
  
  return dist;
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
               const string region_name,
               const size_t len_sgRNA,
               size_t &n_forward,
               size_t &n_rev_comp,
               vector<sgRNA> &possible_sgRNAs){

 // string region_seq_rev_comp = reverse_complement(region_seq);
  string rev_PAM = reverse_complement(PAM);

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
        possible_sgRNAs.push_back(sgRNA(region_seq.substr(i - len_sgRNA, len_sgRNA + PAM_len),
                                        region_name, false));
      }
    }
    for(size_t i = 0; i < region_seq.size() - len_sgRNA - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(test_seq == rev_PAM){
        possible_sgRNAs.push_back(sgRNA(region_seq.substr(i, len_sgRNA + PAM_len),
                                        region_name, true));
      }
    }
  }
  
  else{
    if(VERBOSE){
      cerr << "wildcard detected" << endl;
    }
    // size_t last_wildcard_pos = PAM.find_last_not_of("ACGT");
    // test for more than one wildcard
    // if(first_wildcard_pos != last_wildcard_pos){
    //   throw SMITHLABException("more than one wildcard found in PAM");
    // }
    
    for(size_t i = len_sgRNA; i < region_seq.size() - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(test_seq.find_first_of("N") == std::string::npos){
        if(wildcard_seq_match(test_seq.c_str(), PAM.c_str())){
          possible_sgRNAs.push_back(sgRNA(region_seq.substr(i - len_sgRNA, len_sgRNA + PAM_len),
                                          region_name, false));
          ++n_forward;
        }
      }
    }
    for(size_t i = 0; i < region_seq.size() - len_sgRNA - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(test_seq.find_first_of("N") == std::string::npos){
        if(wildcard_seq_match(test_seq.c_str(), rev_PAM.c_str())){
          possible_sgRNAs.push_back(sgRNA(region_seq.substr(i, len_sgRNA + PAM_len),
                                          region_name, true));
          ++n_rev_comp;
        }
      }
    }
  }
  
}


// idea from http://www.cplusplus.com/forum/general/14268/:
// build a hash table with key as sgRNA string and
// value as the occurence count
// once done building the hash table, iterate through
// and delete the sgRNAs that have count > 1
void
remove_duplicate_sgRNAs(const bool VERBOSE,
                        const size_t len_sgRNA,
                        vector<sgRNA> &possible_sgRNAs){
  if(VERBOSE){
    cerr << "starting number of sgRNAs: "
         << possible_sgRNAs.size() << endl;
  }
  unordered_map<string, size_t> sgRNA_counts;
  for(size_t i = 0; i < possible_sgRNAs.size(); i++){
    // need to cut out PAM
    string test_seq;
    // reverse complement: PAM is at start
    if(possible_sgRNAs[i].rev_comp){
      test_seq = possible_sgRNAs[i].seq.substr(3, len_sgRNA);
    }
    // not reverse complement, PAM is end
    else{
      test_seq = possible_sgRNAs[i].seq.substr(0, len_sgRNA);
    }
    unordered_map<string, size_t>::iterator it = sgRNA_counts.find(test_seq);
    // if test_seq is not in hash table, add it
    if(it == sgRNA_counts.end()){
      sgRNA_counts.insert(std::make_pair<string, size_t>(test_seq, 1));
    }
    // if test_seq is in hash table, update counts
    else{
      it->second++;
    }
  }
  
  // now we have the hash table with counts,
  // iterate through and delete the sgRNA with count > 1
  for(size_t i = 0; i < possible_sgRNAs.size(); i++){
    string test_seq;
    // reverse complement: PAM is at start
    if(possible_sgRNAs[i].rev_comp){
      test_seq = possible_sgRNAs[i].seq.substr(3, len_sgRNA);
    }
    // not reverse complement, PAM is end
    else{
      test_seq = possible_sgRNAs[i].seq.substr(0, len_sgRNA);
    }
    unordered_map<string, size_t>::iterator it = sgRNA_counts.find(test_seq);
    // if more than 1 copy, delete
    if(it->second > 1){
      possible_sgRNAs.erase(possible_sgRNAs.begin() + i);
    }
  }
  if(VERBOSE){
    cerr << "ending number of sgRNAs: "
         << possible_sgRNAs.size() << endl;
  }
  
}


// return true if we encounter a new chromosome
// else return false
bool
load_chrom(const bool VERBOSE,
           std::ifstream &in,
           string &chrom_seq,
           string &next_chrom_name){
  // read in fasta file iteratively
  string buffer;
  chrom_seq.clear();
  next_chrom_name.clear();
  while(getline(in, buffer)){
    std::istringstream is(buffer);
    if(buffer[0] == '>'){
      // new chromosome, return name
      next_chrom_name.assign(buffer.substr(1, buffer.length() - 1));
      // break from while loop by returning
      return true;
    }
    else{
      std::transform(buffer.begin(), buffer.end(), buffer.begin(), to_upper);
      chrom_seq.append(buffer);
    }
  }
  return false;
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
                           "", "-i input fasta file of proposed region");
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
    
   // const size_t PAM_len = PAM_seq.size();
    const string PAM_rev_comp = reverse_complement(PAM_seq);
    if(VERBOSE){
      cerr << "PAM = " << PAM_seq << endl;
    }
    vector<string> seqs, names;
    size_t n_seqs = read_fasta_batch(input_file_name, seqs, names);
    
    if(VERBOSE){
      cerr << "sequences read in: " << n_seqs << endl;
  
      for(size_t i = 0; i < min(seqs.size(), names.size()); i++){
        cerr << names[i] << endl;
        cerr << seqs[i] << endl;
      }
    }
    
    vector<sgRNA> possible_sgRNAs;
    size_t n_forward = 0;
    size_t n_rev_comp = 0;
    for(size_t i = 0; i < seqs.size(); i++)
      propose_sgRNAs(false, PAM_seq, seqs[i], names[i], len_sgRNA, n_forward, n_rev_comp, possible_sgRNAs);
 
    cerr << "# possible sgRNAs = " << possible_sgRNAs.size() << endl;
    /*
    if(VERBOSE){
      cerr << "proposed sgRNAs:" << endl;
      for(size_t i = 0; i < possible_sgRNAs.size(); i++){
        cerr << possible_sgRNAs[i] << endl;
      }
    }
     */
    
    cerr << "removing duplicates" << endl;
    remove_duplicate_sgRNAs(VERBOSE, len_sgRNA, possible_sgRNAs);
    // remove duplicates
    sort(possible_sgRNAs.begin(), possible_sgRNAs.end(), compare_sgRNA);
    // vector<sgRNA>::iterator unique_sgRNAs = unique(possible_sgRNAs.begin(), possible_sgRNAs.end(), equal_sgRNA_seq);
    // possible_sgRNAs.resize(std::distance(possible_sgRNAs.begin(), unique_sgRNAs));

    cerr << "# possible sgRNAs after removing duplicates = " << possible_sgRNAs.size() << endl;
    
    if(VERBOSE){
      cerr << "proposed sgRNAs:" << endl;
      for(size_t i = 0; i < possible_sgRNAs.size(); i++){
        cerr << possible_sgRNAs[i] << endl;
      }
    }
 
    std::ofstream of;
    if (!output_file_name.empty()) of.open(output_file_name.c_str());
    std::ostream out(output_file_name.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    for(size_t i = 0; i < possible_sgRNAs.size(); i++){
      out << ">" << possible_sgRNAs[i].target_name << endl;
      out << possible_sgRNAs[i].seq << endl;
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