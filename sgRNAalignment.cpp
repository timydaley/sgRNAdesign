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


struct sgRNA {
  sgRNA () {}
  sgRNA(const string s, const bool rc) {seq = s; rev_comp = rc;}
  sgRNA(const string s, const size_t k, const bool rc) {seq = s; key = k; rev_comp = rc;}
  
  string seq;
  int key;
  bool rev_comp;
  vector<MappedRead> matches;
};

std::ostream&
operator<<(std::ostream& the_stream, const sgRNA &sgrna){

  the_stream << sgrna.seq << endl;
  
  for(size_t i = 0; i < sgrna.matches.size(); i++){
    the_stream << sgrna.matches[i] << endl;
  }
  
  
  return the_stream;
}

typedef unordered_multimap<size_t, sgRNA> SeedHash;

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

//
inline size_t
updateHashValForward(const size_t prev_val,
                     const size_t seed_length,
                     const char start_char,
                     const char next_char){
  // next_val = (prev_val - (seq[0] * base^(seq.size() - 1))) * base + base2int(next_char) mod modulus
  // precompute = base^(seq.size() - 1)
  size_t next_val = prev_val - base2int(start_char);
  // make sure bit shifting will work correctly
  if(!(next_val % 4 == 0)){
    cerr << "possible issue with updating hash val" << endl;
    cerr << "prev_val = " << prev_val << endl;
    cerr << "start_char = " << start_char << endl;
    cerr << "next_val = " << next_val << endl;
  }
  assert(next_val % 4 == 0);
  // divide by 4 using bit shift
  next_val >>= 2;
  // add 4^(seed_length - 1)*base2int(next_char)
  next_val += (1 << 2*(seed_length - 1))*base2int(next_char);
  //  next_val = next_val % modulus; // if modulus = max size_t, comment this out
  return next_val;
}

/*
inline size_t
updateHashValReverseComp(const size_t prev_val,
                         const size_t seed_length,
                         const char prev_char,
                         const char next_char){
  const char prev_char_comp = complement(to_upper(prev_char));
  const char next_char_comp = complement(to_upper(next_char));

  // next_val = (prev_val - 4^(seed_length - 1)*base2int(start_char_comp))*4 + base2int(next_char_comp)
  size_t next_val = prev_val - (1 << 2*(seed_length - 1))*base2int(prev_char_comp);
  
  // multiply by 4 using bit shift
  next_val <<= 2;
  // add next char
  next_val += base2int(next_char_comp);
  //  next_val = next_val % modulus; // if modulus = max size_t, comment this out
  return next_val;
}
 */

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
        possible_sgRNAs.push_back(sgRNA(region_seq.substr(i - len_sgRNA, len_sgRNA + PAM_len), false));
      }
    }
    for(size_t i = 0; i < region_seq.size() - len_sgRNA - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(test_seq == rev_PAM){
        possible_sgRNAs.push_back(sgRNA(region_seq.substr(i, len_sgRNA + PAM_len), true));
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
      if(wildcard_seq_match(test_seq.c_str(), PAM.c_str())){
        possible_sgRNAs.push_back(sgRNA(region_seq.substr(i - len_sgRNA, len_sgRNA + PAM_len), false));
        ++n_forward;
      }
    }
    for(size_t i = 0; i < region_seq.size() - len_sgRNA - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(wildcard_seq_match(test_seq.c_str(), rev_PAM.c_str())){
        possible_sgRNAs.push_back(sgRNA(region_seq.substr(i, len_sgRNA + PAM_len), true));
        ++n_rev_comp;
      }
    }
  }
}


void
build_seed_hash(const bool VERBOSE,
                const size_t seed_length,
                const string PAM,
                const vector<sgRNA> &possible_sgRNAs,
                SeedHash &seed_hash){
  string rev_PAM = reverse_complement(PAM);
  
  for(size_t i = 0; i < possible_sgRNAs.size(); i++){
    string seed_seq;
    size_t key;
    if(possible_sgRNAs[i].rev_comp){
      seed_seq = possible_sgRNAs[i].seq.substr(3, seed_length);
    }
    else{
      seed_seq = possible_sgRNAs[i].seq.substr(possible_sgRNAs[i].seq.size()
                                              - PAM.size() - seed_length,
                                               seed_length);

    }
    key = string2int(seed_seq);
    if(VERBOSE){
      cerr << "seed = " << seed_seq << endl;
      cerr << "key  = " << key << endl;
    }
    seed_hash.insert(std::make_pair<size_t, sgRNA>(key, possible_sgRNAs[i]));
  }
}

/*
void
remove_duplicate_sgRNAs(const bool VERBOSE,
                        const size_t len_sgRNA,
                        const size_t seed_length,
                        SeedHash &seed_hash){
  // need to iterate over unique keys
  // upper_bound is the position of the next unique key
  SeedHash::iterator upper_bound = seed_hash.begin();
  size_t hash_count = 0;
  size_t n_entries = 0;
  
  for(SeedHash::iterator hashit = seed_hash.begin();
      hashit != seed_hash.end(); hashit = upper_bound){
    //cerr << "hash load = " << seed_hash.count(hashit->first) << endl;
    if(seed_hash.count(hashit->first) > 1){
      // more than one entry with the same key hashit->first
      vector<string> duplicate_hashes;
      std::pair<SeedHash::iterator, SeedHash::iterator>
        duplicate_iterator = seed_hash.equal_range(hashit->first);
      for(SeedHash::iterator it = duplicate_iterator.first; it != duplicate_iterator.second; it++){
        duplicate_hashes.push_back(it->second.seq.substr(0, len_sgRNA));
        n_entries++;
      }
      //cerr << "number of duplicate hashes = " << duplicate_hashes.size() << endl;
      upper_bound = duplicate_iterator.second;
    }
    else{
      upper_bound++;
      n_entries++;
    }
    hash_count++;
  }
  //cerr << "went through " << hash_count << " keys and " << n_entries << " entries" << endl;
}
*/

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


void
add_match(const bool VERBOSE,
          const string chrom_name,
          const int edit_dist,
          const size_t start_pos,
          const size_t end_pos,
          const char strand,
          SeedHash::iterator &it){
  MappedRead match;
  match.seq = it->second.seq;
  match.scr = edit_dist;
  GenomicRegion gr = GenomicRegion(chrom_name,
                                   start_pos, end_pos);
  gr.set_strand(strand);
  match.r = gr;
  it->second.matches.push_back(match);
}

int
main(const int argc, const char **argv) {
  
  try {
    // option variables
    size_t edit_dist = 2;
    size_t len_sgRNA = 20;
    size_t seed_length = 8;
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
    opt_parse.add_opt("seed_length", 's', "length of seed region distal to PAM"
                      " (default: " + toa(seed_length) + ")",
                      false, seed_length);
    opt_parse.add_opt("input", 'i', "input file of the DNA sequence of the "
                      "regions to search in fasta format",
                      true, input_file_name);
    opt_parse.add_opt("genome", 'g', "genome file in fasta format",
                      true, genome_file_name);
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
    
    const size_t PAM_len = PAM_seq.size();
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
      propose_sgRNAs(false, PAM_seq, seqs[i], len_sgRNA, n_forward, n_rev_comp, possible_sgRNAs);
 
    cerr << "# possible sgRNAs = " << possible_sgRNAs.size() << endl;
    if(VERBOSE){
      cerr << "proposed sgRNAs:" << endl;
      for(size_t i = 0; i < possible_sgRNAs.size(); i++){
        cerr << possible_sgRNAs[i] << endl;
      }
    }

    
    SeedHash seed_hash;
    build_seed_hash(false, seed_length, PAM_seq, possible_sgRNAs, seed_hash);
    if(VERBOSE){
      cerr << "seed hash table size = " << seed_hash.size() << endl;
    }
    
    // do we need to remove duplicate sgRNAs?
  //  remove_duplicate_sgRNAs(VERBOSE, len_sgRNA, seed_length, seed_hash);
    
    vector<string> chroms;
    vector<string> chrom_names;
    std::ifstream in(genome_file_name.c_str());
    read_fasta_batch(genome_file_name, chroms, chrom_names);
    if(VERBOSE){
      cerr << "name" << '\t' << "size" << endl;
      for(size_t i = 0; i < chroms.size(); i++)
        cerr << chrom_names[i] << '\t' << chroms[i].length() << endl;
    }
    assert(chroms.size() == chrom_names.size());

    size_t n_forward_hits = 0;
    size_t n_rev_comp_hits = 0;

    // loop over chroms
    for(size_t i = 0; i < chroms.size(); i++){
      size_t iter = chroms[i].find_first_not_of("Nn");
      iter += PAM_len;
      size_t hash_val = string2int(chroms[i].substr(iter, seed_length));
      if(VERBOSE){
        cerr << "iter = " << iter << endl;
        cerr << "hash_val = " << hash_val << endl;
        cerr << "seq = " << chroms[i].substr(iter, seed_length) << endl;
      }
      
      do{

        if(seed_hash.find(hash_val) != seed_hash.end()){
          if(VERBOSE){
            cerr << "match for " << hash_val << " found at " << iter << endl;
            cerr << "seq = " << chroms[i].substr(iter, seed_length) << endl;
            cerr << "int2string = " << int2string(hash_val, seed_length) << endl;
            cerr << "surrounding seq = " << chroms[i].substr(iter - PAM_len, seed_length + 2*PAM_len) << endl;
          }
          // test PAM for forward match
          if( (iter >= len_sgRNA) &&
              (MismatchWildcardMetric(chroms[i].substr(iter + seed_length, PAM_len),
                                    PAM_seq) == 0)){
            n_forward_hits++;
            
            std::pair<SeedHash::iterator, SeedHash::iterator>
              matches = seed_hash.equal_range(hash_val);
            for(SeedHash::iterator it = matches.first;
                it != matches.second; it++){
              string test_seq = chroms[i].substr(iter - len_sgRNA + seed_length, len_sgRNA);
              cerr << "proposed sgRNA = " << it->second.seq << endl;
              cerr << "full sequence  = " << chroms[i].substr(iter - len_sgRNA + seed_length, len_sgRNA + PAM_len) << endl;
              int d = LevenshteinWildcardMetric(test_seq,
                                                it->second.seq.substr(0, len_sgRNA));
              cerr << "edit distance  = " << d << endl;
              if(d <= edit_dist){
                // remove sgRNA if there is more than one match
                if(it->second.matches.size() > 1){
                  seed_hash.erase(it);
                }
                else{
                  cerr << "matches before:" << endl;
                  for(size_t j = 0; j < it->second.matches.size(); j++){
                    cerr << it->second.matches[i] << endl;
                  }
                  cerr << endl;
                  add_match(VERBOSE, chrom_names[i], d, iter + seed_length - len_sgRNA,
                            iter + PAM_len + seed_length, '-', it);
                  cerr << "matches after:" << endl;
                  for(size_t j = 0; j < it->second.matches.size(); j++){
                    cerr << it->second.matches[i] << endl;
                  }
                  cerr << endl;
                }
              }
              // if d > edit_dist, keep sgRNA
              /*
              else{
                cerr << "matches before:" << endl;
                for(size_t j = 0; j < it->second.matches.size(); j++){
                  cerr << it->second.matches[i] << endl;
                }
                cerr << endl;
                add_match(VERBOSE, chrom_names[i], d, iter + seed_length - len_sgRNA,
                          iter + PAM_len + seed_length, '-', it);
                cerr << "matches after:" << endl;
                for(size_t j = 0; j < it->second.matches.size(); j++){
                  cerr << it->second.matches[i] << endl;
                }
                cerr << endl;
              }
               */
            }
          }
        // test PAM for reverse match
          if( (iter > PAM_len) &&
               (MismatchWildcardMetric(chroms[i].substr(iter - PAM_len, PAM_len),
                                       PAM_rev_comp) == 0)){
            if(VERBOSE){
              cerr << "reverse complement match for " << hash_val << " found at " << iter << endl;
              cerr << "seq = " << chroms[i].substr(iter, seed_length) << endl;
            }
            n_rev_comp_hits++;
            std::pair<SeedHash::iterator, SeedHash::iterator>
            matches = seed_hash.equal_range(hash_val);
            for(SeedHash::iterator it = matches.first;
                it != matches.second; it++){
              string test_seq = reverse_complement(chroms[i].substr(iter + PAM_len, len_sgRNA));
              cerr << "proposed sgRNA = " << it->second.seq << endl;
              cerr << "full seq =       " << chroms[i].substr(iter - PAM_len, len_sgRNA + PAM_len) << endl;
              
              int d = LevenshteinWildcardMetric(test_seq,
                                                it->second.seq.substr(0, len_sgRNA));
              cerr << "edit distance  = " << d << endl;
              if(d <= edit_dist){
                // remove sgRNA if there is more than one match
                if(it->second.matches.size() > 1){
                  seed_hash.erase(it);
                }
                else{
                  if(VERBOSE){
                    cerr << "matches before:" << endl;
                    for(size_t j = 0; j < it->second.matches.size(); j++){
                      cerr << it->second.matches[i] << endl;
                    }
                    cerr << endl;
                  }
                  add_match(VERBOSE, chrom_names[i], d, iter - PAM_len,
                            iter + len_sgRNA, '-', it);
                  if(VERBOSE){
                    cerr << "matches after:" << endl;
                    for(size_t j = 0; j < it->second.matches.size(); j++){
                      cerr << it->second.matches[i] << endl;
                    }
                    cerr << endl;
                  }
                }
              }
              // if d > edit_dist, keep sgRNA
              /*
              else{
                if(VERBOSE){
                  cerr << "matches before:" << endl;
                  for(size_t j = 0; j < it->second.matches.size(); j++){
                    cerr << it->second.matches[i] << endl;
                  }
                  cerr << endl;
                }
                add_match(VERBOSE, chrom_names[i], d, iter - PAM_len,
                          iter + len_sgRNA, '-', it);
                if(VERBOSE){
                  cerr << "matches after:" << endl;
                  for(size_t j = 0; j < it->second.matches.size(); j++){
                    cerr << it->second.matches[i] << endl;
                  }
                  cerr << endl;
              }
               */

            }
          }
        }

        // update hash vals for RabinKarp after checking
        hash_val =
        updateHashValForward(hash_val, seed_length,
                             chroms[i][iter],
                             chroms[i][iter + seed_length]);
        /*
        hash_val =
        updateHashValReverseComp(hash_val, seed_length,
                                 chroms[i][iter + PAM_len],
                                 chroms[i][iter + PAM_len + seed_length]);
         */
        iter++;

      //  cerr << "iter = " << iter << endl;
        // what to do when you reach an N?
      }while(iter < chroms[i].length() - len_sgRNA - PAM_len && chroms[i][iter + len_sgRNA] != 'N');
      cerr << "number of bases iterated over = " << iter << endl;

    }
    
    cerr << "the number of proposed sgRNAs was " << possible_sgRNAs.size() << endl;
    cerr << "the number of filtered sgRNAs is  " << seed_hash.size() << endl;
    cerr << "number of forward sgRNAs = " << n_forward << endl;
    cerr << "number of rev comp sgRNAs = " << n_rev_comp << endl;
    cerr << "number of foward hits = " << n_forward_hits << endl;
    cerr << "number of backward hits = " << n_rev_comp_hits << endl;
    cerr << "hash table:" << endl;
 /*
    for(SeedHash::iterator i = seed_hash.begin(); i != seed_hash.end(); i++){
      cerr << i->first << '\t' << i->second << endl;
    }
*/
    
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