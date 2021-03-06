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
#include <iterator>
#include <cstdlib>



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
using std::atoi;


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
  sgRNA(const string s, const string n, const size_t k, const bool rc, const int p)
    {seq  = s; target_name = n; key = k; rev_comp = rc; pos = p; to_delete = false;}
  sgRNA(const string s, const string n, const bool rc, const int p, const string chr)
    {seq  = s; target_name = n; rev_comp = rc; pos = p; chrom = chr; to_delete = false;}


  
  string seq;
  string target_name;
  int key;
  bool rev_comp;
  bool to_delete;
  int pos;
  string chrom;
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
  if(sgrna.matches.size() > 0){
    the_stream << "matches: " << endl;
    for(size_t i = 0; i < sgrna.matches.size(); i++){
      the_stream << sgrna.matches[i] << endl;
    }
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

// look for string s1 in string s2
// wildcard = 'N'
bool
SeqMatchWithWildcard(const string &s1,
                     const string &s2){
  // string s1 has to be smaller than
  assert(s1.size() < s2.size());
  for(size_t i = 0; i < s2.size() - s1.size(); i++){
    if(s2[i] == s1[0] || s1[0] == 'N' || s2[i] == 'N'){
      // find match at first position, now look for full match
      bool match_test = true;
      // if a mismatch is found, set match_test = false and exit loop
      for(size_t j = 1; j < s1.size(); j++){
        if(!(s2[i + j] == s1[j] || s1[j] == 'N' || s2[i + j] == 'N')){
          // mismatch, stop checking
          match_test = false;
          break;
        }
      }
      if(match_test)
        return match_test;
    }
  }
  return false;
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
read_seqs(const string &input_file_name,
          vector<string> &seqs){
  std::ifstream in(input_file_name.c_str());
  if (!in) //if file doesn't open
    throw SMITHLABException("could not open input file: " + input_file_name);
  size_t line_count = 0ul;
  string buffer;
  string current_seq;
  while(getline(in, buffer)){
    ++line_count;
    std::istringstream is(buffer);
    seqs.push_back(buffer);
  }
  
  return line_count;
}


template<typename Out>
void split(const std::string &s, char delim, Out result) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}


void
split(const std::string &s, char delim, vector<string> &result) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  result = elems;
}


size_t
read_fasta_batch(const string &input_file_name,
                 vector<string> &seqs,
                 vector< vector<string> > &names){
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
      vector<string> name;
      split(buffer.substr(1, buffer.length() - 1), '\t', name);
      names.push_back(name);
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
               const vector<string> &region_name,
               const size_t len_sgRNA,
               size_t &n_forward,
               size_t &n_rev_comp,
               vector<sgRNA> &possible_sgRNAs){

 // string region_seq_rev_comp = reverse_complement(region_seq);
  string rev_PAM = reverse_complement(PAM);

  const size_t PAM_len = PAM.size();
  const string chrom = region_name[0];
  const int start = atoi(region_name[1].c_str());
  const int end = atoi(region_name[2].c_str());
  const string gene = region_name[3];
  
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
        possible_sgRNAs.push_back(sgRNA(region_seq.substr(i - len_sgRNA, len_sgRNA),
                                        gene, false,
                                        start + i - len_sgRNA,
                                        chrom));
      }
    }
    for(size_t i = 0; i < region_seq.size() - len_sgRNA - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(test_seq == rev_PAM){
        possible_sgRNAs.push_back(sgRNA(reverse_complement(region_seq.substr(i + PAM_len, len_sgRNA)),
                                        gene, true, start + i + PAM_len, chrom));
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
          possible_sgRNAs.push_back(sgRNA(region_seq.substr(i - len_sgRNA, len_sgRNA),
                                          gene, false,
                                          start + i - len_sgRNA,
                                          chrom));
          ++n_forward;
        }
      }
    }
    for(size_t i = 0; i < region_seq.size() - len_sgRNA - PAM_len; i++){
      string test_seq = region_seq.substr(i, PAM_len);
      if(test_seq.find_first_of("N") == std::string::npos){
        if(wildcard_seq_match(test_seq.c_str(), rev_PAM.c_str())){
          possible_sgRNAs.push_back(sgRNA(reverse_complement(region_seq.substr(i + PAM_len, len_sgRNA)),
                                          gene, true, start + i + PAM_len, chrom));
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
    string test_seq = possible_sgRNAs[i].seq;
    
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
    string test_seq = possible_sgRNAs[i].seq;
    
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



// need to remove constant tri nucleotides
// from the sgRNAs, but not the PAM
void
remove_trinucleotides(const bool VERBOSE,
                      const size_t len_sgRNA,
                      const size_t PAM_len,
                      vector<sgRNA> &possible_sgRNAs){
  const string aaa = "AAA";
  const string ccc = "CCC";
  const string ggg = "GGG";
  const string ttt = "TTT";
  for(vector<sgRNA>::iterator i = possible_sgRNAs.begin();
      i != possible_sgRNAs.end();){
    string sgRNAseq = (*i).seq;

    // now look for trinucleotides
    if(sgRNAseq.find(aaa) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else if(sgRNAseq.find(ccc) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else if(sgRNAseq.find(ggg) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else if(sgRNAseq.find(ttt) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else
      i++;
  }
}

// need to remove constant quad nucleotides
// from the sgRNAs, but not the PAM
void
remove_quadnucleotides(const bool VERBOSE,
                      const size_t len_sgRNA,
                      const size_t PAM_len,
                      vector<sgRNA> &possible_sgRNAs){
  const string aaaa = "AAAA";
  const string cccc = "CCCC";
  const string gggg = "GGGG";
  const string tttt = "TTTT";
  for(vector<sgRNA>::iterator i = possible_sgRNAs.begin();
      i != possible_sgRNAs.end();){
    string sgRNAseq = (*i).seq;
    
    // now look for trinucleotides
    if(sgRNAseq.find(aaaa) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else if(sgRNAseq.find(cccc) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else if(sgRNAseq.find(gggg) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else if(sgRNAseq.find(tttt) != std::string::npos){
      i = possible_sgRNAs.erase(i);
    }
    else
      i++;
  }
}

void
gc_correction(const bool VERBOSE,
              const size_t len_sgRNA,
              const size_t PAM_len,
              const double gc_content_lower_bound,
              const double gc_content_upper_bound,
              vector<sgRNA> & possible_sgRNAs){
  for(vector<sgRNA>::iterator i = possible_sgRNAs.begin();
      i != possible_sgRNAs.end();){
    string sgRNAseq = (*i).seq;

    double gc_count = 0;
    for(size_t j = 0; j < sgRNAseq.size(); j++){
      if(sgRNAseq[j] == 'C' || sgRNAseq[j] == 'G')
        gc_count++;
    }
    double gc_content = gc_count/sgRNAseq.size();
    if(gc_content < gc_content_lower_bound ||
       gc_content > gc_content_upper_bound){
      i = possible_sgRNAs.erase(i);
    }
    else
      i++;
  }
}




void
remove_enzyme_cut_seqs(const bool VERBOSE,
                       const size_t len_sgRNA,
                       const size_t PAM_len,
                       const vector<string> &enzymes,
                       vector<sgRNA> &possible_sgRNAs){
  // search through guides to find enzymes
  // due to the possibility of N's, need to do inexact matching
  for(vector<sgRNA>::iterator i = possible_sgRNAs.begin();
      i != possible_sgRNAs.end();){
    string sgRNAseq = (*i).seq;

    bool match_found = false;
    for(size_t j = 0; j < enzymes.size(); j++){
      if(SeqMatchWithWildcard(enzymes[j], sgRNAseq)){
        match_found = true;
        i = possible_sgRNAs.erase(i);
        break;
      }
    }
    // iterator is updated if match is found
    // otherwise need to set it manually
    if(!match_found)
      i++;
  }
}

void
add_g_for_U6promoter(vector<sgRNA> &possible_sgRNAs){
  for(size_t i = 0; i < possible_sgRNAs.size(); i++){
    string sgRNAseq = possible_sgRNAs[i].seq;
    if(sgRNAseq[0] != 'G'){
      sgRNAseq.insert(0, "G");
      possible_sgRNAs[i].seq = sgRNAseq;
    }
    // else do nothing
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


void
write_guides(const string output_file_name,
             const vector<sgRNA> &possible_sgRNAs){
  std::ofstream of;
  if (!output_file_name.empty()) of.open(output_file_name.c_str());
  std::ostream out(output_file_name.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  for(size_t i = 0; i < possible_sgRNAs.size(); i++){
    string seq;
    seq = possible_sgRNAs[i].seq;
      
    out << "> " << possible_sgRNAs[i].target_name << '\t'
        << possible_sgRNAs[i].chrom << '\t'
        << possible_sgRNAs[i].pos << '\t';
    if(possible_sgRNAs[i].rev_comp)
      out << "-";
    else
      out << "+";
    out<< endl;
    out << seq << endl;
  }
}


int
main(const int argc, const char **argv) {
  
  try {
    // option variables
//    size_t edit_dist = 2;
    size_t len_sgRNA = 20;
    string input_file_name;
    string genome_file_name;
    string output_file_name;
    string PAM_seq = "NGG";
    double gc_content_lower_bound = 0.1;
    double gc_content_upper_bound = 0.9;
    bool VERBOSE = false;
    bool REMOVE_DUPLICATES = false;
    bool REMOVE_TRINUCLEOTIDES = false;
    bool REMOVE_QUADNUCLEOTIDES = false;
    string enzyme_cut_seq_names;
    
    /********** GET COMMAND LINE ARGUMENTS  FOR C_CURVE ***********/
    OptionParser opt_parse(strip_path(argv[1]),
                           "", "-i input fasta file of proposed region");
    opt_parse.add_opt("output", 'o', "output file",
                      false , output_file_name);
//    opt_parse.add_opt("edit_distance", 'd', "maximum edit distance for matches"
//                      " (default: " + toa(edit_dist) + ")",
//                      false, edit_dist);
    opt_parse.add_opt("PAM", 'P', "PAM sequence (default: "
                      + PAM_seq + ")",
                      false, PAM_seq);
    opt_parse.add_opt("length", 's', "length of sgRNAs (default: "
                      + toa(len_sgRNA) + ")",
                      false, len_sgRNA);
    opt_parse.add_opt("gc_upper", 'u', "upper bound for GC content of sgRNAs"
                      " (default: " + toa(gc_content_upper_bound) + ")",
                      false, gc_content_upper_bound);
    opt_parse.add_opt("gc_lower", 'l', "lower bound for GC content of sgRNAs"
                      " (default: " + toa(gc_content_lower_bound) + ")",
                      false, gc_content_lower_bound);
    opt_parse.add_opt("input", 'i', "input file of the DNA sequence of the "
                      "regions to search in fasta format",
                      true, input_file_name);
    opt_parse.add_opt("VERBOSE", 'V', "verbose mode",
                      false, VERBOSE);
    opt_parse.add_opt("REMOVE_DUPLICATES", 'D', "remove duplicate sgRNAs",
                      false, REMOVE_DUPLICATES);
    opt_parse.add_opt("REMOVE_TRINUCLEOTIDES", 'T', "remove constant trinucleotides "
                      "(aaa, ccc, ggg, ttt)", false, REMOVE_TRINUCLEOTIDES);
    opt_parse.add_opt("REMOVE_QUADNUCLEOTIDES", 'Q', "remove constant quad nucleotides "
                      "(aaaa, cccc, gggg, tttt)", false, REMOVE_QUADNUCLEOTIDES);
    opt_parse.add_opt("cut_seqs", 'c', "filename of enzyme cutting sequences to remove",
                      false, enzyme_cut_seq_names);
    
    
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
    vector<string> seqs;
    vector< vector<string> > names;
    size_t n_seqs = read_fasta_batch(input_file_name, seqs, names);
    assert(names.size() == seqs.size());
    if(VERBOSE){
      cerr << "number of seqs read in = " << n_seqs << endl;
      cerr << "names = " << endl;
      for(size_t i = 0; i < names.size(); i++){
        for(size_t j = 0; j < names[i].size(); j++){
          cerr << names[i][j] << '\t';
        }
        cerr << endl;
//cerr << names[i][0] << '\t' << atoi(names[i][2].c_str()) << '\t'
//            << atoi(names[i][3].c_str()) << '\t' << names[i][1] << endl;
      }
      cerr << "proposing guides" << endl;
    }
    
    vector<sgRNA> possible_sgRNAs;
    size_t n_forward = 0;
    size_t n_rev_comp = 0;
    for(size_t i = 0; i < seqs.size(); i++){
      cerr << "seq " << i + 1 << endl;
      if(seqs[i].length() > len_sgRNA + PAM_seq.length()){
        propose_sgRNAs(false, PAM_seq, seqs[i], names[i], len_sgRNA, n_forward, n_rev_comp, possible_sgRNAs);
      }
    }
 
    cerr << "# possible sgRNAs = " << possible_sgRNAs.size() << endl;
    /*
    if(VERBOSE){
      cerr << "proposed sgRNAs:" << endl;
      for(size_t i = 0; i < possible_sgRNAs.size(); i++){
        cerr << possible_sgRNAs[i] << endl;
      }
    }
     */
    
    if(REMOVE_DUPLICATES){
      remove_duplicate_sgRNAs(VERBOSE, len_sgRNA, possible_sgRNAs);
      if(VERBOSE)
        cerr << "removing duplicates" << endl;

      if(VERBOSE)
        cerr << "# possible sgRNAs after removing duplicates = "
             << possible_sgRNAs.size() << endl;
    }
    if(REMOVE_TRINUCLEOTIDES){
      if(VERBOSE)
        cerr << "removing sgRNAs with constant trinucleotides" << endl;
      remove_trinucleotides(VERBOSE, len_sgRNA, PAM_seq.size(), possible_sgRNAs);
      if(VERBOSE)
        cerr << "# possible sgRNAs after removing trinucleotides = "
             << possible_sgRNAs.size() << endl;
    }
    else if(REMOVE_QUADNUCLEOTIDES){
      if(VERBOSE)
        cerr << "removing sgRNAs with constant quadnucleotides" << endl;
      remove_quadnucleotides(VERBOSE, len_sgRNA, PAM_seq.size(), possible_sgRNAs);
      if(VERBOSE)
        cerr << "# possible sgRNAs after removing trinucleotides = "
        << possible_sgRNAs.size() << endl;

    }
    if(!(enzyme_cut_seq_names.empty())){
      if(VERBOSE)
        cerr << "removing sgRNAs with the specified enzyme cutting sequences" << endl;
      vector<string> enzyme_cutting_seqs;
      read_seqs(enzyme_cut_seq_names, enzyme_cutting_seqs);
      if(VERBOSE){
        cerr << "enzyme cutting seqs = " << endl;
        for(size_t i = 0; i < enzyme_cutting_seqs.size(); i++)
          cerr << enzyme_cutting_seqs[i] << endl;
      }
      remove_enzyme_cut_seqs(VERBOSE, len_sgRNA, PAM_seq.size(),
                             enzyme_cutting_seqs, possible_sgRNAs);
      if(VERBOSE)
        cerr << "# possible sgRNAs after removing cutting sites = "
             << possible_sgRNAs.size() << endl;
      
    }
    if(VERBOSE)
      cerr << "removing guides based on GC content" << endl;
    
    gc_correction(VERBOSE, len_sgRNA, PAM_seq.size(),
                  gc_content_lower_bound, gc_content_upper_bound,
                  possible_sgRNAs);
    
    add_g_for_U6promoter(possible_sgRNAs);
    
    if(VERBOSE)
      cerr << "# possible sgRNAs after removing GC rich and poor guides = "
           << possible_sgRNAs.size() << endl;

    
    write_guides(output_file_name, possible_sgRNAs);

    
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