#pragma once
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace KmerBruteForce
{
   using kmer_vec_type = std::vector<std::pair<std::string, unsigned int>>;

   kmer_vec_type
   findKmerFrequencies(
      const std::string& seq, 
      unsigned int kmer_len, 
      unsigned int num_kmer);
}