#include "KmerBruteForce.h"

#include <algorithm>
#include <exception>
#include <iterator>
#include <map>
#include <vector>

namespace KmerBruteForce
{
   kmer_vec_type
   findKmerFrequencies(
      const std::string& seq,
      unsigned int kmer_len,
      unsigned int num_kmer)
   {
      auto n = seq.size();
      if (kmer_len > n)
         throw std::exception("Invalid aguments to findKmerFrequencies");

      // Max number of kmers in the seq
      auto nkmers = n - kmer_len + 1;
      
      // Map to record each unique kmer and the number of times
      // it occurs in the sequence
      using kmer_map_type = std::map<std::string, unsigned int>;
      using map_value_type = kmer_map_type::value_type;
      kmer_map_type kmer_count_map;

      // Loop through all the kmers and add them to the map.
      // Increment the count for the respective kmer
      for (size_t i = 0; i != nkmers; ++i)      
         ++kmer_count_map[seq.substr(i, kmer_len)];

      // Copy map to vector for sorting by count
      using vec_value_type = kmer_vec_type::value_type;
      kmer_vec_type kmer_count_vec;
      kmer_count_vec.reserve(kmer_count_map.size());
      std::copy(
         begin(kmer_count_map),
         end(kmer_count_map), 
         std::back_inserter(kmer_count_vec));

      // functor for sorting counts in descending order instead 
      // of ascending order
      auto count_is_greater = 
         [](const vec_value_type& p1, const vec_value_type& p2) -> bool
      {
         return p1.second > p2.second;
      };

      // sort by counts
      std::sort(begin(kmer_count_vec), end(kmer_count_vec), count_is_greater);
      if (num_kmer < kmer_count_vec.size())
         return kmer_vec_type(begin(kmer_count_vec), begin(kmer_count_vec) + num_kmer);
      else
         return kmer_count_vec;
   }
}