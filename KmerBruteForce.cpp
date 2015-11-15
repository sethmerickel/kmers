#include "KmerBruteForce.h"

#include <fstream>
#include <algorithm>
#include <exception>
#include <iterator>
#include <map>
#include <vector>

namespace KmerBruteForce
{
   kmer_vec_type
   findKmerFrequencies(
      const std::string& file_path,
      unsigned int kmer_len,
      unsigned int num_kmer)
   {
      std::ifstream fstrm(file_path);
      if (fstrm)
      {
         // ignore first line
         fstrm.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

         std::map<std::string, unsigned int> kmer_count_map;
         std::string line;
         while (!fstrm.eof())
         {
            std::getline(fstrm, line);
            auto n_kmers = line.size() - kmer_len + 1;
            
            // Loop through all the kmers and add them to the map.
            // Increment the count for the respective kmer
            for (size_t i = 0; i != n_kmers; ++i)
               ++kmer_count_map[line.substr(i, kmer_len)];

            // ignore 3 more lines
            fstrm.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            fstrm.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            fstrm.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
         }

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
      else
      {
         throw std::exception("Invalid file path");
      }
   }
}