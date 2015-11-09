#include "KmerBruteForce.h"

#include <exception>

namespace KmerBruteForce
{
   std::vector<std::pair<unsigned int, std::string>>
   findKmerFrequencies(const std::string& seq, unsigned int kmer_len)
   {
      throw std::exception("Function not implemented yet");

      auto l = seq.size();
      if (kmer_len > l)
         throw std::exception("Invalid aguments to findKmerFrequencies");

      auto nkmers = l - kmer_len + 1;
   }
}