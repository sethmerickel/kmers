#pragma once
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace KmerBruteForce
{
   std::vector<std::pair<unsigned int, std::string>>
   findKmerFrequencies(const std::string& seq, unsigned int kmer_len);
}