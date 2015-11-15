#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>

#include "KmerBruteForce.h"

std::string usage =
"Usage: KmerTestDriver [fastq file]";

constexpr unsigned int kmer_len = 4;

std::string getStringFromFile(const char* file)
{
   std::ifstream fstrm(file);
   if (fstrm)
   {
      // ignore first line
      fstrm.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::string second_line;
      std::getline(fstrm, second_line);
      return second_line;
   }
   else
   {
      throw std::exception("Invalid file path");
   }
}

template <typename T>
struct pairPrinter;

template <typename T1, typename T2>
struct pairPrinter<std::pair<T1, T2>>
{
   pairPrinter(std::ostream& strm) : m_strm(strm) {}

   void operator()(const std::pair<T1, T2>& pair)
   {
      m_strm << '[' << pair.first << ", " << pair.second << "]\n";
   }

   std::ostream& m_strm;
};

int main(int argc, char* argv[])
{
   try
   {
      std::string seq; 
      if (argc == 1)
      {
         seq = "AAAAQWERKFJKTAAAAQEWRFLKJAAAA;LJF;LJVAAAA;LKFJA;LKAAAA";
      }
      else if (argc == 2)
      {         
         seq = getStringFromFile(argv[1]);
      }
      else
      {
         std::cout << usage << std::endl;
      }

      auto kmers = KmerBruteForce::findKmerFrequencies(seq, kmer_len);

      using value_type = KmerBruteForce::kmer_vec_type::value_type;
      std::for_each(begin(kmers), end(kmers), pairPrinter<value_type>{std::cout});
   }
   catch (std::exception& e)
   {
      std::cerr << e.what() << std::endl;
   }
   catch (...)
   {
      std::cerr << "Something horrible happened" << std::endl;
   }

   return 0;
}