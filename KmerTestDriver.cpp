#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>

#include "KmerBruteForce.h"

std::string usage =
"Usage: KmerTestDriver FASTQ_FILENAME KMER_LENGTH NUMBER_OF_KMER";

unsigned int kmer_len = 30;
unsigned int num_kmers = 25;

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

// Functor for printing std::pair<T1, T2> 
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
         // Run with no arguments so just use a simple test case
         seq = "123456789123456789123456789123" \
            "123456789123456789123456789123" \
            "123456789123456789123456789123" \
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" \
            "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB";
      }
      else if (argc == 4)
      {         
         // parse fastq file
         seq = getStringFromFile(argv[1]);

         // second argument is the kmer length
         std::stringstream ss_len(argv[2]);
         if (!(ss_len >> kmer_len)) throw std::exception("Bad input args");

         // third argument is the number of kmers to return
         std::stringstream ss_num(argv[3]);
         if (!(ss_num >> num_kmers)) throw std::exception("Bad input args");
      }
      else
      {
         // wrong number of arguments
         throw std::exception(usage.c_str());
      }

      // Find the most frequently occuring kmers
      auto kmers = KmerBruteForce::findKmerFrequencies(seq, kmer_len, num_kmers);

      // Do some simple checks on the output for testing
      if (kmers.size() > num_kmers)
         throw std::exception("kmer algorithm failed.  Bug!");

      if (!kmers.empty() && kmers[0].first.size() != kmer_len)
         throw std::exception("kmer algorithm failed.  Bug!");
      
      // Print the results to stdout
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