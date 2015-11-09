#include <iostream>
#include <fstream>
#include <string>

#include "KmerBruteForce.h"

std::string usage =
"Usage: KmerTestDriver [fastq file]";

constexpr unsigned int kmer_len = 25;

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

int main(int argc, char* argv[])
{
   try
   {
      std::string seq; 
      if (argc == 1)
      {
         seq = "ABCABCABCABCABABABC";
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
   }
   catch (std::exception& e)
   {
      std::cerr << e.what() << std::endl;
   }

   return 0;
}