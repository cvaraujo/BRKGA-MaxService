#include <iostream>
#include "MSDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include <chrono>

int main(int argc, char* argv[]) {
  const unsigned p = 500;	// size of population
  const double pe = 0.15;		// fraction of population to be the elite-set
  const double pm = 0.10;		// fraction of population to be replaced by mutants
  const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
  const unsigned K = 3;		// number of independent populations
  const unsigned MAXT = 3;	// number of threads for parallel decoding
	
  MSDecoder decoder;			// initialize the decoder
  decoder.loadInstance(argv[1], argv[2]); // Load the instance
  srand(time(NULL));
  const unsigned n = decoder.getM();		// size of chromosomes
  const long unsigned rngSeed = 0;	// seed to the random number generator
  MTRand rng(rngSeed);				// initialize the random number generator
	
  // initialize the BRKGA-based heuristic
  BRKGA< MSDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
	
  unsigned generation = 0;	// current generation
  const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
  const unsigned X_NUMBER = 2;	// exchange top 2 best
  const unsigned MAX_GENS = 1000;	// run for 1000 gens
  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  double runtime = 0;
  
  
  do {
    algorithm.evolve();	// evolve the population for one generation
    cout << algorithm.getBestFitness() << " -> " << decoder.getIncumbent() << endl;
    //    getchar();
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    end = chrono::steady_clock::now();
    runtime = chrono::duration_cast<chrono::seconds>(end - start).count();
  } while (generation < MAX_GENS && runtime < 1800);

  cout << algorithm.getBestFitness() << " -> " << decoder.getIncumbent() << endl;
  
  ofstream output;
  output.open(argv[3]);

  output << "UB: " << decoder.getIncumbent() << "\n" << "Runtime: " << runtime << endl;
  output.close();
  std::cout << "Best solution found has objective value = "
	    << decoder.getIncumbent()  << std::endl;
	
  return 0;
}
