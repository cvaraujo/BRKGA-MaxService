#include <iostream>
#include "MSDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include <chrono>
/*
void convertSolution(MSDecoder decoder, vector<double>& chromosome) {
  bool found;
  Edge ed;
  int edNum = 0, count = 0;
  vector<Edge> spanningTree = vector<Edge>();
  vector<int> delayVec = vector<int>(n), jitterVec = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n);
  int alphaD = 0, alphaJ = 0, alphaV = 0;

  // Update the arcs cost
  for (int i = 0; i < n; i++) {
    if (!removed[i])
      for (auto e : arcs[i]) {
	if (!removed[e.j]) {
	  tie(ed, found) = edge(i, e.j, graph);
	  if (found) boost::put(edge_weight_t(), graph, ed, chromosome[edNum++]);
	}
      }
  }
  // Compute the MST
  kruskal_minimum_spanning_tree(graph, back_inserter(spanningTree));

  // Compte the paths from the tree
  BoostGraph paths = BoostGraph(n);
  
  int i, j;
  vector<Vertex> pred = vector<Vertex>(n);
  vector<int> distance = vector<int>(n);

  for (std::vector<Edge>::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) {
    i = source(*ei, graph); j = target(*ei, graph);
    add_edge(i, j, costs[i][j].first, paths);
  }
  
  dijkstra_shortest_paths(paths, root, predecessor_map(make_iterator_property_map(pred.begin(), get(vertex_index, paths))).distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, paths))));

  int actual, losts = 0;
  double meanValue = 0.0;  
  
  for (auto k : terminals) {
    actual = k;
    delayVec[k] = distance[k];
    while(actual != root) {
      jitterVec[k] += costs[pred[actual]][actual].second;
      actual = pred[actual];
    }
    
    if (delayVec[k] > paramDelay || jitterVec[k] > paramJitter) {
      losts++, notAttended[k] = true;
      if (delayVec[k] > paramDelay) alhpaD += (paramDelay - delayVec[k]) / paramDelay;
      if (jitterVec[k] > paramJitter) alphaJ += (paramJitter - jitterVec[k]) / paramJitter;
    } else meanValue += delayVec[k];
  }
  
  meanValue = meanValue / double(terminals.size() - losts);
  int lessThanAvg = 0, greaterThanAvg = 0;

  for (auto k : terminals) {
    if (!notAttended[k]) {
      if (delayVec[k] <= meanValue) lessThanAvg++;
      else greaterThanAvg++;
    }
  }
  
  for (auto k : terminals) {
    if (!notAttended[k]) {
      for (auto l : terminals) {
	if (l != k && !notAttended[l]) {
	  if (abs(delayVec[k] - delayVec[l]) > paramVariation) {
	    alphaV += (abs(delayVec[k] - delayVec[l]) - paramVariation) / paramVariation;
	    if (lessThanAvg < greaterThanAvg) {
	      if (delayVec[k] < delayVec[l]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    } else {
	      if (delayVec[k] > delayVec[l]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    }
	  }
	}
      }
    }
  }

  count = 0;
  for (auto t : terminals) 
    if (notAttended[t]) count++;

}
*/
int main(int argc, char* argv[]) {
  const unsigned p = 100;	// size of population
  const double pe = 0.20;		// fraction of population to be the elite-set
  const double pm = 0.10;		// fraction of population to be replaced by mutants
  const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
  const unsigned K = 3;		// number of independent populations
  const unsigned MAXT = 2;	// number of threads for parallel decoding
	
  MSDecoder decoder;			// initialize the decoder
  decoder.loadInstance(argv[1], argv[2]); // Load the instance
	
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

    //int fitness = decoder.decodeFinal(algorithm.getBestChromosome());
    //  cout << algorithm.getBestFitness() << " -> " << fitness << endl;
    //    getchar();
    // cout << algorithm.getBestFitness() << std::endl;
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    end = chrono::steady_clock::now();
    runtime = chrono::duration_cast<chrono::seconds>(end - start).count();
  } while (generation < MAX_GENS && runtime < 1800);

  int fitness = decoder.decodeFinal(algorithm.getBestChromosome());
  cout << algorithm.getBestFitness() << " -> " << fitness << endl;

  
  ofstream output;
  output.open(argv[3]);

  output << "UB: " << fitness << "\n" << "Runtime: " << runtime << endl;
  output.close();
  std::cout << "Best solution found has objective value = "
	    << fitness  << std::endl;
	
  return 0;
}
