/*
* Created by Carlos - 07/09/2020
*/

#include "MSDecoder.h"

MSDecoder::MSDecoder() {}

MSDecoder::~MSDecoder() {}

void MSDecoder::loadInstance(string instance, string param) {
  int u, v;
  double delay, jitter, bandwidth, ldp, paramDelayToken, paramJitterToken, paramVariationToken, paramBandwidthToken;
  int delayInt, jitterInt, variationInt;
  string token;
  ifstream fileBoostGraph, fileParam;

  fileParam.open(param, fstream::in);

  while (!fileParam.eof()) {
    fileParam >> token;
    if (token == "Delay") {
      fileParam >> token;
      if (token == "variation") {
	fileParam >> token >> paramVariationToken;
	paramVariation = int(1e5 * paramVariationToken);
      } else {
	fileParam >> paramDelayToken;
	paramDelay = int(1e5 * paramDelayToken);
      }
    }
    if (token == "Jitter") {
      fileParam >> token >> paramJitterToken;
      paramJitter = int(1e6 * paramJitterToken);
    }
    if (token == "Bandwidth") {
      fileParam >> token >> paramBandwidthToken;
      paramBandwidth = int(paramBandwidthToken);
    }
  }

  fileBoostGraph.open(instance, fstream::in);

  while (!fileBoostGraph.eof()) {
    fileBoostGraph >> token;
    if (token == "Nodes") {
      fileBoostGraph >> n; n++;
      arcs = vector<vector<Arc>>(n, vector<Arc>());
      costs = vector<vector<pair<int, int>>>(n, vector<pair<int, int>>(n));
      graph = BoostGraph(n);
    }
    
    if (token == "Arcs") fileBoostGraph >> m; m += n;
    if (token == "E") {
      fileBoostGraph >> u >> v >> delay >> jitter >> bandwidth >> ldp;
      if (bandwidth >= paramBandwidth) {
	delayInt = int(1e5 * delay), jitterInt = int(1e6 * jitter);
	arcs[u].push_back(Arc(v, delayInt, jitterInt));
	costs[u][v] = costs[v][u] = make_pair(delayInt, jitterInt);
	add_edge(u, v, 0.0, graph);
      }
    }
    if (token == "Root") fileBoostGraph >> root;
    if (token == "T") fileBoostGraph >> u, terminals.push_back(u);
  }

  bool isTerminal;
  for (int i = 1; i < n; ++i) {
    isTerminal = false;
    if (i != root) {
      for (auto t : terminals) {
	if (i == t) {
	  isTerminal = true;
	  break;
	}
      }
      if (!isTerminal) nonTerminals.push_back(i), DuS.push_back(i);
    }
  }

  add_edge(root, 0, 0.0, graph);
  arcs[root].push_back(Arc(0, paramDelay, paramJitter));
  costs[root][0] = costs[0][root] = make_pair(paramDelay, paramJitter);

  for(auto k : DuS) {
    add_edge(0, k, 0.0, graph);
    arcs[0].push_back(Arc(k, paramDelay, paramJitter));
    costs[k][0] = costs[0][k] = make_pair(paramDelay, paramJitter);  
  }
  cout << "Load graph successfully" << endl;
}

// Revise this function
double MSDecoder::decode(const std::vector<double>& chromosome) {
  bool found;
  Edge ed;
  int edNum = 0, count = 0;
  vector<Edge> spanningTree;
  vector<int> delayVec = vector<int>(n), jitterVec = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n);

  // Update the arcs cost
  for (int i = 0; i < n; i++) {
    for (auto e : arcs[i]) {
      tie(ed, found) = edge(i, e.j, graph);
      if (found) boost::put(edge_weight_t(), graph, ed, chromosome[edNum++]);
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
    
    if (delayVec[k] > paramDelay || jitterVec[k] > paramJitter)
      losts++, notAttended[k] = true;
    else meanValue += delayVec[k];
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
	  if (delayVec[k] - delayVec[l] > paramVariation) {
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
  return count;
}

int MSDecoder::getM() const {
  return m;
}

int MSDecoder::getN() const {
  return n;
}


