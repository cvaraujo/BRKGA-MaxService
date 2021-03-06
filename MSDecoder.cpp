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

  BoostGraph auxGraph;
  
  while (!fileBoostGraph.eof()) {
    fileBoostGraph >> token;
    if (token == "Nodes") {
      fileBoostGraph >> n;
      arcs = vector<vector<Arc>>(n, vector<Arc>());
      costs = vector<vector<pair<int, int>>>(n, vector<pair<int, int>>(n));
      removed = vector<bool>(n);
      auxGraph = BoostGraph(n);
    }
    
    if (token == "Edges") fileBoostGraph >> m;
    if (token == "E") {
      fileBoostGraph >> u >> v >> delay >> jitter >> bandwidth >> ldp;
      u--, v--;
      if (bandwidth >= paramBandwidth) {
	delayInt = int(1e5 * delay), jitterInt = int(1e6 * jitter);
	arcs[u].push_back(Arc(v, delayInt, jitterInt));
	costs[u][v] = costs[v][u] = make_pair(delayInt, jitterInt);
	add_edge(u, v, 1.0, auxGraph);
      }
    }
    if (token == "Root") fileBoostGraph >> root, root--;
    if (token == "T") fileBoostGraph >> u, u--, terminals.push_back(u), DuS.push_back(u);
  }

  vector<Vertex> predecessors = vector<Vertex>(n);
  vector<int> distanceDelay = vector<int>(n);

  dijkstra_shortest_paths(auxGraph, root, predecessor_map(make_iterator_property_map(predecessors.begin(), get(vertex_index, auxGraph))).distance_map(make_iterator_property_map(distanceDelay.begin(), get(vertex_index, auxGraph))));
  
  bool isTerminal;
  for (int i = 0; i < n; ++i) {
    if (distanceDelay[i] >= numeric_limits<int>::max()) {
      removed[i] = true;
    } else {
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
  }

  graph = BoostGraph();
  for (int i = 0; i < n; i++)
    if (!removed[i])
      add_vertex(i, graph);
  
  for (int i = 0; i < n; i++){
    if (!removed[i]) {
      for (auto arc : arcs[i]) {
	if (!removed[arc.j])
	  add_edge(i, arc.j, 0.0, graph);
      }
    }
  }
  incumbent = terminals.size();
  cout << "Load graph successfully" << endl;
}

// Driver function to sort the vector elements 
// by second element of pairs 
bool sortbysec(const pair<int,int> &a, const pair<int,int> &b) { 
    return (a.second > b.second); 
} 

double MSDecoder::decode(const std::vector<double>& chromosome) {
  bool found;
  Edge ed;
  int i, j, edNum = 0, actual, jitter, delay, bestCand, obj;
  int count = 0;
  vector<Edge> spanningTree = vector<Edge>();
  vector<int> delayVec = vector<int>(n), jitterVec = vector<int>(n);
  vector<int> delayPaths = vector<int>(n), jitterPaths = vector<int>(n),
    predecessors = vector<int>(n), delayPathAux = vector<int>(n),
    jitterPathAux = vector<int>(n), predecessorsAux = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n), notAttendedAux = vector<bool>(n);
  vector<int> changed = vector<int>();
  double alpha;
  
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
  vector<int> distance = vector<int>(n);

  for (std::vector<Edge>::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) {
    i = source(*ei, graph); j = target(*ei, graph);
    add_edge(i, j, costs[i][j].first, paths);
  }
  
  dijkstra_shortest_paths(paths, root, predecessor_map(make_iterator_property_map(predecessors.begin(), get(vertex_index, paths))).distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, paths))));

  random_shuffle(terminals.begin(), terminals.end());

  obj = 0;
  for (auto k : DuS) {
    actual = k, jitter = 0, delay = 0;
    while (actual != root) {
      jitter += costs[predecessors[actual]][actual].second,
	delay += costs[predecessors[actual]][actual].first;
      actual = predecessors[actual];
    }
    delayPaths[k] = delay, jitterPaths[k] = jitter;
    // cout << k << " -> " << delay << " -> " << jitter << endl;
    if (delay > paramDelay || jitter > paramJitter) {
      notAttended[k] = true, obj++;
    }
  }
  // cout << paramDelay << " - " << paramJitter << endl;
  // getchar();

  int diffDelay, diffJitter;
  bool canMove;
  int selected = -1, losts;

  for (auto k : terminals) {
    if (!notAttended[k]) {
      selected = k;
      break;
    }
  }
  // cout << selected << endl;
  // getchar();
  if (selected == -1) {
    return numeric_limits<int>::max();
  } else {
    // get the values of this path
    for (auto k : terminals) {
      if (k != selected) {
	bestCand = -1, delay = paramDelay+1, jitter = paramJitter+1; 
            
	if (delayPaths[k] < (delayPaths[selected] - paramVariation) ||
	    delayPaths[k] > (delayPaths[selected] + paramVariation)) {
                
	  for (auto arc : arcs[k]) {
	    i = arc.j;
	    // Get the best candidate to move
	    if (i != predecessors[k] && k != predecessors[i]) {
	      if (delayPaths[i] + arc.delay >= (delayPaths[selected] - paramVariation) &&
		  delayPaths[i] + arc.delay <= (delayPaths[selected] + paramVariation) &&
		  delayPaths[i] + arc.delay <= paramDelay && 
		  jitterPaths[i] + arc.jitter <= paramJitter) {

		bestCand = i;
		delay = delayPaths[i] + arc.delay;
		jitter = jitterPaths[i] + arc.jitter;
		break;
	      }
	    }
	  }

	  if (bestCand != -1) {
	    // create the temporary vectors
	    notAttended[k] = false;
	    canMove = true;
	    losts = 0;
	    for (i = 0; i < n; i++) {
	      delayPathAux[i] = delayPaths[i], jitterPathAux[i] = jitterPaths[i], 
		predecessorsAux[i] = predecessors[i], notAttendedAux[i] = notAttended[i];
	    }

	    // Evaluate the move                    
	    diffDelay = delay - delayPathAux[k], diffJitter = jitter - jitterPathAux[k];
	    predecessorsAux[k] = bestCand;
	    delayPathAux[k] = delay, jitterPathAux[k] = jitter;

	    changed.erase(changed.begin(), changed.end());
	    changed.push_back(k);
	    while (!changed.empty()) {
	      actual = changed.back();
	      changed.pop_back();
	      for (int j : DuS) {
		if (j != actual && predecessorsAux[j] == actual) {
		  delayPathAux[j] += diffDelay, jitterPathAux[j] += diffJitter;
		  if (delayPathAux[j] > paramDelay || 
		      jitterPathAux[j] > paramJitter /*||
		      delayPathAux[j] < (delayPaths[selected] - paramVariation) ||
		      delayPathAux[j] > (delayPaths[selected] + paramVariation)*/) {
		    notAttendedAux[j] = true;
		    losts++;
		  } else notAttendedAux[j] = false;

		  if (losts >= 2 || (j == selected && notAttendedAux[j])) {
		    canMove = false;
		    changed.erase(changed.begin(), changed.end());
		    break;
		  } else changed.push_back(j);
		}
	      }
	    }

	    if (canMove) {
	      for (i = 0; i < n; i++) {
		delayPaths[i] = delayPathAux[i], jitterPaths[i] = jitterPathAux[i], 
		  predecessors[i] = predecessorsAux[i], notAttended[i] = notAttendedAux[i];
	      }
	    }
	  }  
	} 
      }
    }
  }


  vector<int> delays = vector<int>();  
  double auxMetric;
  vector<double> minMetric = vector<double>(3, numeric_limits<int>::max());
  for (auto k : terminals) {
    if (delayPaths[k] > paramDelay) {
      auxMetric = double(delayPaths[k] - paramDelay) / paramDelay;
      if (auxMetric < minMetric[0])
	minMetric[0] = auxMetric;
    } else if (jitterPaths[k] <= paramJitter) delays.push_back(delayPaths[k]);
 
    if (jitterPaths[k] > paramJitter) {
      auxMetric = double(jitterPaths[k] - paramJitter) / paramJitter;
      if (auxMetric < minMetric[1]) minMetric[1] = auxMetric;
    }

        for (auto l : terminals) {
      if (l != k && abs(delayPaths[k] - delayPaths[l]) > paramVariation) {
	auxMetric = double(abs(delayPaths[k] - delayPaths[l])) / paramVariation;
	if (auxMetric < minMetric[2]) minMetric[2] = auxMetric;
      }	
    }
  }

  sort(delays.begin(), delays.end());

  int actMax = 0, p2 = 0, maxi = 0, fd;
  for (i = 0; i < delays.size(); i++) {
    fd = delays[i];
    
    while (delays[p2] <= fd + paramVariation && p2 < delays.size()) p2++;

    if ((p2-i) > maxi) maxi = (p2-i);
  }

  int s = terminals.size();

  count = s - maxi;

  //  cout << "C: " << count << endl;
  /*int lessThanAvg = 0, greaterThanAvg = 0;
  double auxMetric;
  vector<double> minMetric = vector<double>(3, numeric_limits<int>::max());
   for (auto k : terminals) {
    for (auto l : terminals)
      if (k != l) 
	if (abs(delayPaths[k] - delayPaths[l]) > paramVariation) {
	  auxMetric = double(abs(delayPaths[k] - delayPaths[l]) - paramVariation) / paramVariation;
	  if (auxMetric < minMetric[2])
	    minMetric[2] = auxMetric;
	}

    if (k != selected && !notAttended[k]) {
      for (auto l : terminals) {
	if (l != k && l != selected && !notAttended[l]) {
	  if (delayPaths[k] - delayPaths[l] > paramVariation) {
	    if (lessThanAvg < greaterThanAvg) {
	      if (delayPaths[k] < delayPaths[selected]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    } else {
	      if (delayPaths[k] > delayPaths[selected]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    }
	  }
	}
      }
    }
    }*/

  /*  count = 0;
  for (auto t : terminals) 
    if (notAttended[t]) count++;

  int s = terminals.size();*/
  alpha = *min_element(minMetric.begin(), minMetric.end());

  if (count < incumbent){
    incumbent = count;
    cout << "Value: " << (count*s) + alpha << ", FO: " << incumbent << endl;
    /*
    ofstream file;
    file.open("solution.sol");
    for (i = 0; i < n; i++) {
      file << predecessors[i] << " " << i << endl;
    }
    file.close();*/
  }

  return ((count * s) + alpha);
}

int MSDecoder::getM() const {
  return m;
}

int MSDecoder::getN() const {
  return n;
}

int MSDecoder::getIncumbent() const {
  return incumbent;
}
/*
int MSDecoder::decodeFinal(const std::vector<double>& chromosome) {
  bool found;
  Edge ed;
  int edNum = 0, count = 0;
  vector<Edge> spanningTree = vector<Edge>();
  vector<int> delayPaths = vector<int>(n), jitterVec = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n);

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
    delayPaths[k] = distance[k];
    while(actual != root) {
      jitterVec[k] += costs[pred[actual]][actual].second;
      actual = pred[actual];
    }
    
    if (delayPaths[k] > paramDelay || jitterVec[k] > paramJitter)
      losts++, notAttended[k] = true;
    else meanValue += delayPaths[k];
  }
  
  meanValue = meanValue / double(terminals.size() - losts);
  int lessThanAvg = 0, greaterThanAvg = 0;

  for (auto k : terminals) {
    if (!notAttended[k]) {
      if (delayPaths[k] <= meanValue) lessThanAvg++;
      else greaterThanAvg++;
    }
  }
  
  for (auto k : terminals) {
    if (!notAttended[k]) {
      for (auto l : terminals) {
	if (l != k && !notAttended[l]) {
	  if (delayPaths[k] - delayPaths[l] > paramVariation) {
	    if (lessThanAvg < greaterThanAvg) {
	      if (delayPaths[k] < delayPaths[l]) {
		notAttended[k] = true;
		break;
	      } else notAttended[l] = true;
	    } else {
	      if (delayPaths[k] > delayPaths[l]) {
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

  for (auto k : terminals) {
    if (delayPaths[k] > paramDelay || jitterVec[k] > paramJitter) {
      if (!notAttended[k]) cout << "Erro no caminho" << endl;
    }
    if (!notAttended[k]) {
      for (auto l : terminals) {
	if (k != l) {
	  if (!notAttended[l]) {
	    if (delayPaths[k] - delayPaths[l] > paramVariation) cout << "Erro na variancia" << endl;
	  }
	}
      }
    }
  }
  return count;
}
*/
