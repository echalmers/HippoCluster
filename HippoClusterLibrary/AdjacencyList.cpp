#include "stdafx.h"
#include "AdjacencyList.h"
#include <string>
#include <fstream>
#include <time.h>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <ppl.h>

namespace HippoClusterLibrary
{
	// A class representing an adjacency list representation of a graph
	AdjacencyList::AdjacencyList()
	{
		srand(time(NULL));
	}

	// ____fromTSV____
	// load graph information from a tab-separated values file.
	// each row in the file defines one directed edge, and is formatted as follows:
	// vertexA \t vertexeB \t weight
	// where vertexA and vertexB are the (string) names of the vertices that the edge goes from and to, respectively.
	// weight is the numeric weight for the edge, and is optional.
	// The AdjacencyList object will store the string vertex names, but will assign each vertex an integer ID for internal use.
	void AdjacencyList::fromTSV(std::string filename)
	{
		adjList.clear();
		vertexNames.clear();
		vertexNumbers.clear();

		std::ifstream inFile;
		inFile.open(filename);

		// read in data from file
		if (inFile.is_open())
		{
			std::string line;
			while (getline(inFile, line))
			{
				// get the first two string values from the line
				std::vector<std::string> neighbors;
				int pos;
				while ((pos = line.find('\t')) != std::string::npos)
				{
					neighbors.push_back(line.substr(0, pos));
					line.erase(0, pos + 1);
				}
				// get the last numerical value
				double frequency;
				
				// if this text file doesn't include weights
				if (neighbors.size() == 1)
				{
					neighbors.push_back(line);
					frequency = 1;
				}
				else
					frequency = stod(line);

				addEdge(neighbors[0], neighbors[1], frequency, false);
			}
			inFile.close();
		}
		else
			throw("error opening file");

		// calculate edge probabilities
		for (int i = 0; i < adjList.size(); i++)
		{
			calculateNeighborProbabilities(i);
		}

	}

	// ____toTSV____
	// export the graph information to a tab-separated values file.
	void AdjacencyList::toTSV(std::string filename)
	{
		std::ofstream file;
		file.open(filename);
		for (int i = 0; i < adjList.size(); i++)
		{
			std::string start = vertexNames[i];
			for (int j = 0; j < adjList[i].size(); j++)
			{
				file << start << "\t" << vertexNames[std::get<0>(adjList[i][j])] << "\t1" << std::endl;
			}
		}
		file.close();
	}

	// ____addEdge____
	// add a directed adge to the graph
	// parameters:
	// vert1 - string name of the vertex the edge departs from (the string name will be stored, but this vertex will be given an integer ID for internal use)
	// vert2 - string name of the vertex the edge is incident to (the string name will be stored, but this vertex will be given an integer ID for internal use)
	// count - weight for this edge
	// updateProbabilities - whether to update probabilities associated with a transition out of vert1 after adding this edge info (default true)
	void AdjacencyList::addEdge(std::string vert1, std::string vert2, double count, bool updateProbabilities)
	{
		// add first vertex to vertex list if necessary
		if (vertexNumbers.emplace(vert1, vertexNames.size()).second == true)
		{
			vertexNames.push_back(vert1);
			adjList.push_back(std::vector<std::tuple<int, double, double>>());
		}

		// add second vertex to vertex list if necessary
		bool neighborExists = false;
		if (vertexNumbers.emplace(vert2, vertexNames.size()).second == true)
		{
			vertexNames.push_back(vert2);
			adjList.push_back(std::vector<std::tuple<int, double, double>>());
		}
		else // if this neighbor already existed in the adjacency list, modify the count or edge weight associated with it
		{
			for (int i = 0; i < adjList[vertexNumbers[vert1]].size(); i++)
			{
				if (vertexNumbers[vert2] == std::get<0>(adjList[vertexNumbers[vert1]][i]))
				{
					neighborExists = true;
					std::get<1>(adjList[vertexNumbers[vert1]][i]) += count;
					break;
				}
			}
		}

		// add to adjacency list
		if (!neighborExists)
		{
			std::tuple<int, double, double> t(vertexNumbers[vert2], count, 0);
			adjList[vertexNumbers[vert1]].push_back(t);
		}

		// recalculate edge probabilities
		if (updateProbabilities)
			calculateNeighborProbabilities(vertexNumbers[vert1]);
	}

	// ____calculateNeighborProbabilities(int startVertex)____
	// calculate probabilities associated with transitions from a specified vertex to its neighbors
	// parameters:
	// startVertex - the integer ID for the vertex. This ID can be retrieved using the getVertexNumber method.
	void AdjacencyList::calculateNeighborProbabilities(int startVertex)
	{
		double total = 0;
		for (int j = 0; j < adjList[startVertex].size(); j++)
		{
			total += std::get<1>(adjList[startVertex][j]);
		}
		for (int j = 0; j < adjList[startVertex].size(); j++)
		{
			std::get<2>(adjList[startVertex][j]) = std::get<1>(adjList[startVertex][j]) / total;
		}
	}

	// ____calculateNeighborProbabilities()____
	// calculate all probabilities associated with transitions between vertices, based on edge weights
	// this should be done after adding edge info to the graph, before actually using the graph.
	// it is not necessary if the graph was loaded using the fromTSV method
	void AdjacencyList::calculateNeighborProbabilities()
	{
		for (int i = 0; i < adjList.size(); i++)
		{
			calculateNeighborProbabilities(i);
		}
	}
	
	// ____randNeighbor____
	// probabilistically select a random neighbor of a specified vertex, based on the edge weights
	// parameters:
	//    vertex: the integer number of the source vertex. This number can be retrieved using the getVertexNumber method.
	// returns:
	//    the integer number of the selected neighbor.
	int AdjacencyList::randNeighbor(int vertex)
	{
		double p = (double)rand() / RAND_MAX;
		double total = 0;
		for (int i = 0; i < adjList[vertex].size(); i++)
		{
			total += std::get<2>(adjList[vertex][i]);
			if (total >= (p - epsilon))
			{
				return std::get<0>(adjList[vertex][i]);
			}
		}

		//std::cout << "error" << std::endl;
		return -1;
	}

	// ____allNeighbors____
	// returns all neighbors of a specified vertex (all vertices which can be reached from the vertex via a directed edge)
	// parameters:
	//    vertex: the integer number of the source vertex. This number can be retrieved using the getVertexNumber method.
	// returns:
	//     a vector of tuple<int, double, double>. The first element of each pair is the integer index of one neighbor of vertex.
	//     the second element is the edge weight between the two vertices
	//     the third element is the probability associated with a transition to that neighbor, given the weights.
	std::vector<std::tuple<int, double, double>> AdjacencyList::allNeighbors(int vertex)
	{
		std::vector<std::tuple<int, double, double>> neighbors;  
		for (int i = 0; i < adjList[vertex].size(); i++)
		{
			//std::pair<int, double> thisNeighbor(std::get<0>(adjList[vertex][i]), std::get<2>(adjList[vertex][i]));
			neighbors.push_back(adjList[vertex][i]);// thisNeighbor);
		}
		return neighbors;
	}

	// ____isNeighborOf____
	// returns all vertices that a specified vertex is a neighbor of (all vertices which can reach the specified vertex via a directed edge)
	// parameters:
	//    vertex: the integer number of the destination vertex. This number can be retrieved using the getVertexNumber method.
	// returns:
	//     a vector of tuple<int, double, double>. The first element of each pair is the integer ID of one vertex for which the specified vertex is a neighbor.
	//     the second element is the edge weight between the two vertices
	//     the third element is the probability associated with a transition to the specified neighbor, given edge weights.
	std::vector<std::tuple<int, double, double>> AdjacencyList::isNeighborOf(int vertex)
	{
		std::vector<std::tuple<int, double, double>> vertices;

		for (int v = 0; v < adjList.size(); v++)
		{
			for (int n = 0; n < adjList[v].size(); n++)
			{
				int sourceID = v;
				int destinationID = std::get<0>(adjList[v][n]);
				if (destinationID == vertex)
				{
					double weight = std::get<1>(adjList[v][n]);
					double prob = std::get<2>(adjList[v][n]);
					vertices.push_back(std::tuple<int, double, double>(sourceID, weight, prob));
				}
			}
		}
		return vertices;
	}


	// ____randPath____
	// generate a random path through the graph, starting at a specified vertex
	// parameters:
	//    startVertex - the integer index of the start vertex. This number can be retrieved using the getVertexNumber method.
	//    length - the length of the path
	//    path - reference to a vector<int> where the integer IDs of vertices along the path will be stored.
	//    pathVector - reference to a vector<double> with the same length as the number of vertices in the graph.
	//                 elements of this vector corresponding to the integer indices in 'path' will be set to 1.
	void AdjacencyList::randPath(int startVertex, int length, std::vector<int>& path, std::vector<double>& pathVector)
	{
		pathVector.resize(adjList.size(), 0);
		pathVector[startVertex] = 1;

		path.clear();
		path.push_back(startVertex);

		int currentVertex = startVertex;
		for (int i = 1; i < length; i++)
		{
			int nextVertex = randNeighbor(currentVertex);
			if (nextVertex == -1)
			{
				break;
			}

			path.push_back(nextVertex);
			pathVector[nextVertex] = 1;
			currentVertex = nextVertex;
		}
	}

	// ____getVertexName____
	// given the integer ID of a vertex, returns that vertex's string name
	std::string AdjacencyList::getVertexName(int vertexNumber)
	{
		return vertexNames[vertexNumber];
	}

	// ____vertexExists____
	// checks to see if a vertex with a given string name exists in the graph
	bool AdjacencyList::vertexExists(std::string vertexName)
	{
		return vertexNumbers.count(vertexName) > 0;
	}

	// ____getVertexNumber____
	// given the string name of a vertex, returns that vertex's integer ID
	int AdjacencyList::getVertexNumber(std::string vertexName)
	{
		return vertexNumbers[vertexName];
	}

	// ____numVertices____
	// returns the number of vertices in the graph
	int AdjacencyList::numVertices()
	{
		return adjList.size();
	}

	// ____numEdges____
	// returns the number of edges in the graph
	int AdjacencyList::numEdges()
	{
		int e = 0;
		for (int i = 0; i < adjList.size(); i++)
		{
			e += adjList[i].size();
		}
		return e;
	}

	// ___sumWeights____
	// returns the sum of all edge weights in the graph
	double AdjacencyList::sumWeights()
	{
		double sum = 0;
		for (int v = 0; v < adjList.size(); v++)
		{
			for (int n = 0; n < adjList[v].size(); n++)
			{
				sum += std::get<1>(adjList[v][n]);
			}
		}
		return sum;
	}

	// ____countAbsorbingVertices____
	// returns the number of absorbing vertices in the graph (vertices with no outgoing edges)
	int AdjacencyList::countAbsorbingVertices()
	{
		std::ofstream file;
		file.open("absorbingVertices.csv");

		int numAbsorb = 0;
		for (int i = 0; i < adjList.size(); i++)
		{
			if (adjList[i].size() == 0)
			{
				numAbsorb++;
				file << vertexNames[i] << std::endl;
			}
		}

		file.close();

		return numAbsorb;
	}

	// depreciated
	void AdjacencyList::countSubgraphs(std::unordered_map<int, int>& subgraphSizes, std::vector<int>& subgraphAssignments)
	{
		int currentSubgraphNumber = 0;
		subgraphAssignments.resize(adjList.size(), -1);

		std::queue<int> Q;
		
		for (int startV = 0; startV < adjList.size(); startV++)
		{
			std::cout << startV << " (" << currentSubgraphNumber << ")" << std::endl;
			if (startV>0 && subgraphAssignments[startV] != -1)
				continue;

			Q.push(startV);
			subgraphAssignments[startV] = currentSubgraphNumber++;

			while (Q.size() > 0)
			{
				int current = Q.front();
				//std::cout << "current: " << current << " (graph: " << subgraphNumber[current] << ")" << std::endl;
				Q.pop();

				for (int neighborIndex = 0; neighborIndex < adjList[current].size(); neighborIndex++)
				{
					int neighborNumber = std::get<0>(adjList[current][neighborIndex]);
					//std::cout << "   neighbor: " << neighborNumber << std::endl;
					if (subgraphAssignments[neighborNumber] == -1)
					{
						subgraphAssignments[neighborNumber] = subgraphAssignments[current];
						//std::cout << "   assign to graph: " << subgraphNumber[neighborNumber] << std::endl;
						Q.push(neighborNumber);
					}
					else if (subgraphAssignments[neighborNumber] != subgraphAssignments[current])
					{
						int mapFrom = subgraphAssignments[current];
						int mapTo = subgraphAssignments[neighborNumber];
						std::cout << "map all members of " << mapFrom << " to " << mapTo << std::endl;

						for (int i = 0; i < subgraphAssignments.size(); i++)
						{
							if (subgraphAssignments[i] == mapFrom)
								subgraphAssignments[i] = mapTo;
							if (subgraphAssignments[i] > min(mapFrom, mapTo))
								subgraphAssignments[i]--;
							//std::cout << i << ": " << subgraphNumber[i] << std::endl;
						}

					}
				}
			}
		}

		std::cout << "writing " << subgraphAssignments.size() << " subgraph assignments" << std::endl;
		std::ofstream file;
		file.open("graphAssignments.csv");
		for (int i = 0; i < subgraphAssignments.size(); i++)
		{
			file << subgraphAssignments[i] << std::endl;
		}
		file.close();

		// calculate subgraph sizes
		for (int i = 0; i < subgraphAssignments.size(); i++)
		{
			if (subgraphSizes.emplace(subgraphAssignments[i], 1).second == false)
			{
				subgraphSizes[subgraphAssignments[i]]++;
			}
		}
		file.open("subgraphSizes.csv");
		for (auto& i:subgraphSizes)
		{
			file << i.first << "," << i.second << std::endl;
		}
		file.close();

	}

	
	// ____findSCCs____
	// runs Tarjan's algorithm to identify all strongly-connected components in the graph
	// parameters:
	//    SCCs - reference to a vector of vectors that contain integer IDs of vertices in one strongly connected component
	// returns:
	//    the index in SCCs corresponding to the largest SCC
	int AdjacencyList::findSCCs(std::vector<std::vector<int>*>& SCCs)
	{
		std::stack<int> S;
		long index = 0;
		std::vector<int> indices(adjList.size(), -1);
		std::vector<int> lowLinks(adjList.size(), -1);
		std::vector<bool> onStack(adjList.size(), false);

		for (int v = 0; v < adjList.size(); v++)
		{
			//std::cout << index << std::endl;
			if (indices[v] == -1)
			{
				strongConnect(v, index, S, indices, lowLinks, onStack, SCCs);
			}
		}

		std::ofstream file;
		file.open("SCCsizes.txt");

		int maxSize = SCCs[0]->size(); int maxIndex = 0;
		for (int i = 1; i < SCCs.size(); i++)
		{
			file << SCCs[i]->size() << std::endl;

			if (SCCs[i]->size()>maxSize)
			{
				maxSize = SCCs[i]->size();
				maxIndex = i;
			}
		}
		file.close();

		return maxIndex;

	}

	// ____strongConnect____
	// private function for use in finding strongly connected components
	void AdjacencyList::strongConnect(int v, long& index, std::stack<int>& S, std::vector<int>& indices, std::vector<int>& lowLinks, std::vector<bool>& onStack, std::vector<std::vector<int>*>& SCCs)
	{
		// Set the depth index for v to the smallest unused index
		indices[v] = index;
		lowLinks[v] = index;
		index++;
		S.push(v);
		onStack[v] = true;

		// Consider successors of v
		for (int n = 0; n < adjList[v].size(); n++)
		{
			int w = std::get<0>(adjList[v][n]);
			if (indices[w] == -1)
			{
				// Successor w has not yet been visited; recurse on it
				strongConnect(w, index, S, indices, lowLinks, onStack, SCCs);
				lowLinks[v] = min(lowLinks[v], lowLinks[w]);
			}
			else if (onStack[w])
			{
				// Successor w is in stack S and hence in the current SCC
				lowLinks[v] = min(lowLinks[v], indices[w]);
			}
		}

		// If v is a root node, pop the stack and generate an SCC
		if (lowLinks[v] == indices[v])
		{
			SCCs.push_back(new std::vector<int>());
			int w = -1;
			do
			{
				w = S.top();
				S.pop();
				onStack[w] = false;
				SCCs[SCCs.size() - 1]->push_back(w);
			} while (w != v);
		}
	}

	// depreciated
	void AdjacencyList::removeSubgraphs(std::unordered_map<int, int>& subgraphSizes, std::vector<int>& subgraphAssignments, int minSize)
	{
		std::ofstream file;
		file.open("reducedGraph.tsv");

		for (int i = 0; i < adjList.size(); i++)
		{
			int thisAssignment = subgraphAssignments[i];
			if (subgraphSizes[thisAssignment] < minSize)
				continue;

			for (int n = 0; n < adjList[i].size(); n++)
			{
				std::string v1 = vertexNames[i];
				std::string v2 = vertexNames[std::get<0>(adjList[i][n])];
				file << v1 << "\t" << v2 << "\t" << std::get<1>(adjList[i][n]) << std::endl;
			}
		}

		file.close();
		fromTSV("reducedGraph.tsv");
	}

	// ____subgraph____
	// reduces the adjacency list to a specified subgraph
	// parameters: 
	//    verticesToKeep - pointer to a vector of integer vertex IDs to keep
	// returns:
	//    a new AdjacencyList object representing the subgraph
	AdjacencyList AdjacencyList::subgraph(std::vector<int>* verticesToKeep)
	{
		AdjacencyList newList;

		std::unordered_map<int, bool> keep;
		for (int i = 0; i < verticesToKeep->size(); i++)
		{
			keep[verticesToKeep->at(i)] = true;
		}

		for (int i = 0; i < verticesToKeep->size(); i++)
		{
			int v = verticesToKeep->at(i);
			for (int j = 0; j < adjList[v].size(); j++)
			{
				int n = std::get<0>(adjList[v][j]);
				if (keep[n])
					newList.addEdge(vertexNames[v], vertexNames[n], std::get<1>(adjList[v][j]), false);
			}
		}

		newList.calculateNeighborProbabilities();
		return newList;
	}

	// ____toSparseMatrix____
	// convert the adjacency list to a sparse adjacency matrix
	// parameters:
	//    mat - reference to an unordered map of unordered maps that will store the sparse matrix
	void AdjacencyList::toSparseMatrix(std::unordered_map<int, std::unordered_map<int, double>>& mat)
	{
		mat.clear();

		for (int i = 0; i < adjList.size(); i++)
		{
			for (int j = 0; j < adjList[i].size(); j++)
			{
				int neighbor = std::get<0>(adjList[i][j]);
				mat[i][neighbor] = std::get<2>(adjList[i][j]);
			}
		}
	}

	// ____getClusterStats____
	// parameters:
	//    clusterAssignments - reference to a vector<int>, in which the ith element indicates the cluster assignment for the vertex with integer ID i.
	// returns:
	//    a ClusterStats object containing some stats for the specified clustering
	ClusterStats AdjacencyList::getClusterStats(std::vector<int>& clusterAssignments)
	{
		CRITICAL_SECTION criticalSection;
		InitializeCriticalSection(&criticalSection);

		ClusterStats stats;

		// identify number of unique clusters
		std::unordered_set<int> uniqueClusters;
		std::unordered_map<int, int> renumberedAssignments;
		for (int i = 0; i < clusterAssignments.size(); i++)
		{
			if (uniqueClusters.emplace(clusterAssignments[i]).second == true)
			{
				renumberedAssignments[clusterAssignments[i]] = uniqueClusters.size() - 1;
			}

		}
		int numClusters = uniqueClusters.size();
		int numVertices = adjList.size();

		// compute internal and external degrees for each cluster
		std::vector<double> degInt(numClusters, 0);
		std::vector<double> degExt(numClusters, 0);
		std::vector<double> clusterCnts(numClusters, 0);

		concurrency::parallel_for(0, numVertices, [&](int v)
		//	for (int v = 0; v < numVertices; v++)
		{
			int thisRenumberedAssignment = renumberedAssignments[clusterAssignments[v]];
			clusterCnts[thisRenumberedAssignment]++;

			std::vector<std::tuple<int, double, double>> neighbors = allNeighbors(v);
			for (int n = 0; n < neighbors.size(); n++)
			{
				EnterCriticalSection(&criticalSection);
				if (clusterAssignments[v] == clusterAssignments[std::get<0>(neighbors[n])])
					degInt[thisRenumberedAssignment] += 1;// neighbors[n].second;
				else
					degExt[thisRenumberedAssignment] += 1;// neighbors[n].second;
				LeaveCriticalSection(&criticalSection);
			}
		}
		);
		
		

		// calculate average relative density
		stats.avgRelativeDensity = 0;
		stats.weightedAvgRelativeDensity = 0;
		int total = 0;
		for (int n = 0; n < numClusters; n++)
		{
			if (clusterCnts[n] > 0 && (degInt[n] + degExt[n])>0)
			{
				double temp = (degInt[n] / (degInt[n] + degExt[n]));
				stats.avgRelativeDensity += temp;
				stats.weightedAvgRelativeDensity += clusterCnts[n] * temp;
				total += clusterCnts[n];
			}
		}
		stats.avgRelativeDensity /= numClusters;
		stats.weightedAvgRelativeDensity /= total;

		// calculate cluster size stats
		stats.clusterSizeMax = 0;
		stats.clusterSizeSD = 0;
		stats.clusterSizeMean = 0;
		for (int i = 0; i < numClusters; i++)
		{
			stats.clusterSizeMean += clusterCnts[i];
			stats.clusterSizeMax = max(stats.clusterSizeMax, clusterCnts[i]);
		}
		stats.clusterSizeMean /= numClusters;

		stats.clusterSizeSD = 0;
		for (int i = 0; i < numClusters; i++)
		{
			stats.clusterSizeSD += (clusterCnts[i] - stats.clusterSizeMean) * (clusterCnts[i] - stats.clusterSizeMean);
		}
		stats.clusterSizeSD /= numClusters;
		double temp = sqrt(stats.clusterSizeSD);
		stats.clusterSizeSD = temp;

		DeleteCriticalSection(&criticalSection);

		return stats;
	}


	AdjacencyList::~AdjacencyList()
	{
	}






	ClusterStats::ClusterStats()
	{
	}

	HIPPOCLUSTERLIBRARY_API double ClusterStats::getAvgRelativeDensity()
	{
		return avgRelativeDensity;
	}

	HIPPOCLUSTERLIBRARY_API double ClusterStats::getWeightedAvgRelativeDensity()
	{
		return weightedAvgRelativeDensity;
	}

	HIPPOCLUSTERLIBRARY_API double ClusterStats::getClusterSizeMean()
	{
		return clusterSizeMean;
	}

	HIPPOCLUSTERLIBRARY_API double ClusterStats::getClusterSizeSD()
	{
		return clusterSizeSD;
	}

	HIPPOCLUSTERLIBRARY_API double ClusterStats::getClusterSizeMax()
	{
		return clusterSizeMax;
	}


	ClusterStats::~ClusterStats()
	{
	}
}
