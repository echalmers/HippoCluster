#include "stdafx.h"
#include "ChineseWhispers.h"
#include <unordered_map>
#include <vector>
#include <fstream>

using namespace std;

namespace HippoClusterLibrary
{
	ChineseWhispers::ChineseWhispers(AdjacencyList* AdjList)
	{
		adjList = AdjList;
		numVertices = adjList->numVertices();

		clusterAssignments.clear();
		for (int i = 0; i < adjList->numVertices(); i++)
		{
			clusterAssignments.push_back(i);
		}
	}

	void ChineseWhispers::step()
	{
		// perform the operation on all nodes in a random order
		vector<int> nodeNumbers;
		for (int i = 0; i < numVertices; i++)
		{
			nodeNumbers.push_back(i);
		}

		for (int st = 0; st < numVertices; st++)
		{
			// select random vertex
			//int v = (rand() * (RAND_MAX + 1) + rand()) % numVertices;
			int index = (rand() * (RAND_MAX + 1) + rand()) % nodeNumbers.size();
			int v = nodeNumbers[index];
			nodeNumbers.erase(nodeNumbers.begin() + index);


			// get most popular class assignment among this vertex's neighbors
			vector<tuple<int, double, double>> neighbors = adjList->allNeighbors(v);
			unordered_map<int, double> counts;
			for (int i = 0; i < neighbors.size(); i++)
			{
				int thisNeighbor = std::get<0>(neighbors[i]);
				int neighborsCluster = clusterAssignments[thisNeighbor];
				double neighborsWeight = std::get<2>(neighbors[i]);
				if (counts.emplace(neighborsCluster, neighborsWeight).second == false)
				{
					counts[neighborsCluster] += neighborsWeight;
				}
			}
			// find out the total weighting of the most popular cluster
			//int mostPopularCluster = -1; 
			double mostPopularCount = -1;
			for (auto i = counts.begin(); i != counts.end(); i++)
			{
				if (i->second > mostPopularCount)
				{
					mostPopularCount = i->second;
					//mostPopularCluster = i->first;
				}
			}
			vector<int> tiedWinners;
			for (auto i = counts.begin(); i != counts.end(); i++)
			{
				if (i->second == mostPopularCount)
					tiedWinners.push_back(i->first);
			}
			clusterAssignments[v] = tiedWinners[rand() % tiedWinners.size()];//mostPopularCluster;
		}
	}

	int ChineseWhispers::getCluster(string vertexName)
	{
		int vertexNumber = adjList->getVertexNumber(vertexName);
		
		return clusterAssignments[vertexNumber];
	}

	void ChineseWhispers::getAllClusterAssignments(std::vector<int>& assignments)
	{
		assignments = std::vector<int>(clusterAssignments);
	}

	void ChineseWhispers::exportClusterAssignments(std::string filename)
	{
		std::ofstream file;
		file.open(filename);
		for (int i = 0; i < numVertices; i++)
		{
			std::string name = adjList->getVertexName(i);
			file << name << '\t' << clusterAssignments[i] << std::endl;
		}
		file.close();
	}

	ChineseWhispers::~ChineseWhispers()
	{
	}
}
