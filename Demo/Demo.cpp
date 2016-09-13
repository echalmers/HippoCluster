// Demo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "../HippoClusterLibrary/AdjacencyList.h"
#include "../HippoClusterLibrary/HippoCluster.h"
#include <iostream>
#include <vector>
using namespace std;
using namespace HippoClusterLibrary;

int main()
{
	// read in the test graph
	AdjacencyList adjacencyList;
	adjacencyList.fromTSV("testGraph.tsv");

	// another way to create the graph is to add edge info one-by-one
	//AdjacencyList adjacencyList;
	//adjacencyList.addEdge("A", "B", 1);
	//adjacencyList.addEdge("B", "A", 1);
	//adjacencyList.addEdge("D", "A", 0.1);
	//etc...

	// list all vertex names and assigned id numbers
	cout << "vertex name	id" << endl;
	for (int i = 0; i < adjacencyList.numVertices(); i++)
	{
		cout << adjacencyList.getVertexName(i) << "\t\t" << i << endl;
	}

	// show some info about the graph
	cout << "number of vertices:  " << adjacencyList.numVertices() << endl;
	cout << "number of edges:     " << adjacencyList.numEdges() << endl;
	cout << "sum of edge weights: " << adjacencyList.sumWeights() << endl << endl;

	// find strongly connected components
	vector<vector<int>*> sccs;
	int maxIndex = adjacencyList.findSCCs(sccs);

	// keep only the biggest SCC
	cout << "reducing to biggest SCC" << endl;
	AdjacencyList mainSCC = adjacencyList.subgraph(sccs[maxIndex]);

	// list all vertex names and assigned id numbers in the main SCC
	cout << "vertex name	id" << endl;
	for (int i = 0; i < mainSCC.numVertices(); i++)
	{
		cout << mainSCC.getVertexName(i) << "\t\t" << i << endl;
	}

	// run a clustering algorithm on the full graph
	cout << endl << "running the new clustering algorith for 1000 iterations" << endl;
	HippoCluster hc(&adjacencyList);
	hc.setNumNodes(10);
	hc.setTrajectoryLength(adjacencyList.numVertices());
	for (int i = 0; i < 1000; i++)
	{
		hc.step();
		hc.dropUnusedNodes();
	}
	cout << "vertex name	cluster assignment" << endl;
	vector<int> assignments;
	hc.getAllClusterAssignments(assignments);
	for (int i = 0; i < assignments.size(); i++)
	{
		cout << adjacencyList.getVertexName(i) << "\t\t" << assignments[i] << endl;
	}

	// display some stats for the clustering
	ClusterStats stats = adjacencyList.getClusterStats(assignments);
	cout << "mean cluster size: " << stats.getClusterSizeMean() << " +/- " << stats.getClusterSizeSD() << endl;
	cout << "weighted average relative density: " << stats.getWeightedAvgRelativeDensity() << endl;

	
	getchar();
    return 0;
}

