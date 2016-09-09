// TestClient.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "../HippoClusterLibrary/AdjacencyList.h"
#include "../HippoClusterLibrary/HippoCluster.h"
#include "../HippoClusterLibrary/ChineseWhispers.h"
#include "../HippoClusterLibrary/MCL.h"
#include <iostream>
#include <fstream>
#include <stack>
#include <unordered_map>
#include <ctime>
#include <algorithm>

using namespace HippoClusterLibrary;
using namespace std;

int main()
{

	AdjacencyList al;//
	al.fromCSV("C:\\Users\\Eric\\Google Drive\\Lethbridge Projects\\Fuzzy Place Field Test\\SimpleTestGraph.tsv");
	cout << "vertices:  " << al.numVertices() << endl;
	cout << "edges:     " << al.numEdges() << endl;
	cout << "abs. Vert: " << al.countAbsorbingVertices() << endl;

	unordered_map<int, unordered_map<int, double>> m;
	al.toSparseMatrix(m);

	for (int i = 0; i < al.numVertices(); i++)
	{
		for (int j = 0; j < al.numVertices(); j++)
		{
			cout << round(100*m[i][j])/100 << "\t";
		}
		cout << endl;
	}
	cout << endl << endl << endl;

	MCL mcl(&al);
	for (int i = 0; i < 50; i++)
	{
		mcl.step();
	}
	
	for (int i = 0; i < al.numVertices(); i++)
	{
		for (int j = 0; j < al.numVertices(); j++)
		{
			cout << round(100*mcl.T[i][j])/100 << "\t";
		}
		cout << endl;
	}

	vector<int> assign;
	mcl.getAllClusterAssignments(assign);

	for (int i = 0; i < assign.size(); i++)
	{
		cout << assign[i] << endl;
	}
	cout << "A: " << mcl.getCluster("C");

	//ChineseWhispers cw(&al);
	////HippoCluster cw(&al);
	////cw.setNumNodes(al.numVertices());
	////cw.setTrajectoryLength(9);
	//
	//for (int i = 0; i <= 100; i++)
	//{
	//	cw.step();
	//	
	//	vector<int> clusters;
	//	cw.getAllClusterAssignments(clusters);
	//	ClusterStats stats = al.getClusterStats(clusters);

	//	cout << stats.getAvgRelativeDensity() << endl;

	//	cout << cw.getCluster("A") << ',';
	//	cout << cw.getCluster("B") << ',';
	//	cout << cw.getCluster("C") << ',';
	//	cout << cw.getCluster("D") << ',';
	//	cout << cw.getCluster("E") << ',';
	//	cout << cw.getCluster("F") << ',';
	//	cout << cw.getCluster("G") << ',';
	//	cout << cw.getCluster("H") << ',';
	//	cout << cw.getCluster("I") << endl;
	//}





	//AdjacencyList al;
	////al.fromCSV("C:\\Users\\Eric\\Downloads\\enwiki-20110722-csv.tar\\enwiki-20110722-csv\\enwiki-20110722-csv\\SccGraphData.tsv");
	////al.fromCSV("C:\\Users\\Eric\\Downloads\\2016_04_en_clickstream.tsv\\reduced.tsv");
	////al.fromCSV("C:\\Users\\Eric\\Google Drive\\Lethbridge Projects\\Fuzzy Place Field Test\\SimpleTestGraph.txt");
	////al.fromCSV("C:\\Users\\Eric\\Downloads\\web-Google.txt\\reduced.tsv");
	//al.fromCSV("C:\\Users\\Eric\\Downloads\\roadNet-TX.txt\\roadNet-TX.txt");
	////al.fromCSV("reduced.tsv");
	//
	////std::vector<std::vector<int>*> scc;
	////int maxI = al.findSCCs(scc);
	////al = al.subgraph(scc[maxI]);
	////al.toCSV("reduced.tsv");

	//cout << "vertices:  " << al.numVertices() << endl;
	//cout << "edges:     " << al.numEdges() << endl;
	//cout << "abs. Vert: " << al.countAbsorbingVertices() << endl;
	//
	//int numClus = 2000, numRuns = 3000;

	////cout << "clusters: ";
	////cin >> numClus;
	////cout << "runs: ";
	////cin >> numRuns;

	////ChineseWhispers clusterer(&al);
	//HippoCluster clusterer(&al);
	//clusterer.setNumNodes(numClus);
	//clusterer.setTrajectoryLength(al.numVertices());

	//unsigned int start = clock();
	//for (int i = 0; i <= numRuns; i++)
	//{
	//	//cout << "step " << i << endl;
	//		clusterer.step();

	//	if (i % 100 == 0)
	//	{
	//		cout << "step: " << i << ", " << (clock() - start) / CLOCKS_PER_SEC << "s" << endl;
	//		
	//		cout << "  dropped " << clusterer.dropUnusedNodes() << " nodes" << endl;
	//		//cout << "  entropy: " << clusterer.entropy() << endl;

	//		vector<int> assigns;
	//		clusterer.getAllClusterAssignments(assigns);
	//		ClusterStats stats = al.getClusterStats(assigns);

	//		cout << "  relative density: " << stats.getAvgRelativeDensity() << endl;
	//		cout << "  weighted relative density: " << stats.getWeightedAvgRelativeDensity() << endl;
	//		cout << "  mean cluster size: " << stats.getClusterSizeMean() << " ± " << stats.getClusterSizeSD() << " (max " << stats.getClusterSizeMax() << ")" << endl;
	//	}

	//}
	//cout << (clock() - start) / CLOCKS_PER_SEC;


	//clusterer.exportClusterAssignments("clusterAssignments.tsv");
	//cout << "done";

	getchar();
    return 0;
}

