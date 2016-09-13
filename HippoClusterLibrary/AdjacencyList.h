#pragma once
#ifdef HIPPOCLUSTERLIBRARY_EXPORTS
#define HIPPOCLUSTERLIBRARY_API __declspec(dllexport) 
#else
#define HIPPOCLUSTERLIBRARY_API __declspec(dllimport) 
#endif

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>
#include <stack>

namespace HippoClusterLibrary
{
	class ClusterStats
	{
		friend class AdjacencyList;
	private:
		double avgRelativeDensity = NAN;
		double weightedAvgRelativeDensity = NAN;
		double clusterSizeMean = NAN;
		double clusterSizeSD = NAN;
		double clusterSizeMax = NAN;
		std::vector<double> relativeDensities;

	public:
		HIPPOCLUSTERLIBRARY_API ClusterStats();

		HIPPOCLUSTERLIBRARY_API double getAvgRelativeDensity();
		HIPPOCLUSTERLIBRARY_API double getWeightedAvgRelativeDensity();
		HIPPOCLUSTERLIBRARY_API double getClusterSizeMean();
		HIPPOCLUSTERLIBRARY_API double getClusterSizeSD();
		HIPPOCLUSTERLIBRARY_API double getClusterSizeMax();
		HIPPOCLUSTERLIBRARY_API std::vector<double>& getRelativeDensities();

		HIPPOCLUSTERLIBRARY_API ~ClusterStats();
	};

	class AdjacencyList
	{
	private:
		double epsilon = 0.000001;

		std::vector<std::vector<std::tuple<int, double, double>>> adjList;
		std::vector<std::string> vertexNames;
		std::unordered_map<std::string, int> vertexNumbers;

		void calculateNeighborProbabilities(int startVertex);
		void strongConnect(int v, long& index, std::stack<int>& S, std::vector<int>& indices, std::vector<int>& lowLinks, std::vector<bool>& onStack, std::vector<std::vector<int>*>& SCCs);

		

	public:
		HIPPOCLUSTERLIBRARY_API AdjacencyList();

		HIPPOCLUSTERLIBRARY_API void fromTSV(std::string filename);
		HIPPOCLUSTERLIBRARY_API void toTSV(std::string filename);
		HIPPOCLUSTERLIBRARY_API void addEdge(std::string vert1, std::string vert2, double count, bool updateProbabilities = true);
		HIPPOCLUSTERLIBRARY_API void calculateNeighborProbabilities();

		HIPPOCLUSTERLIBRARY_API int randNeighbor(int vertex);
		HIPPOCLUSTERLIBRARY_API std::vector<std::tuple<int, double, double>> allNeighbors(int vertex);
		HIPPOCLUSTERLIBRARY_API std::vector<std::tuple<int, double, double>> isNeighborOf(int vertex);
		HIPPOCLUSTERLIBRARY_API void randPath(int startVertex, int length, std::vector<int>& path, std::vector<double>& pathVector);

		HIPPOCLUSTERLIBRARY_API std::string getVertexName(int vertexNumber);

		HIPPOCLUSTERLIBRARY_API bool vertexExists(std::string vertexName);
		HIPPOCLUSTERLIBRARY_API int getVertexNumber(std::string vertexName);

		HIPPOCLUSTERLIBRARY_API int numVertices();
		HIPPOCLUSTERLIBRARY_API int numEdges();
		HIPPOCLUSTERLIBRARY_API double sumWeights();
		HIPPOCLUSTERLIBRARY_API int countAbsorbingVertices();
		HIPPOCLUSTERLIBRARY_API void countSubgraphs(std::unordered_map<int, int>& subgraphSizes, std::vector<int>& subgraphAssignments);
		HIPPOCLUSTERLIBRARY_API void removeSubgraphs(std::unordered_map<int, int>& subgraphSizes, std::vector<int>& subgraphAssignments, int minSize);
		HIPPOCLUSTERLIBRARY_API int findSCCs(std::vector<std::vector<int>*>& SCCs);
		HIPPOCLUSTERLIBRARY_API AdjacencyList subgraph(std::vector<int>* verticesToKeep);

		HIPPOCLUSTERLIBRARY_API void toSparseMatrix(std::unordered_map<int, std::unordered_map<int, double>>& mat);
		
		HIPPOCLUSTERLIBRARY_API ClusterStats getClusterStats(std::vector<int>& clusterAssignments);

		HIPPOCLUSTERLIBRARY_API ~AdjacencyList();
	};

	

}

