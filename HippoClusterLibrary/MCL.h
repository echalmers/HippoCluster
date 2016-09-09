#pragma once
#ifdef HIPPOCLUSTERLIBRARY_EXPORTS
#define HIPPOCLUSTERLIBRARY_API __declspec(dllexport) 
#else
#define HIPPOCLUSTERLIBRARY_API __declspec(dllimport) 
#endif


#include "AdjacencyList.h"
#include <unordered_map>

using namespace std;

namespace HippoClusterLibrary
{
	class MCL
	{
	private:
		double eta = 0.000001;
		AdjacencyList* adjList;
		int expansionFactor = 20;
		double inflationFactor = 1.5;
		double pruningThreshold = 0.0;
		

		void sparseMatrixMultiply(unordered_map<int, unordered_map<int, double>>& A, unordered_map<int, unordered_map<int, double>>& B, unordered_map<int, unordered_map<int, double>>& result, int n);
		void rescaleColumns(unordered_map<int, unordered_map<int, double>>& mat);
		void inflate();
		void prune_threshold(double threshold);
		void prune_resouces(int resources);
		int getCluster(int vertexNumber);

	public:
		HIPPOCLUSTERLIBRARY_API MCL(AdjacencyList* AdjList);
unordered_map<int, unordered_map<int, double>> Tg;
unordered_map<int, unordered_map<int, double>> T;
HIPPOCLUSTERLIBRARY_API int totalElements();
		HIPPOCLUSTERLIBRARY_API void step();
		HIPPOCLUSTERLIBRARY_API int getCluster(std::string vertexName);
		HIPPOCLUSTERLIBRARY_API void getAllClusterAssignments(std::vector<int>& assignments);

		HIPPOCLUSTERLIBRARY_API ~MCL();
	};
}

