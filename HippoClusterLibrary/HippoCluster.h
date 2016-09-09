#pragma once
#pragma once
#ifdef HIPPOCLUSTERLIBRARY_EXPORTS
#define HIPPOCLUSTERLIBRARY_API __declspec(dllexport) 
#else
#define HIPPOCLUSTERLIBRARY_API __declspec(dllimport) 
#endif

#include "stdafx.h"
#include <unordered_map>
#include <vector>
#include "AdjacencyList.h"

//#include <mutex>

namespace HippoClusterLibrary
{
	class HippoCluster
	{
	private:
		int numVertices;
		int numNodes = 10;
		std::vector<std::vector<float>> nodeWeights;
		double alpha = 0.01;
		double inhibition = 0.01;
		int trajectoryLength = 1;
		AdjacencyList* adjList;
		//std::mutex mtx;

		void initialize();

		void cosineDist(std::vector<int>& trajectory, double *output);
		void manhatDist(std::vector<double>& trajVector, double* output);
		int getCluster(int vertexNumber);

	public:
		HIPPOCLUSTERLIBRARY_API HippoCluster(AdjacencyList* AdjList);

		int GetCluster(char* vertexName);
		HIPPOCLUSTERLIBRARY_API int getCluster(std::string vertexName);
		HIPPOCLUSTERLIBRARY_API void getAllClusterAssignments(std::vector<int>& assignments);
		HIPPOCLUSTERLIBRARY_API void exportClusterAssignments(std::string fileName);
		HIPPOCLUSTERLIBRARY_API void step();
		void updateWeights(int node, int minIndex, std::vector<int>& traj, std::vector<double>& trajVector);
		HIPPOCLUSTERLIBRARY_API int dropUnusedNodes();

		HIPPOCLUSTERLIBRARY_API double entropy();
		HIPPOCLUSTERLIBRARY_API void setTrajectoryLength(int length);
		HIPPOCLUSTERLIBRARY_API void setNumNodes(int number);

		HIPPOCLUSTERLIBRARY_API int getRepresentativeVertex(int clusterNumber);

		HIPPOCLUSTERLIBRARY_API ~HippoCluster();
	};
}

