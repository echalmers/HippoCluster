#pragma once
#pragma once
#ifdef HIPPOCLUSTERLIBRARY_EXPORTS
#define HIPPOCLUSTERLIBRARY_API __declspec(dllexport) 
#else
#define HIPPOCLUSTERLIBRARY_API __declspec(dllimport) 
#endif

#include <vector>
#include "AdjacencyList.h"

namespace HippoClusterLibrary
{
	class ChineseWhispers
	{
	private:
		std::vector<int> clusterAssignments;
		AdjacencyList* adjList;
		int numVertices;

	public:
		HIPPOCLUSTERLIBRARY_API ChineseWhispers(AdjacencyList* AdjList);

		HIPPOCLUSTERLIBRARY_API void step();
		HIPPOCLUSTERLIBRARY_API int getCluster(std::string vertexName);
		HIPPOCLUSTERLIBRARY_API void getAllClusterAssignments(std::vector<int>& assignments);

		HIPPOCLUSTERLIBRARY_API void exportClusterAssignments(std::string filename);

		HIPPOCLUSTERLIBRARY_API ~ChineseWhispers();
	};
}
