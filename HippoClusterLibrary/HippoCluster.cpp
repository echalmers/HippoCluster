#include "stdafx.h"
#include "HippoCluster.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <ppl.h>
//#include <mutex>



namespace HippoClusterLibrary
{
	HippoCluster::HippoCluster(AdjacencyList* AdjList)
	{
		adjList = AdjList;
		numVertices = adjList->numVertices();
		initialize();
	}

	void HippoCluster::setNumNodes(int number)
	{
		numNodes = number;
		initialize();
	}

	void HippoCluster::initialize()
	{
		// initialize with random weights
		nodeWeights.resize(numNodes);
		concurrency::parallel_for(0, numNodes, [&](int i)
			//for (int i = 0; i < numNodes; i++)
		{
			int startPt = numVertices / numNodes * i;
			nodeWeights[i].resize(numVertices);
			for (int j = 0; j < numVertices; j++)
			{
				nodeWeights[i][(startPt + j) % numVertices] = ((double)rand() / RAND_MAX);
			}
		}
		);
	}

	void HippoCluster::exportClusterAssignments(std::string fileName)
	{
		std::ofstream file;
		file.open(fileName);
		for (int i = 0; i < numVertices; i++)
		{
			std::string name = adjList->getVertexName(i);
			int assignment = getCluster(i);
			file << name << '\t' << assignment << std::endl;
		}
		file.close();
	}


	int HippoCluster::getRepresentativeVertex(int clusterNumber)
	{
		int vertexNumber = 0;
		double maxValue = 0;
		for (int w = 0; w < nodeWeights[clusterNumber].size(); w++)
		{
			if (nodeWeights[clusterNumber][w] > maxValue)
			{
				maxValue = nodeWeights[clusterNumber][w];
				vertexNumber = w;
			}
		}
		return vertexNumber;

	}


	void HippoCluster::step()
	{
		// generate a random trajectory
		int startState = (rand() * (RAND_MAX + 1) + rand()) % numVertices;
		//int startState = getRepresentativeVertex(rand() % numNodes);

		std::vector<int> traj;
		std::vector<double> trajVector;
		adjList->randPath(startState, trajectoryLength, traj, trajVector);
		//std::cout << "  trajectory length: " << traj.size() << std::endl;

		// find the closest match neuron
		double *d = new double[numNodes];
		cosineDist(traj, d);
		//manhatDist(trajVector, d);
		int minIndex = 0; double minDist = d[0];
		for (int i = 1; i < numNodes; i++)
		{
			if (d[i] < minDist)
			{
				minDist = d[i];
				minIndex = i;
			}
		}

		// update weights
		/*for (int node = 0; node < numNodes; node++)
		{
			updateWeights(node, minIndex, traj, trajVector);
		}*/
		Concurrency::parallel_for(0, numNodes,
			[&](int node)
			{
				updateWeights(node, minIndex, traj, trajVector);
			}
		);


		delete d;
	}

	void HippoCluster::updateWeights(int node, int minIndex, std::vector<int>& traj, std::vector<double>& trajVector)
	{
		// for the winning node, update all weights
		if (node == minIndex)
		{
			for (int w = 0; w < numVertices; w++)
			{
				nodeWeights[node][w] = nodeWeights[node][w] + alpha * (trajVector[w] - nodeWeights[node][w]);
			}
		}
		// for other nodes, inhibit connections corresponding to this trajectory
		else
		{
			for (int w = 0; w < traj.size(); w++) // should this operate on the list of states (which could include duplicates) or the vector? *********
			{
				nodeWeights[node][traj[w]] = nodeWeights[node][traj[w]] * (1 - inhibition);
			}
		}
	}


	int HippoCluster::dropUnusedNodes()
	{
		CRITICAL_SECTION criticalSection;
		InitializeCriticalSection(&criticalSection);

		int dropped = 0;

		double *clusterCnts = new double[numNodes];
		for (int n = 0; n < numNodes; n++)
		{
			clusterCnts[n] = 0;
		}

		concurrency::parallel_for(0, numVertices, [&](int v)
		//for (int v = 0; v < numVertices; v++)
		{
			int n = getCluster(v);
			
			EnterCriticalSection(&criticalSection);
			clusterCnts[n]++;
			LeaveCriticalSection(&criticalSection);
		}
		);

		for (int n = numNodes - 1; n >= 0; n--)
		{
			if (clusterCnts[n] == 0)
			{
				nodeWeights.erase(nodeWeights.begin() + n);
				numNodes--;
				dropped++;
			}
		}

		delete clusterCnts;
		DeleteCriticalSection(&criticalSection);

		return dropped;
	}

	void HippoCluster::setTrajectoryLength(int length)
	{
		trajectoryLength = length;
	}

	void HippoCluster::cosineDist(std::vector<int>& trajectory, double *output)
	{
		double trajMag = sqrt(trajectory.size());


		concurrency::parallel_for(0,numNodes,
			[&] (int node)
		//for (int node = 0; node < numNodes; node++)
		{
			double weightMag = 0;
			int index = 0;
			for (int v = 0; v < numVertices; v++)
			{
				weightMag += nodeWeights[node][v] * nodeWeights[node][v];
			}
			weightMag = sqrt(weightMag);

			double dot = 0;
			for (int index = 0; index < trajectory.size(); index++)
			{
				dot += nodeWeights[node][trajectory[index]];
			}

			output[node] = 1 - dot / weightMag / trajMag;
		}
		);
	}

	void HippoCluster::manhatDist(std::vector<double>& trajVector, double* output)
	{
		
		concurrency::parallel_for(0, numNodes,
			[&](int node)
		{
			double d = 0;
			for (int i = 0; i < trajVector.size(); i++)
			{
				d += abs(trajVector[i] - nodeWeights[node][i]);
			}
			output[node] = d;
		}
		);
	}

	int HippoCluster::getCluster(int vertexNumber)
	{
		/* should we be getting the max weight, or dividing by magnitude as in cosine dist? ***/
		int clusterNum = 0; double maxWeight = nodeWeights[0][vertexNumber];
		for (int i = 1; i < numNodes; i++)
		{
			if (nodeWeights[i][vertexNumber] > maxWeight)
			{
				maxWeight = nodeWeights[i][vertexNumber];
				clusterNum = i;
			}
		}
		return clusterNum;
	}

	int HippoCluster::GetCluster(char* vertexName)
	{
		return getCluster(std::string(vertexName));
	}

	int HippoCluster::getCluster(std::string vertexName)
	{
		int vertexNumber = adjList->getVertexNumber(vertexName);
		return getCluster(vertexNumber);
	}

	double HippoCluster::entropy()
	{
		std::vector<double> h(numNodes, 0);

		concurrency::parallel_for(0, numNodes,
			[&](int i)
		{
				for (int j = 0; j < nodeWeights[i].size(); j++)
				{
					double p = nodeWeights[i][j];
					h[i] += -3.5727 * p * p + 3.5727 * p + 0.127;
				}
		}
		);

		double H = 0;
		for (int i = 0; i < nodeWeights.size(); i++)
		{
			H += h[i];
		}
		return H / numNodes / numVertices;
	}

	void HippoCluster::getAllClusterAssignments(std::vector<int>& assignments)
	{
		// compute all cluster assignments
		assignments.resize(numVertices);
		concurrency::parallel_for(0, numVertices, [&](int v)
			//for (int v = 0; v < numVertices; v++)
		{
			assignments[v] = getCluster(v);
		}
		);
	}

	HippoCluster::~HippoCluster()
	{
		delete adjList;
	}
}