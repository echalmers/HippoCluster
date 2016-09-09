#include "stdafx.h"
#include "MCL.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <ppl.h>
#include <iostream>

using namespace std;

namespace HippoClusterLibrary
{
	MCL::MCL(AdjacencyList* AdjList)
	{
		adjList = AdjList;

		// get the sparse adjacency matrix
		adjList->toSparseMatrix(Tg);

		// add self loops
		vector<double> rowMax(adjList->numVertices(), 0);
		for (auto row = Tg.begin(); row != Tg.end(); row++)
		{
			int r = row->first;
			for (auto col = Tg[r].begin(); col != Tg[r].end(); col++)
			{
				int c = col->first;
				rowMax[r] = max(rowMax[r], Tg[r][c]);
			}
		}
		for (int i = 0; i < adjList->numVertices(); i++)
		{
			Tg[i][i] = rowMax[i];
		}

		// rescale
		rescaleColumns(Tg);


		T = Tg;
	}


	void MCL::sparseMatrixMultiply(unordered_map<int, unordered_map<int, double>>& A, unordered_map<int, unordered_map<int, double>>& B, unordered_map<int, unordered_map<int, double>>& result, int n)
	{
		concurrency::critical_section cs;

		//concurrency::parallel_for_each(begin(A), end(A), [&](pair<int,unordered_map<int,double>> row)
		for (auto row = A.begin(); row != A.end(); row++)
		{
			int i = row->first;
			for (auto col = A[i].begin(); col != A[i].end(); col++)
			{
				int j = col->first;
				
				for (auto col2 = B[j].begin(); col2 != B[j].end(); col2++)
				//for (int k = 0; k < n; k++)
				{
					int k = col2->first;
					cs.lock();
					result[i][k] += A[i][j] * B[j][k];
					cs.unlock();
				}

			}
		}//);

	}

	void MCL::rescaleColumns(unordered_map<int, unordered_map<int, double>>& mat)
	{
		vector<double> columnSum(adjList->numVertices(), 0);
		for (auto row = mat.begin(); row != mat.end(); row++)
		{
			int r = row->first;
			for (auto col = mat[r].begin(); col != mat[r].end(); col++)
			{
				columnSum[col->first] += col->second;
			}
		}

		for (auto row = mat.begin(); row != mat.end(); row++)
		{
			int r = row->first;
			for (auto col = mat[r].begin(); col != mat[r].end(); col++)
			{
				mat[r][col->first] /= columnSum[col->first];
			}
		}
	}

	void MCL::step()
	{
		unordered_map<int, unordered_map<int, double>> newT;

		// expand
		for (int i = 1; i < expansionFactor; i++)
		{
			sparseMatrixMultiply(T, Tg, newT, adjList->numVertices());
			T = newT;
		}

		// inflate
		inflate();

		// prune
		prune_threshold(pruningThreshold);

		// rescale
		rescaleColumns(T);
	}

	void MCL::prune_resouces(int resources)
	{
		throw "not implemented properly";
		vector<vector<double>> columnVals(adjList->numVertices());

		// get lists of values in each column
		for (auto row = T.begin(); row != T.end(); row++)
		{
			int r = row->first;
			for (auto col = T[r].begin(); col != T[r].end(); col++)
			{
				int c = col->first;
				columnVals[c].push_back(T[r][c]);
			}
		}

		// sort the lists
		for (int i = 0; i < columnVals.size(); i++)
		{
			sort(columnVals[i].begin(), columnVals[i].end());
		}

		// keep only the 'resources' biggest entries in each column
		for (auto row = T.begin(); row != T.end(); row++)
		{
			int r = row->first;
			for (int c = T[r].size() - 1; c >= 0; c--)
			{
				int targetIndex = columnVals[c].size() - resources;
				targetIndex = max(targetIndex, 0);
				if (T[r][c] < columnVals[c][targetIndex])
					T[r].erase(c);
			}
		}
	}

	void MCL::prune_threshold(double threshold)
	{
		// flag entries for removal
		unordered_map<int, vector<int>> removeIndices;
		for (auto row = T.begin(); row != T.end(); row++)
		{
			int r = row->first;
			for (auto col = T[r].begin(); col != T[r].end(); col++)
			{
				int c = col->first;
				if (T[r][c] < threshold)
				{
					removeIndices[r].push_back(c);
				}
			}
		}

		// perform removal
		for (auto row = removeIndices.begin(); row != removeIndices.end(); row++)
		{
			int r = row->first;
			for (int c = 0; c < removeIndices[r].size(); c++)
			{
				T[r].erase(removeIndices[r][c]);
			}
		}
	}


	int MCL::totalElements()
	{
		int tot = 0;
		for (auto row = T.begin(); row != T.end(); row++)
		{
			tot += row->second.size();
		}
		return tot;
	}

	HIPPOCLUSTERLIBRARY_API int MCL::getCluster(std::string vertexName)
	{
		return getCluster(adjList->getVertexNumber(vertexName));
	}

	int MCL::getCluster(int vertexNumber)
	{

		double maxVal = 0; int maxIndex = -1;
		for (int i = 0; i < adjList->numVertices(); i++)
		{
			if (T[i][vertexNumber] > maxVal)
			{
				maxVal = T[i][vertexNumber];
				maxIndex = i;
			}
		}

		return maxIndex;
	}

	HIPPOCLUSTERLIBRARY_API void MCL::getAllClusterAssignments(std::vector<int>& assignments)
	{
		assignments.resize(adjList->numVertices(), adjList->numVertices());

		for (int i = 0; i < assignments.size(); i++)
		{
			assignments[i] = getCluster(i);
		}
	}

	void MCL::inflate()
	{
		for (auto row = T.begin(); row != T.end(); row++)
		{
			int r = row->first;
			for (auto col = T[r].begin(); col != T[r].end(); col++)
			{
				int c = col->first;
				T[r][c] = pow(T[r][c], inflationFactor);
			}
		}
	}

	MCL::~MCL()
	{
	}
}
