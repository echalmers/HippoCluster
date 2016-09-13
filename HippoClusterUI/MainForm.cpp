#include "MainForm.h"
#include "../HippoClusterLibrary/AdjacencyList.h"
#include "../HippoClusterLibrary/HippoCluster.h"
#include <msclr\marshal_cppstd.h>
#include <ctime>
#include "../HippoClusterLibrary/ChineseWhispers.h"
#include "../HippoClusterLibrary/MCL.h"

using namespace HippoClusterLibrary;
using namespace System;
using namespace System::Windows::Forms;

[STAThread]
void Main()
{
	//Application::EnableVisualStyles();
	//Application::SetCompatibleTextRenderingDefault(false);

	HippoClusterUI::MainForm form;
	Application::Run(%form);
}

void HippoClusterUI::MainForm::loadFromBmp()
{
	adjList = new AdjacencyList();

	OpenFileDialog of;
	of.Filter = "Bitmap Files (*.bmp)|*.bmp";
	of.ShowDialog();

	Bitmap^ bmpOrig = gcnew Bitmap(of.FileName);

	Bitmap^ bmp = gcnew Bitmap(bmpOrig->Width*(int)scaleUpDown->Value, bmpOrig->Height*(int)scaleUpDown->Value);
	Graphics^ g = Graphics::FromImage(bmp);
	g->InterpolationMode = System::Drawing::Drawing2D::InterpolationMode::NearestNeighbor;
	g->PixelOffsetMode = System::Drawing::Drawing2D::PixelOffsetMode::Half;
	g->DrawImage(bmpOrig, 0, 0, bmp->Width, bmp->Height);


	for (int i = 1; i < bmp->Width - 1; i++)
	{
		for (int j = 1; j < bmp->Height - 1; j++)
		{
			Color c = bmp->GetPixel(i, j);
			if (bmp->GetPixel(i, j).R == 0 && bmp->GetPixel(i, j).G == 0 && bmp->GetPixel(i, j).B == 255)
				continue;

			std::string V1 = std::to_string(i) + "," + std::to_string(j);
			std::string V2;


			if (bmp->GetPixel(i + 1, j).R != 0 || bmp->GetPixel(i + 1, j).G != 0 && bmp->GetPixel(i + 1, j).B != 255)
			{
				V2 = std::to_string(i + 1) + "," + std::to_string(j);
				adjList->addEdge(V1, V2, 1);
			}
			if (bmp->GetPixel(i - 1, j).R != 0 || bmp->GetPixel(i - 1, j).G != 0 && bmp->GetPixel(i - 1, j).B != 255)
			{
				V2 = std::to_string(i - 1) + "," + std::to_string(j);
				adjList->addEdge(V1, V2, 1);
			}
			if (bmp->GetPixel(i, j + 1).R != 0 || bmp->GetPixel(i, j + 1).G != 0 && bmp->GetPixel(i, j + 1).B != 255)
			{
				V2 = std::to_string(i) + "," + std::to_string(j + 1);
				adjList->addEdge(V1, V2, 1);
			}
			if (bmp->GetPixel(i, j - 1).R != 0 || bmp->GetPixel(i, j - 1).G != 0 && bmp->GetPixel(i, j - 1).B != 255)
			{
				V2 = std::to_string(i) + "," + std::to_string(j - 1);
				adjList->addEdge(V1, V2, 1);
			}
		}
	}


	mainTextBox->AppendText(of.SafeFileName + Environment::NewLine);
	mainTextBox->AppendText("Vertex Count: " + adjList->numVertices().ToString() + Environment::NewLine + "Edge Count: " + adjList->numEdges() + Environment::NewLine);
}

void HippoClusterUI::MainForm::loadFromCsv()
{
	OpenFileDialog of;
	of.Filter = "csv (*.csv, *.txt, *.tsv)|*.csv;*.txt;*.tsv";
	of.ShowDialog();

	adjList = new AdjacencyList;
	adjList->fromTSV(msclr::interop::marshal_as<std::string>(of.FileName));

	mainTextBox->AppendText(of.SafeFileName + Environment::NewLine);
	mainTextBox->AppendText("Vertex Count: " + adjList->numVertices().ToString() + Environment::NewLine + "Edge Count: " + adjList->numEdges() + Environment::NewLine);
}

void HippoClusterUI::MainForm::run()
{
	chart1->Series["rho"]->Points->Clear();
	chart1->Series["H"]->Points->Clear();
	chart1->Series["weightedRho"]->Points->Clear();

	hippoCluster = new HippoCluster(adjList);
	hippoCluster->setNumNodes(100);
	hippoCluster->setTrajectoryLength(adjList->numVertices());

	int colors[1000][3];
	for (int i = 0; i < 1000; i++)
	{
		colors[i][0] = rand() % 255;
		colors[i][1] = rand() % 255;
		colors[i][2] = rand() % 255;
	}

	std::clock_t start;
	double duration;
	start = std::clock();

	for (int stepNumber = 0; stepNumber <= 1000; stepNumber++)
	{

		hippoCluster->step();

		if (stepNumber % 100 == 0)
		{
			mainTextBox->AppendText(hippoCluster->dropUnusedNodes() + Environment::NewLine);
		}

		if (stepNumber % 100 == 0)
		{
			double xStep = 48 * (double)scaleUpDown->Value / (double)pictureBox1->Width;
			double yStep = 16 * (double)scaleUpDown->Value / (double)pictureBox1->Height;
			Bitmap^ img = gcnew Bitmap(pictureBox1->Width, pictureBox1->Height);
			for (int i = 0; i < img->Width; i++)
			{
				for (int j = 0; j < img->Height; j++)
				{
					std::string vertexName = std::to_string((int)(i*xStep)) + "," + std::to_string((int)(j*yStep));

					if (adjList->vertexExists(vertexName))
					{
						int cluster = hippoCluster->getCluster(vertexName);
						img->SetPixel(i, j, Color::FromArgb(colors[cluster][0], colors[cluster][1], colors[cluster][2]));
					}
					else
						img->SetPixel(i, j, Color::Black);
				}
			}
			pictureBox1->Image = img;
			pictureBox1->Refresh();

			double entropy = hippoCluster->entropy();

			std::vector<int> clusters;
			hippoCluster->getAllClusterAssignments(clusters);
			ClusterStats stats = adjList->getClusterStats(clusters);

			chart1->Series["rho"]->Points->AddXY(stepNumber, stats.getAvgRelativeDensity());
			chart1->Series["weightedRho"]->Points->AddXY(stepNumber, stats.getWeightedAvgRelativeDensity());
			chart1->Series["H"]->Points->AddXY(stepNumber, entropy);
			chart1->Refresh();

			//Bitmap^ imgResize = gcnew Bitmap(pictureBox1->Width, pictureBox1->Height);
			//Graphics^ g = Graphics::FromImage(imgResize);
			//g->InterpolationMode = System::Drawing::Drawing2D::InterpolationMode::NearestNeighbor;
			//g->DrawImage(img, 0, 0, pictureBox1->Width, pictureBox1->Height);
			
		}
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	mainTextBox->AppendText(duration.ToString() + "s" + Environment::NewLine);
}

void HippoClusterUI::MainForm::runCW()
{
	chart1->Series["rho"]->Points->Clear();
	chart1->Series["H"]->Points->Clear();
	chart1->Series["weightedRho"]->Points->Clear();

	int colors[1000][3];
	for (int i = 0; i < 1000; i++)
	{
		colors[i][0] = rand() % 255;
		colors[i][1] = rand() % 255;
		colors[i][2] = rand() % 255;
	}

	ChineseWhispers cw(adjList);

	std::clock_t start;
	double duration;
	start = std::clock();

	for (int stepNumber = 0; stepNumber <= 50; stepNumber++)
	{
		cw.step();


		if (stepNumber % 1 == 0)
		{
			double xStep = 48 * (double)scaleUpDown->Value / (double)pictureBox1->Width;
			double yStep = 16 * (double)scaleUpDown->Value / (double)pictureBox1->Height;
			Bitmap^ img = gcnew Bitmap(pictureBox1->Width, pictureBox1->Height);
			for (int i = 0; i < img->Width; i++)
			{
				for (int j = 0; j < img->Height; j++)
				{
					std::string vertexName = std::to_string((int)(i*xStep)) + "," + std::to_string((int)(j*yStep));

					if (adjList->vertexExists(vertexName))
					{
						int cluster = cw.getCluster(vertexName)%1000;
						img->SetPixel(i, j, Color::FromArgb(colors[cluster][0], colors[cluster][1], colors[cluster][2]));
					}
					else
						img->SetPixel(i, j, Color::Black);
				}
			}
			pictureBox1->Image = img;
			pictureBox1->Refresh();

			std::vector<int> clusters;
			cw.getAllClusterAssignments(clusters);
			ClusterStats stats = adjList->getClusterStats(clusters);

			chart1->Series["rho"]->Points->AddXY(stepNumber, stats.getAvgRelativeDensity());
			chart1->Series["weightedRho"]->Points->AddXY(stepNumber, stats.getWeightedAvgRelativeDensity());
			chart1->Refresh();
		}
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	mainTextBox->AppendText(duration.ToString() + "s" + Environment::NewLine);
}

void HippoClusterUI::MainForm::runMCL()
{
	chart1->Series["rho"]->Points->Clear();
	chart1->Series["weightedRho"]->Points->Clear();
	chart1->Series["H"]->Points->Clear();

	MCL mcl(adjList);

	int colors[1000][3];
	for (int i = 0; i < 1000; i++)
	{
		colors[i][0] = rand() % 255;
		colors[i][1] = rand() % 255;
		colors[i][2] = rand() % 255;
	}

	std::clock_t start;
	double duration;
	start = std::clock();

	for (int stepNumber = 0; stepNumber <= 20; stepNumber++)
	{
		mcl.step();


		if (stepNumber % 1 == 0)
		{
			double xStep = 48 * (double)scaleUpDown->Value / (double)pictureBox1->Width;
			double yStep = 16 * (double)scaleUpDown->Value / (double)pictureBox1->Height;
			Bitmap^ img = gcnew Bitmap(pictureBox1->Width, pictureBox1->Height);
			for (int i = 0; i < img->Width; i++)
			{
				for (int j = 0; j < img->Height; j++)
				{
					std::string vertexName = std::to_string((int)(i*xStep)) + "," + std::to_string((int)(j*yStep));

					if (adjList->vertexExists(vertexName))
					{
						int cluster = mcl.getCluster(vertexName)%1000;
						img->SetPixel(i, j, Color::FromArgb(colors[cluster][0], colors[cluster][1], colors[cluster][2]));
					}
					else
						img->SetPixel(i, j, Color::Black);
				}
			}
			pictureBox1->Image = img;
			pictureBox1->Refresh();

			std::vector<int> clusters;
			mcl.getAllClusterAssignments(clusters);
			ClusterStats stats = adjList->getClusterStats(clusters);

			chart1->Series["rho"]->Points->AddXY(stepNumber, stats.getAvgRelativeDensity());
			chart1->Series["weightedRho"]->Points->AddXY(stepNumber, stats.getWeightedAvgRelativeDensity());
			chart1->Refresh();

			mainTextBox->AppendText(mcl.totalElements() + Environment::NewLine);

			//Bitmap^ imgResize = gcnew Bitmap(pictureBox1->Width, pictureBox1->Height);
			//Graphics^ g = Graphics::FromImage(imgResize);
			//g->InterpolationMode = System::Drawing::Drawing2D::InterpolationMode::NearestNeighbor;
			//g->DrawImage(img, 0, 0, pictureBox1->Width, pictureBox1->Height);

		}
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	mainTextBox->AppendText(duration.ToString() + "s" + Environment::NewLine);
}