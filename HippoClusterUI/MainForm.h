#pragma once
#include "../HippoClusterLibrary/AdjacencyList.h"
#include "../HippoClusterLibrary/HippoCluster.h"

namespace HippoClusterUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for MainForm
	/// </summary>
	public ref class MainForm : public System::Windows::Forms::Form
	{
		
		HippoClusterLibrary::AdjacencyList* adjList;
	private: System::Windows::Forms::NumericUpDown^  scaleUpDown;
	private: System::Windows::Forms::ContextMenuStrip^  contextMenuStrip1;
	private: System::Windows::Forms::ToolStripMenuItem^  saveToolStripMenuItem;
	private: System::Windows::Forms::Button^  runCWbutton;
	private: System::Windows::Forms::Button^  runMCLbutton;
			 HippoClusterLibrary::HippoCluster* hippoCluster;

		void loadFromBmp();
		void loadFromCsv();
		void run();
		void runCW();
		void runMCL();

	public:
		MainForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MainForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Button^  fromBmpButton;
	protected:
	private: System::Windows::Forms::Button^  fromCsvButton;
	private: System::Windows::Forms::Button^  runButton;
	private: System::Windows::Forms::TextBox^  mainTextBox;
	private: System::Windows::Forms::PictureBox^  pictureBox1;
	private: System::Windows::Forms::DataVisualization::Charting::Chart^  chart1;
	private: System::ComponentModel::IContainer^  components;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			System::Windows::Forms::DataVisualization::Charting::ChartArea^  chartArea2 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^  legend2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^  series4 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^  series5 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^  series6 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			this->fromBmpButton = (gcnew System::Windows::Forms::Button());
			this->fromCsvButton = (gcnew System::Windows::Forms::Button());
			this->runButton = (gcnew System::Windows::Forms::Button());
			this->mainTextBox = (gcnew System::Windows::Forms::TextBox());
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			this->contextMenuStrip1 = (gcnew System::Windows::Forms::ContextMenuStrip(this->components));
			this->saveToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->scaleUpDown = (gcnew System::Windows::Forms::NumericUpDown());
			this->runCWbutton = (gcnew System::Windows::Forms::Button());
			this->runMCLbutton = (gcnew System::Windows::Forms::Button());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			this->contextMenuStrip1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->scaleUpDown))->BeginInit();
			this->SuspendLayout();
			// 
			// fromBmpButton
			// 
			this->fromBmpButton->Location = System::Drawing::Point(12, 12);
			this->fromBmpButton->Name = L"fromBmpButton";
			this->fromBmpButton->Size = System::Drawing::Size(75, 23);
			this->fromBmpButton->TabIndex = 0;
			this->fromBmpButton->Text = L"From bmp";
			this->fromBmpButton->UseVisualStyleBackColor = true;
			this->fromBmpButton->Click += gcnew System::EventHandler(this, &MainForm::fromBmpButton_Click);
			// 
			// fromCsvButton
			// 
			this->fromCsvButton->Location = System::Drawing::Point(12, 41);
			this->fromCsvButton->Name = L"fromCsvButton";
			this->fromCsvButton->Size = System::Drawing::Size(75, 23);
			this->fromCsvButton->TabIndex = 1;
			this->fromCsvButton->Text = L"From csv";
			this->fromCsvButton->UseVisualStyleBackColor = true;
			this->fromCsvButton->Click += gcnew System::EventHandler(this, &MainForm::fromCsvButton_Click);
			// 
			// runButton
			// 
			this->runButton->Location = System::Drawing::Point(12, 96);
			this->runButton->Name = L"runButton";
			this->runButton->Size = System::Drawing::Size(75, 23);
			this->runButton->TabIndex = 2;
			this->runButton->Text = L"Run HC";
			this->runButton->UseVisualStyleBackColor = true;
			this->runButton->Click += gcnew System::EventHandler(this, &MainForm::runButton_Click);
			// 
			// mainTextBox
			// 
			this->mainTextBox->Location = System::Drawing::Point(12, 230);
			this->mainTextBox->Multiline = true;
			this->mainTextBox->Name = L"mainTextBox";
			this->mainTextBox->ScrollBars = System::Windows::Forms::ScrollBars::Both;
			this->mainTextBox->Size = System::Drawing::Size(187, 155);
			this->mainTextBox->TabIndex = 3;
			// 
			// pictureBox1
			// 
			this->pictureBox1->Location = System::Drawing::Point(228, 230);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(387, 155);
			this->pictureBox1->TabIndex = 4;
			this->pictureBox1->TabStop = false;
			// 
			// chart1
			// 
			chartArea2->Name = L"ChartArea1";
			this->chart1->ChartAreas->Add(chartArea2);
			this->chart1->ContextMenuStrip = this->contextMenuStrip1;
			legend2->Name = L"Legend1";
			this->chart1->Legends->Add(legend2);
			this->chart1->Location = System::Drawing::Point(228, 12);
			this->chart1->Name = L"chart1";
			series4->ChartArea = L"ChartArea1";
			series4->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series4->Legend = L"Legend1";
			series4->Name = L"rho";
			series5->ChartArea = L"ChartArea1";
			series5->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series5->Color = System::Drawing::Color::Red;
			series5->Legend = L"Legend1";
			series5->Name = L"H";
			series6->ChartArea = L"ChartArea1";
			series6->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series6->Color = System::Drawing::Color::FromArgb(static_cast<System::Int32>(static_cast<System::Byte>(0)), static_cast<System::Int32>(static_cast<System::Byte>(192)),
				static_cast<System::Int32>(static_cast<System::Byte>(0)));
			series6->Legend = L"Legend1";
			series6->Name = L"weightedRho";
			this->chart1->Series->Add(series4);
			this->chart1->Series->Add(series5);
			this->chart1->Series->Add(series6);
			this->chart1->Size = System::Drawing::Size(387, 212);
			this->chart1->TabIndex = 5;
			this->chart1->Text = L"chart1";
			// 
			// contextMenuStrip1
			// 
			this->contextMenuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->saveToolStripMenuItem });
			this->contextMenuStrip1->Name = L"contextMenuStrip1";
			this->contextMenuStrip1->Size = System::Drawing::Size(99, 26);
			// 
			// saveToolStripMenuItem
			// 
			this->saveToolStripMenuItem->Name = L"saveToolStripMenuItem";
			this->saveToolStripMenuItem->Size = System::Drawing::Size(98, 22);
			this->saveToolStripMenuItem->Text = L"Save";
			this->saveToolStripMenuItem->Click += gcnew System::EventHandler(this, &MainForm::saveToolStripMenuItem_Click);
			// 
			// scaleUpDown
			// 
			this->scaleUpDown->Location = System::Drawing::Point(93, 12);
			this->scaleUpDown->Name = L"scaleUpDown";
			this->scaleUpDown->Size = System::Drawing::Size(88, 20);
			this->scaleUpDown->TabIndex = 6;
			this->scaleUpDown->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 0 });
			// 
			// runCWbutton
			// 
			this->runCWbutton->Location = System::Drawing::Point(12, 125);
			this->runCWbutton->Name = L"runCWbutton";
			this->runCWbutton->Size = System::Drawing::Size(75, 23);
			this->runCWbutton->TabIndex = 7;
			this->runCWbutton->Text = L"Run CW";
			this->runCWbutton->UseVisualStyleBackColor = true;
			this->runCWbutton->Click += gcnew System::EventHandler(this, &MainForm::runCWbutton_Click);
			// 
			// runMCLbutton
			// 
			this->runMCLbutton->Location = System::Drawing::Point(12, 154);
			this->runMCLbutton->Name = L"runMCLbutton";
			this->runMCLbutton->Size = System::Drawing::Size(75, 23);
			this->runMCLbutton->TabIndex = 8;
			this->runMCLbutton->Text = L"Run MCL";
			this->runMCLbutton->UseVisualStyleBackColor = true;
			this->runMCLbutton->Click += gcnew System::EventHandler(this, &MainForm::runMCLbutton_Click);
			// 
			// MainForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(643, 419);
			this->Controls->Add(this->runMCLbutton);
			this->Controls->Add(this->runCWbutton);
			this->Controls->Add(this->scaleUpDown);
			this->Controls->Add(this->chart1);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->mainTextBox);
			this->Controls->Add(this->runButton);
			this->Controls->Add(this->fromCsvButton);
			this->Controls->Add(this->fromBmpButton);
			this->Name = L"MainForm";
			this->Text = L"MainForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			this->contextMenuStrip1->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->scaleUpDown))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void fromBmpButton_Click(System::Object^  sender, System::EventArgs^  e) {
		loadFromBmp();
	}
private: System::Void fromCsvButton_Click(System::Object^  sender, System::EventArgs^  e) {
	loadFromCsv();
}
private: System::Void runButton_Click(System::Object^  sender, System::EventArgs^  e) {
	run();
}
private: System::Void saveToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) {
	SaveFileDialog sf;
	sf.ShowDialog();

	pictureBox1->Image->Save(sf.FileName, System::Drawing::Imaging::ImageFormat::Jpeg);
}
private: System::Void runCWbutton_Click(System::Object^  sender, System::EventArgs^  e) {
	runCW();
}
private: System::Void runMCLbutton_Click(System::Object^  sender, System::EventArgs^  e) {
	runMCL();
}
};
}
