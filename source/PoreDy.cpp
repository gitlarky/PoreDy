//============================================================================
// Name        : T201_UniformFilm.cpp
// Author      : David
// Version     :
// Copyright   : Your copyright notice
// Description : "Invasion-Percolation Drying + 2D Navier-Stokes + Film Effect" in C++, Ansi-style
//============================================================================



#include "GlobalDeclaration.h"

using namespace std;
using namespace Eigen;

std::string caseName;
std::ofstream oflg;

std::vector<point_c> Pt;
std::vector< line_c> Ln;
std::vector<block_c> Bk;

std::vector<  uCell_c> U;
std::vector<  vCell_c> V;
std::vector<  cCell_c> C;
std::vector<   pore_c> P;
std::vector< throat_c> T;
std::vector<cluster_c> Ct;

std::vector<  fCell_c> F;

size_t UVsize(0), CPsize(0), ExtPtSize(0), NetworkPtSize(0);

bitset<4> toBits(const string & str);
/*====================================================================================================
 *
 *                                        main Program
 *
  ====================================================================================================*/
// int main(char* argc[], bitset<4> argv) {
int main(int argc, char* argv[]) {
	cout<<scientific<<setprecision(3)<<setw(10);
	cout<<"Please enter the name of the case: "<<endl;
	if (argc>1) {
		caseName=argv[1];
		cout<<caseName<<endl;
	} else {
		cin>>caseName;
	}
	string logName(caseName+".log");
	oflg.open(logName.c_str());
	oflg<<scientific<<setprecision(3)<<setw(10);

	readCase();

	bitset<4> stage;

	cout<<"What are you planning to do? Please enter a 4 digits binary number to select:"<<endl<<endl;

	cout<<"1st Digit: 1: Generate Mesh        ; 0: Read Mesh      "<<endl;
	cout<<"2nd Digit: 1: Initialize Flow Field; 0: Read Flow Field"<<endl;
	cout<<"3rd Digit: 1: Calculate  Flow Field; 0: Set Flow Field "<<endl;
	cout<<"3rd Digit: 1: Calculate Evaporation; 0: Do Nothing     "<<endl;

	cout<<endl<<"Please choose from the above options (Default=1111): "<<endl;
	if (argc>2) {
		stage=toBits(argv[2]);
		cout<<stage[3]<<stage[2]<<stage[1]<<stage[0]<<endl;
	} else {
		cin>>stage;
	}

	for(size_t i=0; i<stage.size(); ++i) {
		cout<<"stage["<<i<<"]="<<stage[i]<<endl;
	}

	if(stage==0) Testing()             ;

	if(stage[3]) meshing()             ; else readMesh();
	if(stage[2]) initializeFlowfield() ; else readFlowfield();
	if(stage[1]) calculateFlowfield()  ; else setFlowfield();
	if(stage[0]) calculateEvaporation();

	oflg.close();

	return 0;
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                        function readCase
 *
  ----------------------------------------------------------------------------------------------------*/
bool readCase() {
	cout<<"Reading case......"<<endl;

	ifstream caseIFS(caseName.c_str());
	string myLine, myWord;

	while(getline(caseIFS, myLine)) {
		istringstream iss(myLine);
		while(iss>>myWord) {
			if(myWord=="Physics:") {
				iss>>cd_c::phy;
				cd_c::phy.nondimensionalize();
			} else if(myWord=="Point:") {
				point_c pt;
				iss>>pt;
				pt.nondimensionalize(cd_c::phy.RefLength);//All geometry related items are based on these nondimensionalized values
				Pt.push_back(pt);
			} else if(myWord=="Block:") {
				block_c rb;
				iss>>rb;
				rb.nondimensionalize(cd_c::phy);//All BC are based on these nondimensionalized values
				Bk.push_back(rb);
			} else {
				cerr<<"Unrecognizable Input!"<<endl;
			}
		}
	}
	cout<<"Physical and Numerical Setting Read in"<<endl;
	oflg<<cd_c::phy<<endl;
	cout<<"Points Read in and Non-dimensionalized"<<endl;
	cout<<"Blocks Read in and Non-dimensionalized"<<endl;

	return true;
}

bitset<4> toBits(const string & str) {
    return bitset<4>(str);
}
