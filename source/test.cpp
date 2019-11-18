#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;



int main() {
    string filename("test.at");

    std::ofstream of;
    of.open(filename.c_str());
    for(int j=0; j<=9; ++j) {
        for(int i=0; i<=9; ++i) {
            of<<j*10+i<<"\t";
        }
        of<<"\n";
    }
    of.close();

    std::ifstream inf;
    inf.open(filename.c_str());
    std::vector<std::vector<int> > PI;
    for(int j=0; j<=9; ++j) {
    	std::vector<int> line;
        for(int i=0; i<=9; ++i) {
        	int v(0);
            inf>>v;
            line.push_back(v);
        }
        PI.push_back(line);
    }
    for(int j=0; j<=9; ++j) {
    	for(int i=0; i<=9; ++i) {
    		std::cout<<PI[j][i]<<"\t";
    	}
    	std::cout<<"\n";
    }



    return 0;
}