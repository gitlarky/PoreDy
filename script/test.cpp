#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

// int main(char* caseName[], bitset<4> & stage) {
//     cout<<caseName<<endl;
//     cout<<stage[3]<<stage[2]<<stage[1]<<stage[0]<<endl;

//     return 0;
// }
int main(string & caseName) {
    cout<<caseName<<endl;

    return 0;
}
// int main() {
//     string filename("");
//     cout<<"Please enter the file name:"<<endl;
//     cin>>filename;
//     cout<<"File name is: "<<filename<<endl;

//     std::ofstream of;
//     of.open(filename.c_str());
//     for(int j=0; j<=9; ++j) {
//         for(int i=0; i<=9; ++i) {
//             of<<j*10+i<<"\t";
//         }
//         of<<"\n";
//     }
//     of.close();

//     // std::ifstream inf;
//     // inf.open(filename.c_str());
//     std::ifstream inf(filename.c_str());
//     std::vector<std::vector<int> > PI;
//     for(int j=0; j<=9; ++j) {
//     	std::vector<int> line;
//         for(int i=0; i<=9; i+=2) {
//         	int v(0);
//             int u(0);
//             inf>>v>>u;
//             line.push_back(v);
//             line.push_back(u);
//         }
//         PI.push_back(line);
//     }
//     for(int j=0; j<=9; ++j) {
//     	for(int i=0; i<=9; ++i) {
//     		std::cout<<PI[j][i]<<"\t";
//     	}
//     	std::cout<<"\n";
//     }

//     return 0;
// }