#include <iostream>
#include <fstream>
#include <string>
#include <new>
#include <boost/algorithm/string.hpp>
#include <sstream>

using namespace std;
using namespace boost;

void input(const string& Input_file_path, int* steps, double* starting, double* step_length, int* gs)
{
    vector <string> fields;
    ifstream input_file(Input_file_path);
    string para_name;
    unsigned int i;

    while(getline(input_file,para_name)) {

        split(fields, para_name, is_any_of(", ="), token_compress_on );

        if (fields.at(0)=="n_lrzcs") {
            for(i = 0; i < 5; i++) {
                steps[i] = stoi(fields.at(i + 1)) + 1;
            }
        }
        else if (fields.at(0)=="lrzcs0"){
            for(i=0; i<5; i++) {
                starting[i] = stof(fields.at(i + 1));
            }
        }
        else if (fields.at(0)=="dlrzcs"){
            for(i=0; i<5; i++){
                step_length[i] = stof(fields.at(i+1));
            }
        }
        else if (fields.at(0)=="gs"){
            for(i=0; i<5; i++){
                gs[i] = stof(fields.at(i+1));
            }
        }
        else{
            input_file.close();
            return;
        }
    }
    input_file.close();
}

void store(const string& PES_file_path, const int* steps, double storation[]){
    vector <string> fields;
    ifstream PES_file(PES_file_path);

    double cm12, Edef, Esh, Eld, AH, rlrzcs[5];
    int Ilrzcs[5];

    while(!PES_file.eof()){
        PES_file>>cm12>>Ilrzcs[0]>>Ilrzcs[1]>>Ilrzcs[2]>>Ilrzcs[3]>>Ilrzcs[4]>>rlrzcs[0]
                >>rlrzcs[1]>>rlrzcs[2]>>rlrzcs[3]>>rlrzcs[4]>>Esh>>Eld>>Edef>>AH;
        storation[Ilrzcs[0]*steps[0]+Ilrzcs[1]*steps[1]+Ilrzcs[2]*steps[2]+Ilrzcs[3]*steps[3]+Ilrzcs[4]] = Edef;
    }
    PES_file.close();
}

//int main()
//{
//    int steps[5];
//    input("fort.112",steps);
//    double* storation = new double[steps[0]*steps[1]*steps[2]*steps[3]*steps[4]];
//    store("/home/dell/U236.txt",steps,storation);
//    return 0;
//}

