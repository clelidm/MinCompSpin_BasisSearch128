//g++ -std=c++11 main.cpp ReadDataFile.cpp Tool_BitOperations.cpp BestBasis_Init_SetOp.cpp

#include <iostream>
#include <sstream>

#include <map>
#include <set>
#include <list>
#include <vector>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

#include "data.h"

using namespace std;

/********************************************************************/
/**************************    PARAMETERS    ************************/
/********************************************************************/
// number of binary (spin) variables:
unsigned int n = 50;

// INPUT DATA FILES (optional):  
// the input datafile can also be specified directly in the main() function, as an argument of the function "read_datafile()":
std::string input_datafile = "INPUT/Big5-IPC1_VS3_Ne5.dat"; //Big5PT.sorted_Ne5";  //"INPUT/Big5PT.sorted";  //"INPUT/SCOTUS_n9_N895_Data.dat"; //

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
vector<pair<__int128_t, unsigned int>> read_datafile128_vect(string datafilename, unsigned int *N, unsigned int r);

/******************************************************************************/
/**************************     Basis  Tools  *********************************/
/******************************************************************************/
bool Is_Basis(vector<Operator128> Basis, unsigned int n);

void PrintTerm_OpBasis(vector<Operator128> OpVect_Basis, unsigned int n, unsigned int N);
void PrintFile_OpBasis(vector<Operator128> OpVect_Basis, unsigned int n, unsigned int N, string filename);

map<unsigned int, unsigned int> Histo_BasisOpOrder(vector<Operator128> Basis);

/******************************************************************************/
/************************   BASIS SEARCH TOOLS    *****************************/
/******************************************************************************/
// Exhaustive Search:
vector<Operator128> BestBasis_ExhaustiveSearch(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool bool_print = false);

// Fixed Representation up to order `k_max``:
vector<Operator128> BestBasisSearch_FixedRepresentation(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, unsigned int B_it, bool bool_print = false, unsigned int m_max=1000);

// Changing representation up to order `k_max``:
vector<Operator128> BestBasisSearch_Final(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, bool bool_print = false, unsigned int m_max=1000);

/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/

const unsigned int one = 1;

__int128_t UpdateOp_R0(vector<__int128_t> BestBasis_R0, __int128_t Opbin_Ri)
{
  __int128_t Op_R0 = 0;

  unsigned int i = 0;
  while(Opbin_Ri)
  {
    if(Opbin_Ri & one) 
        {   Op_R0 ^= (BestBasis_R0[i]); }
    Opbin_Ri >>= 1;
    i++;
  }

  return Op_R0;
}

std::string int_to_bstring(__int128_t bool_nb, unsigned int r);

void int_to_digits(__int128_t bool_nb, unsigned int r);

int main(int argc, char *argv[])
{
    unsigned int k_max = 2;

	/**********************     READ ARGUMENTS    *********************************/

    // argv[0] contains the name of the datafile, from the current folder (i.e. from the folder containing "data.h");
    // argv[1] contains the number of variables to read;
    if (argc == 3)
    {
      	string input_datafile_buffer = argv[1];
        string n_string_buffer = argv[2];

        input_datafile = "INPUT/" + input_datafile_buffer;
        //n_string_buffer = argv[2];
        n = stoul(n_string_buffer);
    }
    else if (argc == 4)
    {
        string input_datafile_buffer = argv[1];
        string n_string_buffer = argv[2];

        input_datafile = "INPUT/" + input_datafile_buffer;
        //n_string_buffer = argv[2];
        n = stoul(n_string_buffer);
        k_max=stoul(argv[3]);
    }
    else if (argc != 1)
    {
        cout << "The number of arguments must be either 0, 2, or 3" << endl;
        return 0;
    }

    cout << "--->> Create the \"OUTPUT\" Folder: (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    // chrono variables:
	auto start = chrono::system_clock::now(); 
	auto end = chrono::system_clock::now(); 
    chrono::duration<double> elapsed = end - start; 

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  READ THE DATA:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;

	unsigned int N=0;  // will contain the number of datapoints in the dataset

    vector<pair<__int128_t, unsigned int>> Nvect = read_datafile128_vect(input_datafile, &N, n); 

	if (N == 0) { return 0; } // Terminate program if the file can't be found or is empty

/*
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  SEARCH IN THE ORIGINAL BASIS:  *******************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool bool_print = false;
    unsigned int k_max = 2;  // largest order of operators to take into account in each representation
    unsigned int B_it = 0;   // initial basis

    vector<Operator128> BestBasis = BestBasisSearch_FixedRepresentation(Nvect, n, N, k_max, B_it, bool_print);
*/
/*
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************************  EXHAUSTIVE SEARCH:  ***********************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************" << endl;

    cout << "The Exhaustive Search is not recommended for dataset with more than n~20 variables." << endl << endl;

    bool bool_print = false;
    vector<Operator128> BestBasis = BestBasis_ExhaustiveSearch(Nvect, n, N, bool_print);
    //PrintTerm_OpBasis(BestBasis, n, N);  
*/

    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************  SEARCH IN DIFFERENT REPRESENTATIONS:  *****************************";
    cout << endl << "********************  ! Stops when the found basis is Identity !  **************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool bool_print = false;
    unsigned int m_max = 50000;

    vector<Operator128> BestBasis = BestBasisSearch_Final(Nvect, n, N, k_max, bool_print, m_max);
    //PrintTerm_OpBasis(BestBasis, n, N); 

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************  PRINT TERMINAL/FILE FINAL OPERATOR SET:  *************************";
    cout << endl << "*******************************************************************************************" << endl;

    PrintTerm_OpBasis(BestBasis, n, N);  
    //Is_Basis(BestBasis_k2, n);   // this function can check if a set of Operators is in independent set

    Histo_BasisOpOrder(BestBasis);


/*
    cout << "--> Compute and rank all the observables of order 1 (fields).. " << endl; 

    __int128_t one_i = 1;
    for (int i=0; i<n; i++) // All Fields:
    { 
        //Op = Value_Op(un_i, Nvect, Nd);
        //OpSet.insert(Op);
        //if (Op.bias < (*lowest_bias)) { (*lowest_bias) = Op.bias; }
        cout <<  int_to_bstring(one_i, n) << "\t";
        int_to_digits(one_i, n); 
        one_i = one_i << 1;
    }
*/
/*
    vector<__int128_t> BestBasis_R0;

    BestBasis_R0.push_back(1);
    BestBasis_R0.push_back(2);
    BestBasis_R0.push_back(7);

    for(auto& Opbin:BestBasis_R0) { cout << ((uint64_t) Opbin) << endl;}

    cout << endl;

    cout << 3 << ": \t" << (uint64_t) UpdateOp_R0(BestBasis_R0, 3) << endl;
    cout << 5 << ": \t" << (uint64_t) UpdateOp_R0(BestBasis_R0, 5) << endl;
    cout << 4 << ": \t" << (uint64_t) UpdateOp_R0(BestBasis_R0, 4) << endl;
*/
    end = chrono::system_clock::now();  
    elapsed = end - start;
    cout << endl << "Elapsed time (in s): " << elapsed.count() << endl << endl;  

    return 0;
}




