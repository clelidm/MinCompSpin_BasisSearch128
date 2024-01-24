//g++ -std=c++11 -O3 src/*.cpp includes/main.cpp -o BestBasis.out

#include <iostream>
#include <sstream>

#include <map>
#include <set>
#include <list>
#include <vector>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

#include "../src/data.h"

using namespace std;

/********************************************************************/
/**************************    PARAMETERS    ************************/
/********************************************************************/
// number of binary (spin) variables:
unsigned int n = 9;

// INPUT DATA FILES (optional):  must be in the INPUT folder
string input_datafile = "INPUT/Shapes_n9_Dataset_N1e5.dat"; 

unsigned int k_max = 3; // default value

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
vector<Operator128> BestBasisSearch_FixedRepresentation(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max = 2, bool bool_print = false, unsigned int B_it = 0, unsigned int m_max=1000);

// Changing representation up to order `k_max``:
vector<Operator128> BestBasisSearch_Final(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max=2, bool bool_print = false, unsigned int m_max=1000);


/******************************************************************************/
/************************ User Interface with Flags ***************************/
/******************************************************************************/

int Read_argument(int argc, char *argv[], string *input_datafile, unsigned int *n, unsigned int *k_max);

/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/

int main(int argc, char *argv[])
{
/**********************     READ ARGUMENTS    *********************************/
    // argv[0] contains the name of the datafile, from the current folder (i.e. from the folder containing "data.h");
    // argv[1] contains the number of variables to read;
    // argv[2] contains flag
    // argv[3] contains kmax

/**********************     CREATE FLAG    ************************************/
// By default:  flag_search = 1 (for the example)
    // 1 = Exhaustive search
    // 2 = Fixed basis search with given choice of k_max
    // 3 = Varying basis search with given choice of k_max

    int flag_search = Read_argument(argc, argv, &input_datafile, &n, &k_max);


/**********************   CREATE OUTPUT DIRECTORY    ***************************/

    cout << "--->> Create the \"OUTPUT\" Folder: (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    // chrono variables:
	auto start = chrono::system_clock::now(); 
	auto end = chrono::system_clock::now(); 
    chrono::duration<double> elapsed = end - start; 

    //Variables:
    bool bool_print = false;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  READ THE DATA:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;

	unsigned int N=0;  // will contain the number of datapoints in the dataset

    vector<pair<__int128_t, unsigned int>> Nvect = read_datafile128_vect(input_datafile, &N, n); 

	if (N == 0) { return 0; } // Terminate program if the file can't be found or is empty


    vector<Operator128> BestBasis;

    if (flag_search == 1)
    {
        cout << endl << "*******************************************************************************************";
        cout << endl << "*******************************************************************************************";
        cout << endl << "**********************************  EXHAUSTIVE SEARCH:  ***********************************";
        cout << endl << "*******************************************************************************************";
        cout << endl << "*******************************************************************************************" << endl;

        cout << endl << "Information: " << endl;
        cout << "  > The following program search for the Best Basis among all the (2^n-1) possible operators." << endl;

        cout << endl << "Important: " << endl;
        cout << "  > The Exhaustive Search is not recommended for dataset with more than n~20 variables." << endl << endl;

        bool_print = false;

        BestBasis = BestBasis_ExhaustiveSearch(Nvect, n, N, bool_print);
    }

    else if (flag_search == 2)
    {
        cout << endl << "*******************************************************************************************";
        cout << endl << "*******************************************************************************************";
        cout << endl << "******************************  SEARCH UP TO ORDER K  *************************************";
        cout << endl << "****************************  IN FIXED REPRESENTATION:  ***********************************";
        cout << endl << "*******************************************************************************************";
        cout << endl << "*******************************************************************************************" << endl;

        bool_print = false;
        // By default, k_max = 3;  // largest order of operators to take into account in each representation

        cout << "Search for the best basis among all operators up to order kmax = " << k_max << "." << endl << endl;

        BestBasis = BestBasisSearch_FixedRepresentation(Nvect, n, N, k_max, bool_print);
    }

    else if (flag_search == 3)
    {
        cout << endl << "*******************************************************************************************";
        cout << endl << "*******************************************************************************************";
        cout << endl << "******************************  SEARCH UP TO ORDER K  *************************************";
        cout << endl << "************************   IN SUCCCESSIVE REPRESENTATIONS  ********************************";
        cout << endl << "********************  ! Stops when the found basis is Identity !  *************************";
        cout << endl << "*******************************************************************************************";
        cout << endl << "*******************************************************************************************" << endl;

        bool_print = false;
        unsigned int m_max = 50000;
        // By default, k_max = 3;  // largest order of operators to take into account in each representation

        cout << "Search for the best basis among all operators up to order kmax = " << k_max << "." << endl << endl;

        cout << "The data is then transformed in the representation given by the best basis." << endl;
        cout << "The search for the best basis up to order kmax is re-iterated in this new representation." << endl;
        cout << "The process is repeated until the basis doesn't change anymore" << endl;
        cout << "(i.e. the basis found in the current representation is identity)." << endl;

        BestBasis = BestBasisSearch_Final(Nvect, n, N, k_max, bool_print, m_max); 
    }

    if (BestBasis.size() == 0)  // Terminate program if the Basis is empty
    {
        cout << "ERROR: No basis were found. Check the argument provided." << endl;
        return 0; 
    }

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************  PRINT TERMINAL/FILE FINAL OPERATOR SET:  *************************";
    cout << endl << "*******************************************************************************************" << endl;

    PrintTerm_OpBasis(BestBasis, n, N);  
    //Is_Basis(BestBasis_k2, n);   // this function can check if a set of Operators is in independent set

    Histo_BasisOpOrder(BestBasis);


    end = chrono::system_clock::now();  
    elapsed = end - start;
    cout << endl << "Elapsed time (in s): " << elapsed.count() << endl << endl;  

    return 0;
}




