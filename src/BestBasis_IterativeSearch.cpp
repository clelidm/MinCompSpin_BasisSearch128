#include <iostream>
#include <sstream>
#include <fstream>

#include <set>
#include <vector>
#include <list>
#include <map>

#include "data.h"

using namespace std;

/******************************************************************************/
/****************     Initial Choice of Operators for Basis    ****************/
/******************    All operators of order k or smaller    *****************/
/******************************************************************************/
set<Operator128> All_Op_k1(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, double *lowest_bias, bool print = false);

void Add_AllOp_kbits_MostBiased(set<Operator128>& OpSet, vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k = 2, double Bias_LowerBound=0, bool print = false);

//void PrintTerm_OpSet(set<Operator128> OpSet);
void PrintTerm_OpSet(set<Operator128> OpSet, unsigned int n);
void PrintFile_OpSet(set<Operator128> OpSet, unsigned int n, string filename);

// Remove the Operator with too small Bias:
void CutSmallBias(set<Operator128>& OpSet, Struct_LowerBound LB);

unsigned int bitset_count(__int128_t bool_nb);
void int_to_digits_file(__int128_t bool_nb, unsigned int r, std::fstream &file);
std::string int_to_bstring(__int128_t bool_nb, unsigned int r);

/******************************************************************************/
/**************************     Select Best Basis    **************************/
/******************************************************************************/
vector<Operator128> BestBasis_inOpSet(set<Operator128> OpSet, unsigned int n, Struct_LowerBound* LowerBound, unsigned int m=1000);

/******************************************************************************/
/**************************     Basis  Tools  *********************************/
/******************************************************************************/
bool Is_Basis(vector<Operator128> Basis, unsigned int n);

void PrintTerm_OpBasis(vector<Operator128> OpVect_Basis, unsigned int n, unsigned int N);
void PrintFile_OpBasis(vector<Operator128> OpVect_Basis, unsigned int n, unsigned int N, string filename);

map<unsigned int, unsigned int> Histo_BasisOpOrder(vector<Operator128> Basis);

bool Check_Basis_Identity(vector<Operator128> Basis)
{
  bool check = true;

  for (auto& Op : Basis) {
    if(bitset_count(Op.bin) != 1)   { check = false; }
  }
  return check;
}

void SaveFile_Basis(vector<Operator128> Basis, unsigned int n, fstream &file)
{
  file << "########### New Basis: \t Total number of operators = " << Basis.size() << endl;  

  for (auto& Op : Basis)
  {
    file <<  int_to_bstring(Op.bin, n) << "\t Bias = " << Op.bias << "\t"; 
    int_to_digits_file(Op.bin, n, file);
  }

  file << endl;
}

/******************************************************************************/
/***************     Search in a Given Representation  Tools  *****************/
/******************************************************************************/

vector<Operator128> BestBasisSearch_FixedRepresentation(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max = 2, bool bool_print = false, unsigned int R_it = 0, unsigned int m_max=1000)
{
  k_max = (k_max<2)?2:k_max;  // k_max must be at least 2;

  cout << endl << "*****************  FIND THE SMALLEST BIAS OF THE CURRENT BASIS (k = 1):  ******************";
  cout << endl << "*******************************************************************************************" << endl;

  // Compute the bias of the basis elements, and find the least biased one:
  // This value will serve as a lower bound for operators that will be kept later on.

  double Bias_LowerBound = 0.;  // Current lower bound (current lowest bias) is 0. --> we accept all possible bias

  set<Operator128> OpSet = All_Op_k1(Nvect, n, N, &Bias_LowerBound, bool_print);

  //PrintTerm_OpSet(OpSet_B0, n);
  PrintFile_OpSet(OpSet, n, "R" + to_string(R_it) +"_k1");

  Struct_LowerBound LB; // Lower Bound Info
  LB.Bias = Bias_LowerBound;
  vector<Operator128> BestBasis;

/*  for (auto& Op: OpSet){
    BestBasis.push_back(Op);
  }
  cout << "Check Basis Identity: k=1: " << Check_Basis_Identity(BestBasis) << endl << endl; 
*/

  string filename_k = "";

  for (unsigned int k = 2; k <= k_max; k++)
  {
      filename_k = "R" + to_string(R_it) +"_k"+to_string(k);

      cout << endl << "****************************  ADD ALL OPERATORS for k = " << k << "  ********************************";
      cout << endl << "*******************************************************************************************" << endl;

      Add_AllOp_kbits_MostBiased(OpSet, Nvect, n, N, k, LB.Bias, bool_print);

      //PrintTerm_OpSet(OpSet, n);
      PrintFile_OpSet(OpSet, n, filename_k);

      cout << endl << "******************************  SEARCH FOR BEST BASIS:  ***********************************" << endl;
//    cout << endl << "*******************************************************************************************" << endl;

      BestBasis.clear();
      BestBasis = BestBasis_inOpSet(OpSet, n, &LB, m_max); // LB will be over-written with the updated values

      PrintFile_OpBasis(BestBasis, n, N, filename_k + "_BestBasis");

      CutSmallBias(OpSet, LB);

      PrintFile_OpSet(OpSet, n, filename_k + "_CutSmallBias");
  }

  cout << endl << "*************************  SEARCH IN GIVEN REPRESENTATION: DONE  **************************"; 
  cout << endl << "*******************************************************************************************" << endl;

  return BestBasis;
}

/******************************************************************************/
/****************     Update REPRESENTATION of the BASIS   ********************/
/******************************************************************************/
const unsigned int one = 1;

__int128_t UpdateOp_inR0(vector<Operator128> BestBasis_R0, __int128_t Opbin_Ri)
{
  __int128_t Op_R0 = 0;

  unsigned int i = 0;
  while(Opbin_Ri)
  {
    if(Opbin_Ri & one) 
        {   Op_R0 ^= (BestBasis_R0[i].bin); }
    Opbin_Ri >>= 1;
    i++;
  }

  return Op_R0;
}

vector<Operator128> UpdateBasis_inR0(vector<Operator128> BestBasis_R0_old, vector<Operator128> BestBasis_Ri)
{
  vector<Operator128> BestBasis_R0_new(BestBasis_Ri); // i-th Basis represented in R0
  unsigned int i=0;

  for(auto& Op_Ri:BestBasis_Ri)
  {
    BestBasis_R0_new[i].bin = UpdateOp_inR0(BestBasis_R0_old, Op_Ri.bin);
    i++;
  }

  return BestBasis_R0_new;
}

/******************************************************************************/
/****************     Search in DIFFERENT REPRESENTATIONS   *******************/
/******************************************************************************/
vector<pair<__int128_t, unsigned int>> build_Kvect(vector<pair<__int128_t, unsigned int>> Nvect, list<__int128_t> Basis);

vector<Operator128> BestBasisSearch_Final(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max = 2, bool bool_print = false, unsigned int m_max=1000)
{
    k_max = (k_max<2)?2:k_max;  // k_max must be at least 2;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  SEARCH IN THE ORIGINAL BASIS:  *******************************";
    cout << endl << "*******************************************************************************************" << endl;

    unsigned int R_it = 0;   // Initial Representation --> R0

    vector<Operator128> BestBasis_R0 = BestBasisSearch_FixedRepresentation(Nvect, n, N, k_max, bool_print, R_it, m_max);

//Save Basis:
    string Basis_filename = OUTPUT_directory + "All_Bases_inRi.dat";
    fstream Basis_file(Basis_filename, ios::out); 

    Basis_file << "### File containing all the successive Basis" << endl;
    Basis_file << "### Note that Bases are given in the successive representation, and not in the original representation" << endl << endl;

    string Basis_filename_R0 = OUTPUT_directory + "All_Bases_inR0.dat";
    fstream Basis_file_R0(Basis_filename_R0, ios::out); 

    Basis_file_R0 << "### File containing all the successive Basis" << endl;
    Basis_file_R0 << "### The Bases are given in the original representation of the data" << endl << endl;

    PrintTerm_OpBasis(BestBasis_R0, n, N);   
    SaveFile_Basis(BestBasis_R0, n, Basis_file);
    SaveFile_Basis(BestBasis_R0, n, Basis_file_R0);

    Histo_BasisOpOrder(BestBasis_R0);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  SEARCH IN THE NEW BASIS:  ************************************";
    cout << endl << "*******************  Stops when the found basis is Identity  ******************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool isBasisIdentity = Check_Basis_Identity(BestBasis_R0);

    if ( !isBasisIdentity )
    {
      // **** BestBasis_Ri = store the Best Basis in the current representation Ri
      // **** BestBasis_R0 = store the Best Basis in the original representation R0;
      vector<Operator128> BestBasis_Ri(BestBasis_R0);
      vector<pair<__int128_t, unsigned int>> Kvect(Nvect); // Kvect = data in the current representation
      list<__int128_t> Basis_li;

      while( !isBasisIdentity ) // if the best basis is not the identity: then continue changing representation
      {
        cout << "-->> Change the representation of the data in the current Best Basis:" << endl;
        Basis_li.clear();
        for(auto& Op:BestBasis_Ri)  { Basis_li.push_back(Op.bin);  }  // extract the integer representation of the basis operators:
    
        Kvect = build_Kvect(Kvect, Basis_li);

        R_it += 1;   // New basis

        BestBasis_Ri.clear();
        BestBasis_Ri = BestBasisSearch_FixedRepresentation(Kvect, n, N, k_max, bool_print, R_it, m_max);

        PrintTerm_OpBasis(BestBasis_Ri, n, N);  
        SaveFile_Basis(BestBasis_Ri, n, Basis_file);
        Histo_BasisOpOrder(BestBasis_Ri);

        isBasisIdentity = Check_Basis_Identity(BestBasis_Ri);

        cout << "Check Basis Identity, Iteration = " << R_it << " : " << Check_Basis_Identity(BestBasis_Ri) << endl << endl; 

        BestBasis_R0 = UpdateBasis_inR0(BestBasis_R0, BestBasis_Ri);
        SaveFile_Basis(BestBasis_R0, n, Basis_file_R0);
      }
    }

    Basis_file.close();

    cout << "-->> All successive Bases are saved in the file: \'" <<  Basis_filename << "\'" << endl;
    cout << "Note that Bases are given in the successive representation, and not in the original representation" << endl << endl;

    return BestBasis_R0;
}





