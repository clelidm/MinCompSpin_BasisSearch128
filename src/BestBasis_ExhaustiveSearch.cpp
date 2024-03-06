#include <iostream>
//#include <sstream>
//#include <fstream>

#include <set>
#include <vector>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;

#include "data.h"

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t one128 = 1;

/******************************************************************************/
/***********************   All Operators with 1 bit only  *********************/
/************************   Find the lowest bias value  ***********************/
/******************************************************************************/

set<Operator128> All_Op_k1(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, double *lowest_bias, bool print = false);

Operator128 Value_Op(__int128_t Op_bin, vector<pair<__int128_t, unsigned int>> Nvect, double Nd);

vector<Operator128> BestBasis_inOpSet(set<Operator128> OpSet, unsigned int n, Struct_LowerBound* LowerBound, unsigned int m=1000);

/******************************************************************************/
/*******************************   All Operators  *****************************/
/****************   Keep only the one with bias larger than LB  ***************/
/******************************************************************************/

set<Operator128> All_Op_LBk1 (vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool print = false)
{
  double lowest_bias = 0;

  set<Operator128> OpSet = All_Op_k1(Nvect, n, N, &lowest_bias, print);
  Operator128 Op;
  double Nd = (double) N;

  cout << "-->> Compute ALL the (2^n-1) Operators" << endl;
  cout << "     Rank the operators with bias larger than lower bound (fixed by the least informative first order operator):" << endl;

  unsigned int Op_bin_max =  (one128 << n) - 1;

  for (unsigned int Op_bin = 1; Op_bin <= Op_bin_max; Op_bin++)
  {
    Op = Value_Op(Op_bin, Nvect, Nd);
    if (Op.bias > lowest_bias) { OpSet.insert(Op); } 
  }

  return OpSet;
}

/******************************************************************************/
/***************************   Exhaustive Search  *****************************/
/******************************************************************************/

vector<Operator128> BestBasis_ExhaustiveSearch(vector<pair<__int128_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool bool_print = false)
{
  auto start = chrono::system_clock::now();

  cout << endl << "*******************************************************************************************";
  cout << endl << "************************  EXHAUSTIVE SEARCH FOR THE BEST BASIS:  **************************";
  cout << endl << "*******************************************************************************************" << endl;

//  cout << "-->> Compute all Operators, with a smallest accepted biased fixed by the least informative first order operator:" << endl;
  cout << endl;
  set<Operator128> OpSet = All_Op_LBk1 (Nvect, n, N, bool_print);

// Time:
  auto end = chrono::system_clock::now();  chrono::duration<double> elapsed = end - start;
  cout << "\t Elapsed time (in s): " << elapsed.count() << endl << endl; 

  cout << "-->> Search for the Best Basis:\t\t"; // << endl;

  Struct_LowerBound LB; // Lower Bound Info
  LB.Bias = 0;
  vector<Operator128> BestBasis = BestBasis_inOpSet(OpSet, n, &LB, 1000); // LB will be over-written with the updated values

// Time:
  end = chrono::system_clock::now();  elapsed = end - start;
  cout << "Total elapsed time (in s): " << elapsed.count() << "\t for Exhaustive Search" << endl << endl; 

  return BestBasis;
}