#include <iostream>
#include <iomanip>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;


/******************************************************************************/
/**********************     CONSTANTS  and TOOLS  *****************************/
/******************************************************************************/
#include "data.h"
const __int128_t one128 = 1;

string int_to_bstring(__int128_t bool_nb, unsigned int n);
void int_to_digits(__int128_t bool_nb, unsigned int n);
void int_to_digits_file(__int128_t bool_nb, unsigned int r, fstream &file);

unsigned int bitset_count(__int128_t bool_nb);

/******************************************************************************/
/*******************   Convert  Basis  to  F2 Matrix   ************************/
/******************************************************************************/
list<unsigned int> RREF_F2(bool** M, int n, int m);

// This function returns a binary matrix that has the operators of the basis 'Basis' as columns:
// Each Operator is a column of the matrix
// n = Number of spins = number of rows --> 1rst index          //     !! We placed the lowest bit (most to the right) in the top row !!
// m = Number of basis operators = number of columns --> 2nd index

bool** Basis_to_MatrixF2(vector<Operator128> Basis, unsigned int n)
{
  unsigned int m = Basis.size();
  __int128_t Op = 0;

  // Create a Boolean Matrix:
  bool** M = (bool**) malloc(n*sizeof(bool*));  // n rows --> 1rst index
  for (int i=0; i<n; i++)
    {   M[i] = (bool*) malloc(m * sizeof(bool));  }  // m columns --> 2nd index

  unsigned int i = 0; //iteration over the n row;
  unsigned int j = 0; //iteration over the m columns;

  // Copy the basis operators:
  for (auto& it_Op : Basis)
  {
    Op = it_Op.bin;

    // *****  Filling in each column with a basis operator: ****************
      // -- Last bit (rightmost) at the bottom; first bit (leftmost) at the top.
      // -- First basis operator to the Left (i.e. placed first in the matrix) 
      //    Last basis operator to the Right (placed last)
      //       --> note that this is the opposite of the convention adopted in writing a state in the new basis
      //           (there sig_1 = Rightmost bit)
      // --> One must be careful when re-invert:Matrix to Basis !!!

    for (i=0; i<n; i++)
    { 
      M[(n-1-i)][j] = Op & one128;
      Op >>= 1;
    }

    j++; // next column
  } 

  return M;
}

/******************************************************************************/
/************   Check if a set of operators are all independent   *************/
/******************************************************************************/
// return 'True' if all the element in the matrix are independent

bool Is_Basis(vector<Operator128> Basis, unsigned int n)
{
  cout << "-->> Check if the set of operators are independent:" << endl;
  cout << "\t Number of operators analysed: " << Basis.size() << endl;
  // Convert the basis to a boolean matrix:
  bool** M_Basis = Basis_to_MatrixF2(Basis, n);

  // Row Reduction procedure:
  list<unsigned int> list_lead = RREF_F2(M_Basis, n, Basis.size());

  // Free memory:
  for (int i=0; i<n; i++)
    {   free(M_Basis[i]);   } 
  free(M_Basis);

  // Return if the operators are independent or not:
  if (list_lead.size() == Basis.size())
    { 
      cout << "All the operators are independents." << endl << endl;
      return true;  
    }

  else
    { 
      cout << "Not all the operators are independents:";
      cout << "   (number of independent operators = " << list_lead.size() << ") < (number of operators = " << Basis.size() << " ) " << endl << endl;
      return false; 
    }
}

/******************************************************************************/
/**********   INVERT a BASIS: RETURN INVERSE GAUGE TRANSFORMATION   ***********/
/******************************************************************************/
vector<Operator128> MatrixF2_to_Basis(bool** M, unsigned int n)
{
  vector<Operator128> Basis;
  Operator128 Op;
  Op.r = n;
  Op.k1 = 0;

  int i = 0; //iteration over the n row;
  int j = 0; //iteration over the m columns;

  __int128_t Op_bin = 1, state = 0;

  // Copy the basis operators from M to Basis:
  for (j=(n-1); j>=0; j--) // Copying each column into an operator:
  // read operators from the right to the left!!!! 
  // as s1 = rightmost bit = first operator of the Basis_list_invert
  //    sn = left most bit
  { 
    Op_bin = 1; // sig_1 (for i=0) = the lowest bit // sig_n (for i=(n-1)) = the highest bit
    //Op_bin = one128 << (n - 1);
    state = 0;
    for (i=0; i<n; i++) // Turn a column to an operator:
    {
      // first element of the column = sig_1 (see function "Basis_to_MatrixF2") 
            // should be placed as the Righmost bit for the user (cf. convention adopted in writing a state in the new basis: with s_1=[lowest bit] )
      if(M[i][j] == 1)  { state += Op_bin; }
      //cout << "Op_bin = " << int_to_bstring(Op_bin, n) << endl; 
      Op_bin = Op_bin << 1;    
    } 
    //j++; // next column

    Op.bin = state;
    Basis.push_back(Op);
    //cout << int_to_bstring(state, n) << endl;
    //string int_to_bstring(__int128_t bool_nb, unsigned int n);
  }

  return Basis;
}

pair<int, bool**> RREF_F2_invert(bool** M, int n);
//bool** RREF_F2_invert(bool** M, int n); // only works for square matrices

// n = number of variables
// m = number of operators in the independent set
// if m > n : the set cannot be independent: stop the procedure;
// if n=m: check if rank=n, then everything is good and return the inverse basis
// if m < n:  even if it is an independent set, it is not a basis, and may not be invertable: stop the procedure.

vector<Operator128> Invert_Basis(vector<Operator128> Basis, unsigned int n)
{
  vector<Operator128> Basis_invert;

  if (Basis.size() == n)
  {
    bool** M_Basis = Basis_to_MatrixF2(Basis, n);

    // Row Reduction procedure to invert:
    pair<int, bool**> M_invert = RREF_F2_invert(M_Basis, n);

    if(M_invert.first == n) // M_invert.first = rank
    {
      cout << "Rank = n = " << n << "\t: this is a basis and can be inverted." << endl << endl;
      Basis_invert = MatrixF2_to_Basis(M_invert.second, n);
    }
    else 
    {
      cout << "The Rank = " << M_invert.first << " is smaller than the number of variables, n = " << n << "\t: this is not a basis." << endl << endl;
      //cout << "Note that the inverse transformation provided is therefore incomplete." << endl;
    } 

    // Free memory:
    for (int i=0; i<n; i++)
    {   
      free(M_Basis[i]);   free(M_invert.second[i]);   
    } 
    free(M_Basis);
    free(M_invert.second);  
  }
  else
  {
    cout << "The number of operators in the basis provided, m=" << Basis.size() << ", is not equal to the number of spin variables, n=" << n << ":" << endl << endl;
  }

  return Basis_invert;
}

/******************************************************************************/
/****************   Print Terminal Vector Best Operators  *********************/
/******************************************************************************/
void PrintTerm_Basis(vector<Operator128> Basis, unsigned int n, unsigned int N)
{
  int i = 1;
  double p1 = 1, LogLi = 0, LogL = 0, Nd = (double) N;

  cout << "-->> Print Basis Operators: \t Number of basis operators = " << Basis.size() << endl << endl;  
  cout << "## 1:i \t 2:bin \t\t 3:bias\t 4:N[Op_i=1] \t 5:p[Op_i=1] \t 6:<Op> \t 7:LogL[Op_i] \t 8:Op_index " << endl << "## " << endl; 

  for (auto& Op : Basis)
  {
    p1 = ((double) Op.k1) / Nd;
    LogLi = (p1!=0 && p1!=1)? p1*log(p1)+(1-p1)*log(1-p1) : 0;
    LogL += LogLi;

    cout << fixed;
    cout << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << setprecision(5) << Op.bias << " \t" << Op.k1  << "\t";
    cout << setprecision(6) << p1 << " \t" << 1-2*p1 << " \t" << LogLi << " \t Indices = "; 
    int_to_digits(Op.bin, n);
    i++;
  }
  cout << endl;
  cout << "##  LogL / N = " << LogL << endl;  //<< setprecision (6) 
  cout << "## -LogL / N / log(2) = " << -LogL/log(2.) << " bits per datapoints "<< endl;  //<< setprecision (6) 
  cout << endl;
}

void PrintFile_Basis(vector<Operator128> Basis, unsigned int n, unsigned int N, string filename)
{
  string OpSet_filename = OUTPUT_directory + filename + ".dat";

  cout << "-->> Print Basis Operators in the file: \'" <<  OpSet_filename << "\'" << endl;
  fstream file_OpBasis(OpSet_filename, ios::out);
  
  int i = 1;
  double p1 = 1, LogLi = 0, LogL = 0, Nd = (double) N;

  file_OpBasis << "## Basis: Total number of operators = " << Basis.size() << endl << "## " << endl; 
  file_OpBasis << "## 1:i \t 2:bin \t\t 3:bias\t 4:N[Op_i=1] \t 5:p[Op_i=1] \t 6:<Op> \t 7:LogL[Op_i] \t 8:Op_index " << endl << "## " << endl; 

  for (auto& Op : Basis)
  {
    p1 = ((double) Op.k1) / Nd;
    LogLi = (p1!=0 && p1!=1)? p1*log(p1)+(1-p1)*log(1-p1) : 0;
    LogL += LogLi;

    file_OpBasis << fixed;
    file_OpBasis << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << setprecision(5) << Op.bias << " \t" << Op.k1  << "\t";
    file_OpBasis << setprecision(6) << p1 << " \t" << 1-2*p1 << " \t" << LogLi << " \t Indices = "; 
    //file_OpBasis << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << Op.bias << "\t" << Op.k1  << "\t" << p1 << "\t" << LogLi << "\t Indices = "; 
    int_to_digits_file(Op.bin, n, file_OpBasis);
    i++;
  }
  file_OpBasis << endl;
  file_OpBasis << "##  LogL / N = " << LogL << endl;  //<< setprecision (6) 
  file_OpBasis << "## -LogL / N / log(2) = " << -LogL/log(2.) << " bits per datapoints "<< endl;  //<< setprecision (6) 
  file_OpBasis.close();

  cout << endl;
}

// ****** SHORT VERSIONS:
void PrintTerm_OpBasis_Short(vector<Operator128> Basis, unsigned int n)
{
  int i = 1;

  cout << "-->> Print Basis Operators: \t Number of basis operators = " << Basis.size() << endl << endl;  
  cout << "## 1:i \t 2:bin \t\t 3:Op_index " << endl << "## " << endl; 

  for (auto& Op : Basis)
  {
    cout << fixed;
    cout << i << "\t" << int_to_bstring(Op.bin, n);
    cout << " \t Indices = "; 
    int_to_digits(Op.bin, n);
    i++;
  }
  cout << endl;
}

void PrintFile_OpBasis_Short(vector<Operator128> Basis, unsigned int n, unsigned int N, string filename)
{
  string OpSet_filename = OUTPUT_directory + filename + ".dat";

  cout << "-->> Print Basis Operators in the file: \'" <<  OpSet_filename << "\'" << endl;
  fstream file_OpBasis(OpSet_filename, ios::out);
  
  int i = 1;

  file_OpBasis << "## Basis: Total number of operators = " << Basis.size() << endl << "## " << endl; 
  file_OpBasis << "## 1:bin \t 2:count \t 3:Op_index " << endl << "## " << endl; 

  for (auto& Op : Basis)
  {
    file_OpBasis << fixed;
    file_OpBasis << int_to_bstring(Op.bin, n) << "\t" << i << " \t Indices = "; 
    //file_OpBasis << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << Op.bias << "\t" << Op.k1  << "\t" << p1 << "\t" << LogLi << "\t Indices = "; 
    int_to_digits_file(Op.bin, n, file_OpBasis);
    i++;
  }
  file_OpBasis.close();

  cout << endl;
}

// ****** FINAL BASIS VERSION:
void PrintTerm_FinalBasis(vector<Operator128> Basis, unsigned int n, unsigned int N)
{
  int i = 1;
  double p1 = 1, LogLi = 0, LogL = 0, Nd = (double) N;

  cout << "-->> Print Basis Operators: \t Number of basis operators = " << Basis.size() << endl << endl;  
  cout << "## 1:i \t 2:bin \t\t 3:bias\t 4:N[Op_i=1] \t 5:p[Op_i=1] \t 6:<Op> \t 7:LogL[Op_i] \t 8:Op_index " << endl << "## " << endl; 

  for (auto& Op : Basis)
  {
    p1 = ((double) Op.k1) / Nd;
    LogLi = (p1!=0 && p1!=1)? p1*log(p1)+(1-p1)*log(1-p1) : 0;
    LogL += LogLi;

    cout << fixed;
    cout << "sig_" << setw(3) << setfill(' ') << left << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << setprecision(5) << Op.bias << " \t" << Op.k1  << "\t";
    cout << setprecision(6) << p1 << " \t" << 1-2*p1 << " \t" << LogLi << " \t Indices = "; 
    int_to_digits(Op.bin, n);
    i++;
  }
  cout << endl;
  cout << "##  LogL / N = " << LogL << endl;  //<< setprecision (6) 
  cout << "## -LogL / N / log(2) = " << -LogL/log(2.) << " bits per datapoints "<< endl;  //<< setprecision (6) 
  cout << endl;

  cout << "Convention for reading this basis:" << endl;
  cout << "## \t   bits are organised in the same order as in the original dataset," << endl;
  cout << "## \t   i.e. rightmost bit of an operator = rightmost bit in the data: labeled s_1 in the inverse basis below." << endl;
  cout << "## \t        leftmost  bit of an operator = leftmost  bit in the data: labeled s_n in the inverse basis below." << endl;
  cout << endl;
}


void PrintFile_FinalBasis(vector<Operator128> Basis, unsigned int n, unsigned int N, string filename)
{
  string OpSet_filename = OUTPUT_directory + filename + ".dat";

  cout << "-->> Print Basis Operators in the file: \'" <<  OpSet_filename << "\'" << endl;
  fstream file_OpBasis(OpSet_filename, ios::out);
  
  int i = 1;
  double p1 = 1, LogLi = 0, LogL = 0, Nd = (double) N;

  file_OpBasis << "## Basis: Total number of operators = " << Basis.size() << endl << "## " << endl; 
  file_OpBasis << "## 1:i \t 2:bin \t\t 3:bias\t 4:N[Op_i=1] \t 5:p[Op_i=1] \t 6:<Op> \t 7:LogL[Op_i] \t 8:Op_index " << endl << "## " << endl; 

  for (auto& Op : Basis)
  {
    p1 = ((double) Op.k1) / Nd;
    LogLi = (p1!=0 && p1!=1)? p1*log(p1)+(1-p1)*log(1-p1) : 0;
    LogL += LogLi;

    file_OpBasis << fixed;
    file_OpBasis << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << setprecision(5) << Op.bias << " \t" << Op.k1  << "\t";
    file_OpBasis << setprecision(6) << p1 << " \t" << 1-2*p1 << " \t" << LogLi << " \t Indices = "; 
    //file_OpBasis << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << Op.bias << "\t" << Op.k1  << "\t" << p1 << "\t" << LogLi << "\t Indices = "; 
    int_to_digits_file(Op.bin, n, file_OpBasis);
    i++;
  }
  file_OpBasis << endl;
  file_OpBasis << "##  LogL / N = " << LogL << endl;  //<< setprecision (6) 
  file_OpBasis << "## -LogL / N / log(2) = " << -LogL/log(2.) << " bits per datapoints "<< endl;  //<< setprecision (6) 
  file_OpBasis << endl;

  file_OpBasis << "## Convention for reading this basis:" << endl;
  file_OpBasis << "## \t   bits are organized in the same order as in the original dataset," << endl;
  file_OpBasis << "## \t   i.e. rightmost bit of an operator = rightmost bit in the data: labeled s_1 in the Inverse Basis file." << endl;
  file_OpBasis << "## \t        leftmost  bit of an operator = leftmost  bit in the data: labeled s_n in the Inverse Basis file." << endl;
  file_OpBasis << endl;

  file_OpBasis.close();
  cout << endl;
}

// ****** INVERT BASIS VERSIONS:
void PrintTerm_Basis_inverse(vector<Operator128> Basis, unsigned int r)
{
  cout << "-->> Print Operators of Inverse Basis: \t Number of basis operators = " << Basis.size() << endl << endl; 
  cout << "## 1:i \t 2:bin \t\t 3:Op_index " << endl << "## " << endl; 
  
  int i = 1;
  for (auto& Op : Basis)
  {
    cout << fixed;
    cout << "s_" << setw(3) << setfill(' ') << left << i << "\t" << int_to_bstring(Op.bin, r);
    cout << " \t Indices = "; 
    int_to_digits(Op.bin, r);
    i++;
  }
  cout << endl;

  cout << "Convention for reading this basis in comparison with the original basis of the data:" << endl;
  cout << "## \t   -- s_1 = Rightmost bit in the original data" << endl;
  cout << "## \t   -- s_n = Leftmost bit in the original data" << endl;
  cout << "##" << endl;
  cout << "## Besides, in the basis above:" << endl;
  cout << "## \t   -- Rightmost bit = first basis operator, labeled sig_1 above;" << endl;
  cout << "## \t   -- Leftmost bit  = last  basis operator, labeled sig_n above." << endl;
  cout << endl;
}

void PrintFile_Basis_inverse(vector<Operator128> Basis, unsigned int n, string filename)
{
  string OpSet_filename = OUTPUT_directory + filename + ".dat";

  cout << "-->> Print Operators of Inverse Basis in the file: \'" <<  OpSet_filename << "\'" << endl;
  fstream file_OpBasis(OpSet_filename, ios::out);
  
  int i = 1;

  file_OpBasis << "## Basis: Total number of operators = " << Basis.size() << endl << "## " << endl; 
  file_OpBasis << "## 1:count \t 2:bin \t 3:Op_index " << endl << "## " << endl; 

  for (auto& Op : Basis)
  {
    file_OpBasis << fixed;
    file_OpBasis << "s_" << setw(3) << setfill(' ') << left << i << " \t" << int_to_bstring(Op.bin, n) << " \t Indices = "; 
    //file_OpBasis << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << Op.bias << "\t" << Op.k1  << "\t" << p1 << "\t" << LogLi << "\t Indices = "; 
    int_to_digits_file(Op.bin, n, file_OpBasis);
    i++;
  }

  file_OpBasis << " " << endl;
  file_OpBasis << "## Convention for reading this basis in comparison with the original basis of the data:" << endl;
  file_OpBasis << "## \t   -- s_1 = Rightmost bit in the original data" << endl;
  file_OpBasis << "## \t   -- s_n = Leftmost bit in the original data" << endl;
  file_OpBasis << "##" << endl;
  file_OpBasis << "## Besides, in the basis above:" << endl;
  file_OpBasis << "## \t   -- Rightmost bit = first basis operator in the Best Basis file;" << endl;
  file_OpBasis << "## \t   -- Leftmost bit  = last  basis operator in the Best Basis file." << endl;

  file_OpBasis.close();
}


/******************************************************************************/
/*********************   Histo Order of basis Operators  **********************/
/******************************************************************************/
map<unsigned int, unsigned int> Histo_BasisOpOrder(vector<Operator128> Basis)
{
  map<unsigned int, unsigned int> histo_order;

  for (auto& Op : Basis) {
    histo_order[bitset_count(Op.bin)]+=1;
  }

  cout << "Number of basis operators of a given order k:" << endl;
  for (auto& histo: histo_order){
    cout << "\t k = " << histo.first << " :\t" << histo.second << " operators" << endl;
  }
  cout << endl;

  return histo_order;
}






