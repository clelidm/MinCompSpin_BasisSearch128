#include <iostream>
#include <set>
#include <list>
#include <vector>
#include <cmath>

using namespace std;


/******************************************************************************/
/**********************     CONSTANTS  and TOOLS  *****************************/
/******************************************************************************/
#include "data.h"
const __int128_t one128 = 1;

//std::string int_to_bstring(__int128_t bool_nb, unsigned int n);
//void int_to_digits(__int128_t bool_nb, unsigned int n);

unsigned int min(unsigned int  a, unsigned int b)
{
  return !(b<a)?a:b;
}

/**************************************************************************************************************************************************/
/************************************************************   GAUSSIAN ELIMINATION:    **********************************************************/
/*************************************************   Return LEAD positions for basis selection   **************************************************/
/**************************************************************************************************************************************************/

/******************************************************************************/
/**********************   Using REF on Matrix over F2   ***********************/
/******************************************************************************/
void print_matrice(bool** M, int n, int m)
{
  for (int i=0; i<n; i++)
    { 
      for (int j=0; j<m; j++)  {  cout << M[i][j]; }
        cout << endl;
    }
}

/********************************************************************/
/*********    OPERATIONS for GAUSSIAN ELIMINATION on F2   ***********/
/********************************************************************/
void swap_row(bool** M, int i1, int i2, int n, int m)   //swap L_i1 and L_i2  in the matrix M
{
  if (i1>=n || i2>=n) { cout << "error swap" << endl; }
  else {
    bool* temp = M[i1]; //(bool*) malloc(m*sizeof(bool));
    M[i1]=M[i2];  M[i2]=temp;
    }
}

void add_row(bool** M, int i1, int i2, int n, int m)   //L_i2 <---- L_i1 XOR L_i2
{
  if (i1>=n || i2>=n) { cout << "error add" << endl; }
  else {
    for (int k = 0; k < m; k++)   {   M[i2][k] = ( M[i2][k] != M[i1][k] );   }
    }
}

/********************************************************************/
/**********************   Matrice REF:    ***************************/
/*********   return LEAD positions for basis selection   ************/
/********************************************************************/
// The matrix M is put in Reduced Row Echelon Form:
// And
// Returns a list of the column index where there is a lead position ==>> set of independent operators starting from the left

list<unsigned int> RREF_F2(bool** M, int n, int m)    // (i_lead, j_lead) = positions of the lead
{
  int j_lead = 0;
  int rank=0;    //rank = final number of leads
  bool test_stop=false;

  list<unsigned int> lead_positions;  // list of the column indices in which there is a lead position

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < m; i_lead++)   // 
    {
      //cout << "j_lead = " << j_lead << endl;
    i = i_lead;

    //search for the 1rst non-zero element in the column
    while (! M[i][j_lead])  // while M[i][j_lead] = 0
      {
        //cout << "test: " << M[i][j_lead] << endl;
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
        {
        //Record No Lead positions:
        //no_lead_positions.push_back(j_lead);
        j_lead++;   //1. go to the next column
          if (j_lead >= m)  { //cout << "break j_lead = " << j_lead << endl; 
                            test_stop=true; break; }  // reduction finished
        i = i_lead; //2. re-start the search for non-zero element from new (i_lead, j_lead)
        }
      }
    if(test_stop)   {  break;  }
    //cout << "real j_lead =" << j_lead << endl;
    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead
    rank++; 
    if (i != i_lead)
      { swap_row(M, i_lead, i, n, m); }

    //put to zero every element under the lead (starting from i = i_lead + 1)
    //for (i=0; i < n; i++)
    for (i=i_lead+1; i < n; i++)
      {
        //if (M[i][j_lead] && i!=i_lead)
        if (M[i][j_lead])
          { add_row(M, i_lead, i, n, m); }  //L_i <---- L_i XOR L_i_lead
      }

    // Record lead position:
    lead_positions.push_back(j_lead);  //cout << i << " ";

    //go to the next column
    j_lead++;
    }

  //for (i=j_lead; i<m; i++)
  //  {
  //    no_lead_positions.push_back(i);  //cout << i << " ";
  //  }

  cout << "\t Final rank = " << rank << endl; 
  return lead_positions;
}

/**************************************************************************************************************************************************/
/**************************************************************************************************************************************************/
/************************************************************   FIND BEST BASIS:    ***************************************************************/
/**************************************************************************************************************************************************/
/**************************************************************************************************************************************************/

/******************************************************************************/
/**********************   Convert  OpSet  to  F2 Matrix   *********************/
/******************************************************************************/
// This function returns a binary matrix that has the 'm' first operators of 'OpSet' as columns

// Each Operator is a column of the matrix
// n = Number of spins = number of rows --> 1rst index          //     !! We placed the lowest bit (most to the right) in the top row !!
// m = Number of selected operators = number of columns --> 2nd index

// This function place the 'm' first operators as column in the matrix:
MatrixF2 OpSet_to_MatrixF2(set<Operator128> OpSet, unsigned int m, unsigned int n)
{
  __int128_t Op = 0;

  // Adjust nb of columns if needed: cannot take more operators than there are in OpSet
  m = min(m, OpSet.size());

  // Create a Boolean Matrix:
  MatrixF2 Mat(n, m);

  unsigned int i = 0; //iteration over the row;
  unsigned int j = 0; //iteration over the columns;

  // Copy the m-first operators:
  set<Operator128>::iterator it_Op = OpSet.begin();

  //cout << "Test Fill: " << endl;

  while(j<m && it_Op != OpSet.end())
  {
    Op = (*it_Op).bin;
    //cout << "\t " << j << "\t" << int_to_bstring(Op, n) << "\t"; // << endl; //<< "\t" << Op_bin
    //int_to_digits(Op, n); 

    // filling in each column:
    for (i=0; i<n; i++)
    { 
      Mat.M[i][j] = Op & one128;
      Op >>= 1;
    }

    j++; // next column
    it_Op++;
  }

  Mat.OpSet_offset = j; // Number of operators analysed so far

  //cout << "Next Operator: " << endl;
  //Op = (*it_Op).bin;
  //cout << "\t " << j << "\t" << int_to_bstring(Op, n) << "\t"; // << endl; //<< "\t" << Op_bin
  //int_to_digits(Op, n); 

  return Mat;
}

/******************************************************************************/
/**********************   Convert  OpSet  to  F2 Matrix   *********************/
/*********************************   REFILL   *********************************/
/******************************************************************************/
void OpSet_to_MatrixF2_Refill(MatrixF2 *Mat, vector<Operator128> BestBasis, set<Operator128> OpSet, unsigned int m)
// Reminder: Each Operator is a column of the matrix
// m = Number of selected operators = number of columns --> 2nd index
// n = Number of spins = number of rows --> 1rst index
// !! We placed the lowest bit (most to the right) in the top row !!
{
  __int128_t Op_bin = 0;

  unsigned int i = 0; //iteration over the row;s
  unsigned int j = 0; //iteration over the columns;

  unsigned int n = (*Mat).n;   //total number of rows;

//  cout << "Test Refill -- First: " << endl; 

// Refill the lead Operators to the left of the Matrix 'M':
  for (auto& it_Op : BestBasis) 
  {
    Op_bin = it_Op.bin;
    //cout << "\t " << j << "\t" << int_to_bstring(Op_bin, n) << "\t";  //<< "\t" << Op_bin
    //int_to_digits(Op_bin, n);

    // filling in each column:
    for (i=0; i<n; i++)
    { 
      (*Mat).M[i][j] = Op_bin & one128;
      Op_bin >>= 1;
    }
    j++; // next column
  }

// Fill M with the next operators from OpSet, until there is m operators in M or until OpSet is empty
  set<Operator128>::iterator it_Op = OpSet.begin();
  advance(it_Op, (*Mat).OpSet_offset);

  //cout << "Test Refill -- Second: " << endl; 

  while(j<m && it_Op != OpSet.end())
  {
    Op_bin = (*it_Op).bin;
    //cout << "\t " << j << "\t" << int_to_bstring(Op_bin, n) << "\t";  //<< "\t" << Op_bin
    //int_to_digits(Op_bin, n);  

    // filling in each column:
    for (i=0; i<n; i++)
    { 
      (*Mat).M[i][j] = Op_bin & one128;
      Op_bin >>= 1;
    }

    j++; // next column
    it_Op++;
  }
  (*Mat).m = j;
  (*Mat).OpSet_offset += j - BestBasis.size();

// Fill rest of M with zeros if needed
  while(j<m)
  {
    for (i=0; i<n; i++)
      {   (*Mat).M[i][j] = 0;  }
    j++;
  }
}


/********************************************************************/
/******************    Print Terminal list lead   *******************/
/********************************************************************/
void PrintTerm_listLeads(list<unsigned int> lead_positions)
{
  for (auto& it_lead : lead_positions)
    {   cout << it_lead << endl;   }
  cout << endl;
}

/********************************************************************/
/******************    Extract Lead Operators    ********************/
/********************************************************************/
void Extract_LeadOp(set<Operator128> OpSet, unsigned int n, list<unsigned int> lead_positions, Struct_LowerBound* LowerBound, vector<Operator128>& BestBasis, unsigned int offset=0)
{
  //cout << "-->> Extract Leads" << endl;

// Current Basis size:  --> elements already saved
  unsigned int r = BestBasis.size();

// Skip the first 'r' leads that were already filled in:
  auto it_lead = lead_positions.begin();
  advance(it_lead, r);
  unsigned int lead_index = r;    // index of the last already extracted lead (this would be '0' if the current basis is empty)

// Skip the operators that were previously analyzed with 'offset': 
  auto it_Op = OpSet.begin();
  advance(it_Op, offset);  

// Extract the first next operators using lead_positions with offset: 
  while(it_lead != lead_positions.end())
  {
    advance(it_Op, (*it_lead)-lead_index);
    BestBasis.push_back(*it_Op);
    lead_index = (*it_lead);
    
    //cout << (*it_lead) << "\t" << int_to_bstring((*it_Op).bin, n) << "\t" << (*it_Op).bias << "\t Indices = "; //<< endl;
    //int_to_digits((*it_Op).bin, n);
    it_lead++;
  }

// Record information on the bias of the last basis operators:
  it_lead--;
  (*LowerBound).Bias = (*it_Op).bias;
  (*LowerBound).Index = offset + (*it_lead)-r;

  //cout << endl;

// Check that the size of the basis is the same as the size of the list of leads:
  if(BestBasis.size() != lead_positions.size())
    { cout << "Error in function 'Extract_LeadOp':   BestBasis.size() != lead_positions.size()" << endl << endl;}
  //cout << endl;
}


/********************************************************************/
/******************   Find Best Basis REF:    ***********************/
/********************************************************************/
/***********    Run multiple times lead search until     ************/
/**********   n independent operators have been found     ***********/
/***************      or  the OpSet is empty      *******************/
/********************************************************************/

/// Important: 'm' must be larger than the number of variables 'n' to find a basis;
///            For 'm < n', the function will look for the m first independent operators

vector<Operator128> BestBasis_inOpSet(set<Operator128> OpSet, unsigned int n, Struct_LowerBound* LowerBound, unsigned int m=1000) 
{
  vector<Operator128> BestBasis;

// ***** Initial search for best basis: ******************************
  MatrixF2 Mat = OpSet_to_MatrixF2(OpSet, m, n);

  cout << "Total number of Operators to analyse = " << OpSet.size() << endl << endl;
  cout << "-->> Search for the Best Basis by step of 'm' = " << Mat.m << " most biased operators:" << endl;

  cout << "\t Nit = " << 0 << ": \t";
  list<unsigned int> list_lead = RREF_F2(Mat.M, Mat.n, Mat.m);

  if (list_lead.size() > 0)
    {  Extract_LeadOp(OpSet, n, list_lead, LowerBound, BestBasis); }

// ***** Refill and search for best basis: ****************************
  unsigned int Nit = 1;                      // Iteration index (starting from 1)
  unsigned int offset = Mat.OpSet_offset;    // Number of operators already analysed in OpSet

  while( (offset < OpSet.size()) && BestBasis.size()<n && BestBasis.size()<m )
  {
    cout << "\t Nit = " << Nit << ": \t"; // "current Basis Size = " << BestBasis.size() << endl;
   
    OpSet_to_MatrixF2_Refill(&Mat, BestBasis, OpSet, m);
    list_lead = RREF_F2(Mat.M, Mat.n, Mat.m);

    if (list_lead.size() > BestBasis.size())
      {  Extract_LeadOp(OpSet, n, list_lead, LowerBound, BestBasis, offset); }


    offset = Mat.OpSet_offset; // New offset
    Nit++;
  }
  cout << endl << "-->> The Final Basis found has " << BestBasis.size() << " independent operators:" << endl;

  if(BestBasis.size()==n) {
    cout << "\t --> this is equal to the number \'n\' of variables: i.e., this is a Basis for the n-dimensional system" << endl;
  }
  else if(m<n && BestBasis.size() == m) {
    cout << "\t --> this is equal to the number 'm' = " << m << " provided for the analysis (which is smaller than the number of variables 'n' = " << n << ");" << endl;
    cout << "\t\t the Basis has the largest number of elements for an m-dimensional system." << endl;
  }
  else {
    cout << "\t --> The basis found has a dimension smaller than the dimension of the system analysed (< n and < m)" << endl;
  }
  cout << "\t --> Smallest Bias among the basis components = " << (*LowerBound).Bias << endl;
  cout << endl;

  return BestBasis;
}



/**************************************************************************************************************************************************/
/**************************************************************************************************************************************************/
/************************************************************    INVERT BASIS:    *****************************************************************/
/**************************************************************************************************************************************************/
/**************************************************************************************************************************************************/

/********************************************************************/
/**********************   Matrice REF:    ***************************/
/********************************************************************/
// The matrix M is put in Reduced Row Echelon Form
// And
// Returns the inverse of M

pair<int, bool**> RREF_F2_invert(bool** M, int n)    // (i_lead, j_lead) = positions of the lead
{
  //cout << "Initial matrix:  " << endl;
  //print_matrice(M, n, n);
  //cout << endl;

  // Create a Identity Matrix:
  bool** M_id = (bool**) malloc(n*sizeof(bool*));  // n rows --> 1rst index
  for (int i=0; i<n; i++)  // n columns --> 2nd index
  {   
    M_id[i] = (bool*) malloc(n * sizeof(bool));  
    for (int j=0; j<n; j++)  // n columns --> 2nd index
    {
      M_id[i][j] = 0;
    } 
    M_id[i][i] = 1;
  }

//  cout << "Identity matrix:  " << endl;
//  print_matrice(M_id, n, n);

  // Reduction:
  int j_lead = 0;
  int rank=0;    //rank = final number of leads
  bool test_stop=false;

  //list<unsigned int> lead_positions;  // list of the column indices in which there is a lead position

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < n; i_lead++)   // 
  {
    //cout << "j_lead = " << j_lead << endl;
    i = i_lead;

    //search for the 1rst non-zero element in the column
    while (! M[i][j_lead])  // while M[i][j_lead] = 0
    {
        //cout << "test: " << M[i][j_lead] << endl;
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
      {
        j_lead++;   //1. go to the next column
        if (j_lead >= n)  
        { 
          test_stop=true; break; 
        }  // reduction finished
        i = i_lead; //2. re-start the search for non-zero element from new (i_lead, j_lead)
      }
    }
    if(test_stop)   {  break;  }
    //cout << "real j_lead =" << j_lead << endl;

    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead
    rank++; 
    if (i != i_lead)
    { 
      //cout << "i_lead = " << i_lead << " and i = " << i << ": swap rows:" << endl;
      swap_row(M, i_lead, i, n, n); 
      swap_row(M_id, i_lead, i, n, n); 
      //print_matrice(M_id, n, n);
      //cout << endl;
    }

    //put to zero all the element with i != i_lead
    for (i=0; i < n; i++)
    {
      if (M[i][j_lead] && i!=i_lead)  //L_i <---- L_i XOR L_i_lead
      { 
        add_row(M, i_lead, i, n, n); 
        add_row(M_id, i_lead, i, n, n); 

        //cout << "Add row: i_lead = " << i_lead << " and i = " << i << " into i" << endl;
        //print_matrice(M_id, n, n);
        //cout << endl;
      }
    }

    // Record lead position:
    //lead_positions.push_back(j_lead);  //cout << i << " ";

    //go to the next column
    j_lead++;
  }

  //for (i=j_lead; i<m; i++)
  //  {
  //    no_lead_positions.push_back(i);  //cout << i << " ";
  //  }   

  return make_pair(rank, M_id);
}





