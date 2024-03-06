#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <list>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
#include "data.h"
const __int128_t one128 = 1;

unsigned int bitset_count(__int128_t bool_nb);

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/

vector<pair<__int128_t, unsigned int>> read_datafile128_vect(string datafilename, unsigned int *N, unsigned int r)    // O(N)  where N = data set size
{
  auto start = chrono::system_clock::now();

  cout << endl << "--->> Read the datafile: \"" << datafilename << "\", \t Build Nset..." << endl;
  cout << "\t Number of variables to read: n = " << r << endl;

  string line, line2;     char c = '1';
  __int128_t state = 0, Op;
  (*N) = 0;            // N = dataset sizes

// ***** The data is stored in Nset as an histogram:  ********************************
  map<__int128_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

  ifstream myfile (datafilename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,r);          //take the r first characters of line
      Op = one128 << (r - 1);
      state = 0;
      for (auto &elem: line2)     //convert string line2 into a binary integer
      {
        if (elem == c) { state += Op; }
        Op = Op >> 1;
      }
      Nset[state] += 1;
      (*N)++;
    }
    myfile.close();
  }
  else cout << endl << "--->> Unable to open file: Check datafilename and location." << endl << endl;

  if ((*N) == 0) 
    { 
    cout << endl << "--->> Failure to read the file, or file is empty:  Terminate." << endl << endl;
    }
  else
    {
    cout << endl << "--->> File has been read successfully:" << endl;
    cout << "\t Data size, N = " << (*N) << endl;
    cout << "\t Number of different states, Nset.size() = " << Nset.size() << endl << endl;
    }

  vector<pair<__int128_t, unsigned int>> Nvect(Nset.size());
  int i=0;
  for (auto& my_pair : Nset)
  {
    Nvect[i]=my_pair;
    i++;
  }

  auto end = chrono::system_clock::now();  
  chrono::duration<double>  elapsed = end - start;
  cout << endl << "Elapsed time (in s): " << elapsed.count() << endl << endl;  

  return Nvect;
}

/****************    PRINT Nset in file:    ************************/
/*void read_Nset (map<uint32_t, unsigned int> Nset, unsigned int N, string OUTPUTfilename)
// map.second = nb of time that the state map.first appears in the data set
{
  map<uint32_t, unsigned int>::iterator it;
  int Ncontrol = 0;

  fstream file(OUTPUTfilename.c_str(), ios::out);
  file << "#N = " << N << endl;
  file << "#Total number of accessible states = " << NOp_tot << endl;
  file << "#Number of visited states, Nset.size() = " << Nset.size() << endl;
  file << "#" << endl;
  file << "#1: state \t #2: nb of pts in state \t #3: Pba state" << endl;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    file << it->first << ":\t" << bitset<n>(it->first) << " => " << it->second; // << endl;
    file << "  \t  P = " << it->second / (float) N << endl;
    Ncontrol += it->second;
  }

  if (Ncontrol != N) { cout << "Error function \'read_Nset\': Ncontrol != N" << endl;  }

  file.close();
}
*/

/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a basis (defined by the m-basis list) --> returns the new m-state (i.e. state in the new m-basis)
// Rem: must have m <= n 

// Basis operators are ordered from Right to Left,
// i.e., the first basis operator Phi_0 corresponds to the lowest bit (rightmost bit)
// the last basis operator Phi_{n-1} corresponds to the n-th highest bit

// mu = old state
// final_mu = new state

__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis)
{
  __int128_t un_i = 1, proj;
  __int128_t final_mu = 0;

  list<__int128_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    proj = (*phi_i) & mu;
    if ( (bitset_count(proj) % 2) == 1) // odd number of 1, i.e. sig_i = 1
    {
      final_mu += un_i;
    }
    un_i = (un_i << 1);
  }

  return final_mu;
}


/******************************************************************************/
/******************** CHANGE of BASIS: build K_SET ****************************/
/******************************************************************************/
// Build Kvect for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:

vector<pair<__int128_t, unsigned int>> build_Kvect(vector<pair<__int128_t, unsigned int>> Nvect, list<__int128_t> Basis)
// sig_m = sig in the new basis and cut on the m first spins 
// Kvect[sig_m] = #of time state mu_m appears in the data set
{
    map<__int128_t, unsigned int > Kvect_map;
    __int128_t sig_m;    // transformed state and to the m first spins

// ***** Build Kvect: *************************************************************************************
    cout << endl << "--->> Build Kvect..." << endl;
    cout << "## Basis elements are ordered from the right (s_1) to the left (s_n)." << endl;

    for (auto const& it : Nvect)
    {
        sig_m = transform_mu_basis((it).first, Basis); // transform the initial state s=(it).first into the new basis
        Kvect_map[sig_m] += ((it).second); // ks = (it).second = number of time state s appear in the dataset
    }
    cout << endl;

// ***** Convert map to a vector:  for faster reading later on ********************************************
    vector<pair<__int128_t, unsigned int>> Kvect(Kvect_map.size());

    int i=0;
    for (auto& my_pair : Kvect_map)
    {
        Kvect[i]=my_pair;
        i++;
    }

    cout << "\t Kvect.size() = " << Kvect.size() << endl;

    return Kvect;
}

/******************************************************************************/
/*********************   TRANSFORM DATASET to a NEW BASIS    ******************/
/******************************************************************************/
string filename_remove_extension(string filename);
string int_to_bstring(__int128_t bool_nb, unsigned int n);

void convert_datafile_to_NewBasis(string input_dir, string datafilename, unsigned int r, vector<Operator128> BestBasis_vect)    // O(N)  where N = data set size
{
  auto start = chrono::system_clock::now();

  list<__int128_t> Basis_li;
  for(auto& Op:BestBasis_vect)  { Basis_li.push_back(Op.bin);  }  // extract the integer representation of the basis operators:

  cout << endl << "--->> Read the datafile: \"" << (input_dir + datafilename) << "\"" << endl;
  cout << "\t Number of variables to read: n = " << r << endl;

  string New_datafilename = OUTPUT_directory + filename_remove_extension(datafilename) + "_inBestBasis.dat";
  fstream file_newdata(New_datafilename, ios::out);

  cout << endl << "--->> Transform to new basis...";
  cout << endl << "\t Write the new dataset in the file: \"" << New_datafilename << "\"" << endl;

  string line, line2;     char c = '1';
  __int128_t state = 0, state_new = 0, Op;
  unsigned int N = 0;            // N = dataset sizes

// ***** Read the original data and convert in new basis:  ********************************

  ifstream file_data ((input_dir + datafilename).c_str());
  if (file_data.is_open())
  {
    while (getline (file_data,line))
    {
      line2 = line.substr (0,r);          //take the r first characters of line
      Op = one128 << (r - 1);
      state = 0;
      for (auto &elem: line2)     //convert string line2 into a binary integer
      {
        if (elem == c) { state += Op; }
        Op = Op >> 1;
      }
      state_new = transform_mu_basis(state, Basis_li);
      file_newdata << int_to_bstring(state_new, r) << endl;
      //int_to_digits_file(state_new, r, file_newdata);
      N++;
    }
    file_data.close();
    file_newdata.close();
  }
  else cout << endl << "--->> Unable to open file: Check datafilename and location." << endl << endl;

  if (N == 0) 
    { 
    cout << endl << "--->> Failure to read the file, or file is empty:  Terminate." << endl << endl;
    }
  else
    {
    cout << endl << "--->> File has been converted to the new basis successfully:" << endl;
    cout << "\t Data size, N = " << N << endl << endl;

    cout << "Operators are ordered from Right to Left: i.e. " << endl; 
    cout << " \t - the first basis operator corresponds to the rightmost bit (lowest bit)" << endl;
    cout << " \t - the n-th basis operator corresponds to the leftmost bit (highest bit) " << endl;
    }

  auto end = chrono::system_clock::now();  
  chrono::duration<double>  elapsed = end - start;
  cout << endl << "Elapsed time (in s): " << elapsed.count() << "\t for converting data" << endl << endl;  
}


