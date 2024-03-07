# MinCompSpin BasisSearch 128

# Search for the Best Independent Minimally Complex Models (MCM)

This program searches for the **Best Basis representation for a chosen binary dataset**, while taking into account all possible **high-order patterns of the data**. It was developed for the paper [*Statistical Inference of Minimally Complex Models*](https://arxiv.org/abs/2008.00520) [1]. More details on the general algorithm can be found in the paper.

The program can be used for datasets with up to $n=128$ random variables.

[1]  C. de Mulatier, P. P. Mazza, M. Marsili, *Statistical Inference of Minimally Complex Models*, [arXiv:2008.00520](https://arxiv.org/abs/2008.00520)

## General information

**Best Basis:** The Best basis of a binary data with `n` variables is the one for which the independent model formed by `n` field operators has the largest log-likelihood (and therefore the largest log-evidence, as all independent models with the same number of operators are equivalent -- see Ref[1]).

There are three main functions that you can use to **search for the best basis** from the `main.cpp`:

 1) **Exhaustive Search:** This function will compute all $2^n-1$ operators and will search for the best basis among them with a Greedy approach (i.e. rank them from the most to least biased and extract the set of the `n` most biased independent operators starting from the most biased one):
    ```c++
    vector<Operator64> BestBasis_ExhaustiveSearch(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool bool_print = false)
    ```
    Note: we advise doing such a search only for small systems (up to ~15 variables).

 2) **Search in a fixed representation up to order `kmax`:** This function searches for the best Basis among all operators up to order `kmax` in a given representation (which is the representation used when storing the data in `Nvect` -- by default, this is the original representation of the data):
    ```c++
    BestBasisSearch_FixedRepresentation(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, unsigned int B_it, bool bool_print = false)
    ```
    If you take the largest order to be equal to the number of variables (i.e., `kmax = n`), then this function will perform an exhaustive search for the best basis among all possible operators (exactly as the algorithm 1 just above).

 3) **Search in varying representations:** This function performs the search procedure described in Ref.[1]. The program first searches for the best basis up to order `k_max`; the data is then successively transformed in the representation given by the previously found best basis, and the program searches for the new best basis in this representation. The algorithm stops when the new basis found is the identity (i.e. the basis has not changed).
    ```c++
    vector<Operator64> BestBasisSearch_Final(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, bool bool_print = false)
    ```
    This is the recommended approach when the number of variables exceeds $n\simeq 15$ - $20$. A priori, this heuristic approach is able to explore possible basis interactions of arbitrary order.

All the functions that can be called from `int main()` are declared at the beginning of the `includes/main.cpp` file. The most useful functions are described just above. You can use these functions without modifying the code in `main.cpp` with direct calls from your terminal using the commands described in the following sections.

## Requirements

The code uses the C++11 version of C++.

## Usage without Makefile:

 - **To compile:**
   ```bash
   g++ -std=c++11 -O3 src/*.cpp includes/main.cpp -o BestBasis.out
   ```
 - **To Execute:** The datafile must be placed in the `INPUT` folder.
   
   In the following commands, replace `[datafilename]` by the datafile name, `[n]` by the number of variables, and `[kmax]` by the value of the highest order of operators to consider (when needed).
   
   | Run  | Command | Comment |
   | --- | --- | --- |
   | Help | `./BestBasis.out -h` | |
   | Example 1 | `./BestBasis.out`| Will run Example 1 (Shape dataset)<br> See "Examples" section below|
   | Search best basis | `./BestBasis.out [datafilename] [n]` | Choose automatically<br>the most appropriate algorithm<br> If used, `kmax = 3`by default |
   | Exhaustive search (*) | `./BestBasis.out [datafilename] [n] --exhaustive`| |
   | Search among all operators<br> up to order kmax | `./BestBasis.out [datafilename] [n] --fix-k [kmax]` | specifying kmax is optional,<br> by default `kmax = 3` |
   | Search among all operators<br> up to order kmax<br> in varying representations (**) | `./BestBasis.out [datafilename] [n] --var-k [kmax]` | specifying kmax is optional,<br> by default `kmax = 3` |

   (*) This program implements the exhaustive search algorithm described in Ref.[1].

   (**) This program implements the heuristic algorithm described in Ref.[1].
 
## Usage with Makefile:

Run the following commands in your terminal, from the root folder (which contains the `makefile` document):

 - **To compile:** `make`

 - **To Execute:**
   The datafile must be placed in the `INPUT` folder.
   Information about your dataset must then be specified at the beginning of the makefile, i.e. the name of the datafile `datafile`, the number of variables `n`, and the value of the highest order of operators to consider `kmax` (if needed).

   Open the makefile and replace the values of the following variables at the very top of the file (an example is provided):
    * `datafile`: name of your datafile; this file must be placed in the folder `INPUT`
    * `n`: number of variables in your file; largest possible value is `n = 128`.
    * [if needed] `kmax`: the highest order of operators to consider (if needed); we advise to take it equal to 3 or 4.
   
   You can then execute the code by running in your terminal one of the following commands (from the root folder):

   | Run  | Command | Comment |
   | --- | --- | --- |
   | Help | `make help` | |
   | Example 1 | `make example`| Will run Example 1 (Shape dataset)<br> See "Examples" section below |
   | Search best basis | `make run` | Choose automatically<br>the most appropriate algorithm<br> If used, `kmax = 3`by default |
   | Exhaustive search (*) | `make run-exhaustive`| `kmax` is not used|
   | Search among all operators<br> up to order kmax | `make run-fix-k` | must specify choice of `kmax` |
   | Search among all operators<br> up to order kmax<br> in varying representations (**) | `make run-var-k` | must specify choice of `kmax` |

   (*) This program implements the exhaustive search algorithm described in Ref.[1].

   (**) This program implements the heuristic algorithm described in Ref.[1].

 - **To clean:** `make clean` (to use only once you're done using the code)

## Output files and format of the returned basis:

 - **Terminal Output:** We provided, in the `OUTPUT` folder, the LOGS returned when running the`./BestBasis.out` code on the example datasets.

 - **List of outputs:**
   - **Best Basis found:** will be printed in the terminal, as well as in the file with the name ending in `_BestBasis.dat`;
   - **Inverse of the Best Basis found:** will be printed in the terminal, as well as in the file with the name ending in `_BestBasis_inverse.dat`;
   - **Dataset converted in the new Basis:** will be printed in the file with the name ending in `_inBestBasis.dat`;

 - **Interpreting the printed Basis:**  How to read the output Basis? and where to find it?
   
 - **Results for the examples:** See Ref.[1] for results and discussions on the best basis obtained for these datasets.

## Examples

In the `INPUT` folder, we provided the following examples:
  - **Example 1: Shape data.** The binary dataset `Shapes_n9_Dataset_N1e5.dat` is one of the datasets used as an example in Ref.[1]. This is an artificial dataset generated from only three states and described in Ref.[1] (see Figure 5).
  - **Example 2: Big 5 data.** The binary dataset `Big5PT.sorted`: this is the binarized version of the Big 5 dataset [2] used as an example in Ref. [1]. The dataset has `50` variables and `N = 1 013 558` datapoints. See Ref.[1] for results and comments on the Best Basis obtained for this dataset.
    Important: this dataset is given in a zip file, which must be decompressed first before being used.

    On a laptop, the analysis in varying representations takes about `15min` with `kmax=3` , and about `3h` with `kmax=4`. 

  - **Example 3: MNIST data.** The binary dataset `MNIST11.sorted`: this is the binarized version of the MNIST dataset [3] used as an example in Ref.[1] (see Fig.~7). The dataset has `n=121` variables and `N=60 000` datapoints.
    
    On a laptop, the analysis in varying representations takes about `3h` with `kmax=3` , and about `30h` with `kmax=4`.
    
  - **Example 4: SCOTUS data** the binary dataset `SCOTUS_n9_N895_Data.dat`, which is the dataset of votes of the US Supreme Court analyzed in Ref.[4] and used as an example in Ref.[1].
    
Each of these datasets can be analyzed by running the program with the `makefile` after commenting/uncommenting the appropriate datafile choices at the beginning of the file.

[1]  C. de Mulatier, P. P. Mazza, M. Marsili, *Statistical Inference of Minimally Complex Models*, [arXiv:2008.00520](https://arxiv.org/abs/2008.00520)

[2] Raw data from [Open-Source Psychometrics Project](https://openpsychometrics.org/_rawdata/) in the line indicated as "Answers to the IPIP Big Five Factor Markers"; [here](https://openpsychometrics.org/_rawdata/IPIP-FFM-data-8Nov2018.zip) is a direct link to the same zip file.

[3] LeCun, L Bottou, Y Bengio, P Haffner, *Gradient-based learning applied to document recognition*. Proc. IEEE 86, 2278–2324 (1998).

[4] E.D. Lee, C.P. Broedersz, W. Bialek, Statistical Mechanics of the US Supreme Court. [J Stat Phys 160, 275–301 (2015)](https://link.springer.com/article/10.1007/s10955-015-1253-6).


## License

This code is an open source project under the GNU GPLv3.
