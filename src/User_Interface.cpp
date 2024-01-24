#include <iostream>
#include <sstream>

using namespace std;

/********************************************************************/
/**************************    PARAMETERS    ************************/
/********************************************************************/
/*
// number of binary (spin) variables:
unsigned int n = 9;

// INPUT DATA FILES (optional):  must be in the INPUT folder
std::string input_datafile = "INPUT/Shapes_n9_Dataset_N1e5.dat"; 

unsigned int k_max = 3; // default value
*/

/******************************************************************************/
/**************************  HELP MESSAGE  ************************************/
/******************************************************************************/

void HELP_message()
{
    cout << endl << "*******************************************************************************************";
    cout << endl << "*****************************  HOW TO RUN THE PROGRAM:  ***********************************";
    cout << endl << "*******************************************************************************************" << endl;
    
    cout << "To perform the Basis search:" << endl;
    cout << "\t 1. Place the datafile in the 'INPUT' folder (datafile must be in the correct" << endl; 
    cout << "\t    binary format -- see README file)" << endl;
    cout << "\t 2. Type in your terminal one of the commands in the options below. Remember to:" << endl;
    cout << "\t \t replace [datafilename] by the name of your datafile" << endl;
    cout << "\t \t replace [n] by the number of variables in your data (must be an integer)" << endl;
    cout << "\t \t replace [kmax] by the largest order of operators you want to use" << endl;

    cout << endl << "***************************************  OPTION 0  ****************************************";
    cout << endl << "***********************************  Default example  *************************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "This will run the example:" << endl << endl;
    cout << "\tRun: "; // << endl;
    cout << "   >> ./BestBasis.out" << endl << endl;
    cout << "\tOR: "; //' << endl;
    cout << "    >> make run" << endl;    
 
    cout << endl << "***************************************  OPTION 1  ****************************************";
    cout << endl << "***********************************  Automatic choice  ************************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "This will automatically choose the Search algorithm that is most appropriate for your dataset:" << endl << endl;
    cout << "\tRun: ";
    cout << "   >> ./BestBasis.out [datafilename] [n]" << endl;

    cout << endl << "***************************************  OPTION 2  ****************************************";
    cout << endl << "**********************************  Exhaustive Search  ************************************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    //cout << "Option 2. Choose which Search algorithm to run:" << endl;
    cout << "This will run the Exhaustive Search for the best basis:" << endl;

    cout << "\tRun: ";
    cout << "   >> ./BestBasis.out [datafilename] [n] --exhaustive" << endl << endl;

    cout << "Important: this is not recommended for dataset with more than n~20 variables" << endl;
    cout << "For datasets with more than 25 variables, the program will automatically switch to option 4" << endl << endl;

    cout << endl << "***************************************  OPTION 3  ****************************************";
    cout << endl << "**********************  Search among all operators up to order \'kmax\'  ********************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    //cout << "Option 2. Choose which Search algorithm to run:" << endl;
    cout << "This will run the Search for the best basis *among all operators up to order kmax*" << endl<< endl;

    cout << "\tRun: ";
    cout << "   >> ./BestBasis.out [datafilename] [n] --fix-k" << endl;
    cout << "\tuses the default value of kmax=3" << endl << endl;

    cout << "\tOR: ";
    cout << "    >> ./BestBasis.out [datafilename] [n] --fix-k [kmax]" << endl;
    cout << "\tto specify your choice of \'kmax\'" << endl;    

    cout << endl << "***************************************  OPTION 4  ****************************************";
    cout << endl << "*******************  Search up to order \'k\' in Various Representations  *******************";
    cout << endl << "*******************************************************************************************" << endl << endl;

    cout << "The program will run the Search for the best basis *among all operators up to order kmax* in varying successive representations" << endl;
    cout << "See README for a description of the algorithm." << endl<< endl;

    cout << "\tRun: ";
    cout << "   >> ./BestBasis.out [datafilename] [n] --var-k" << endl;
    cout << "\tuses the default value of kmax=3" << endl << endl;

    cout << "\tOR: ";
    cout << "    >> ./BestBasis.out [datafilename] [n] --var-k [kmax]" << endl;
    cout << "\tto specify your choice of \'kmax\'" << endl; 

    cout << endl << "*******************************************************************************************" << endl;
    cout << endl;
}


/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/

int Read_argument(int argc, char *argv[], string *input_datafile, unsigned int *n, unsigned int *k_max)
{
    int flag_search = 1;   // 1 by default for the example
    // 1 = Exhaustive search
    // 2 = Fixed basis search with given choice of k_max
    // 3 = Varying basis search with given choice of k_max

	/**********************     READ ARGUMENTS    *********************************/

    // argv[0] contains the name of the datafile, from the current folder (i.e. from the folder containing "data.h");
    // argv[1] contains the number of variables to read;
    if (argc == 2)
    {
        string help = argv[1];
        if (help == "-h")
            { cout << endl << "HELP:" << endl;}
        else
            {cout << endl << "ERROR: The number of arguments is not correct." << endl;}
        HELP_message();
        return 0;        
    }
    else if (argc > 5 || argc == 2)
    {
        cout << endl << "ERROR: The number of arguments is not correct." << endl;
        HELP_message();
        return 0;
    }
    else if (argc >= 3)
    {
        string input_datafile_buffer = argv[1];
        string n_string_buffer = argv[2];

        (*input_datafile) = "INPUT/" + input_datafile_buffer;
        (*n) = stoul(n_string_buffer);

        if (argc == 3)
        {
            if((*n)<15) // Exhaustive search
                { flag_search = 1;  } 
            else     // Varying basis search with default k_max
                { flag_search = 3;  } 
        }
        else //if (argc == 4 || argc == 5)
        {
            string flag = argv[3];
            if (flag == "--exhaustive") // Exhaustive search
            {
                if ((*n) < 25) { flag_search = 1; }
                else {  // Varying basis search with default k_max
                    cout << "The Exhaustive Search is not recommended for dataset with more than n~20-25 variables." << endl;
                    cout << "We Recommend to use the heuristic procedure instead." << endl;
                    flag_search = 3;
                }
            }
            else if (flag == "--fix-k") // Fixed basis search with default k_max
            {
                flag_search = 2;
                if (argc == 5)
                    { (*k_max) = stoul(argv[4]); }
            }
            else if (flag == "--var-k") // Varying basis search with default k_max
            {
                flag_search = 3;
                if (argc == 5)
                    { (*k_max) = stoul(argv[4]); }
            }
            else 
            {
                cout << "ERROR: The arguments are not correct." << endl << endl;
                HELP_message();
                return 0; 
            }
        }
    }

    return flag_search;
}




