########################################################################################################################
####################################      CONSTANT TO SPECIFY     ######################################################
########################################################################################################################
#### If the following two variables are left empty:           ##########################################################
####     then the program will automatically use the filename and number of variables specifed in 'main.cpp'   #########
########################################################################################################################

#### EXAMPLE: Big 5 dataset:
#datafilename := Big5-IPC1_VS3_Ne5.dat  # datafile name ### IMPORTANT: this file must be in the 'INPUT' folder
#n := 50		# number of binary variables in the datafile 
k := 4

#### EXAMPLE 2: Shapes:
datafilename := Shapes_n9_Dataset_N1e5.dat  
n := 9		


########################################################################################################################
######## ENTER THE FOLLOWING IN YOUR TERMINAL:
#### TO COMPILE:  	make
#### TO RUN: 		make run
#### TO CLEAN:  	make clean    --> to use only when you are completely done
########################################################################################################################

### g++ -std=c++11 -O3 src/*.cpp includes/main.cpp -o BestBasis.out

########################################################################################################################
##########################################      DO NOT MODIFY     ######################################################
########################################################################################################################
CC = g++ 	# Flag for implicit rules: used for linker
CXX = g++ 	# Flag for implicit rules: compilation of c++ files
CXXFLAGS = -std=c++11 -O3  #-Wall  #Extra flags to give to the C++ compiler

### Directory for Files:
DIR_Basis = src

### Files:
objects = tools.o ReadDataFile.o User_Interface.o
OBJS := $(objects:%=$(DIR_Basis)/%)

objectsWdataH = Init_OpSet.o ExtractBasis_inOpSet.o BasisTools.o BestBasis_IterativeSearch.o BestBasis_ExhaustiveSearch.o
OBJS_Wdata := $(objectsWdataH:%=$(DIR_Basis)/%)

### Compilation -- Implicite rule:
#%.o : %.c   
#		$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

### Link -- Implicite rule:
#$(CC) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

BestBasis.out: $(OBJS) $(OBJS_Wdata) includes/main.o 
	g++ $(CXXFLAGS) includes/main.o $(OBJS) $(OBJS_Wdata) -o BestBasis.out

main.o: main.cpp src/data.h
	g++ $(CXXFLAGS) -c includes/main.cpp -o includes/main.o

Init_OpSet.o: Init_OpSet.cpp src/data.h
	g++ $(CXXFLAGS) -c Init_OpSet.cpp -o Init_OpSet.o

ExtractBasis_inOpSet.o: ExtractBasis_inOpSet.cpp src/data.h
	g++ $(CXXFLAGS) -c ExtractBasis_inOpSet.cpp -o ExtractBasis_inOpSet.o 

BasisTools.o: BasisTools.cpp src/data.h
	g++ $(CXXFLAGS) -c BasisTools.cpp -o BasisTools.o

BestBasis_IterativeSearch.o: BestBasis_IterativeSearch.cpp src/data.h
	g++ $(CXXFLAGS) -c BestBasis_IterativeSearch.cpp -o BestBasis_IterativeSearch.o

BestBasis_ExhaustiveSearch.o: BestBasis_ExhaustiveSearch.cpp src/data.h
	g++ $(CXXFLAGS) -c BestBasis_ExhaustiveSearch.cpp -o BestBasis_ExhaustiveSearch.o

########################################################################################################################
####################################################      RUN     ######################################################
########################################################################################################################

example:
	time ./BestBasis.out

help:
	./BestBasis.out -h

run:
	./BestBasis.out $(datafilename) $n

run-fix-k:
	time ./BestBasis.out $(datafilename) $n --fix-k $k

run-var-k:
	time ./BestBasis.out $(datafilename) $n --var-k $k

########################################################################################################################
##################################################      CLEAN     ######################################################
########################################################################################################################

clean:
	rm -f includes/main.o $(OBJS) $(OBJS_Wdata) BestBasis.out

