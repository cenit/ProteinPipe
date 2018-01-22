#required OPENMP 4.0 (since gcc 4.9)
HPP = ./cpp/lib
OUT = ./bin
CPP = ./cpp
STD = -std=c++14
fmath = -ffast-math
#fmath = -ffp-mode=fast
OMP  = -fopenmp

ifeq ($(shell uname -o),Cygwin)
    gnu = -std=gnu++14
endif

ifeq ($(OS), Windows_NT)
	MKDIR_P = mkdir $(subst /,\,$(OUT)) > nul 2>&1 || (exit 0)
	INSTALL = ./install.ps1
else
	MKDIR_P = mkdir -p $(OUT) 
	INSTALL = ./install.sh
endif

all: 	install \
		protein_pipe \

install: 	install.sh \
		 	install.ps1
		$(INSTALL) -y

protein_pipe:	Snakefile

		snakemake

pdb2xyz:	$(CPP)/pdb2xyz.cpp \

		$(CXX) $(STD) $(fmath) $(gnu) $(OMP) -O3 -o $(OUT)/pdb2xyz $(TST)/pdb2xyz.cpp

dir_tree:
		$(MKDIR_P)
