
HPP = ./cpp/lib
OUT = ./bin
SRC = ./cpp
OMP  = -fopenmp
STD = -std=c++14
fmath = -ffast-math
snake = ./Snakefile

node_size = 0.005
edge_thickness = 0.0005
pdb = ""

ifeq ($(OS), Windows_NT)
	inst = 'powershell "./install.ps1"'
	remove = del /s 
	empty =  nul
	sep = \\
else
	inst = ./install.sh
	remove = rm
	empty = /dev/null
	sep = /
endif

paper_file = main.tex
paper_out = protein_reconstruction
tex_dir = tex

ifeq ($(OS), Cygwin)
    gnu = -std=gnu++14
endif

install: $(inst)
	$(inst) -y

pipe: $(snake)
	snakemake --cores 8

graph_pipe: $(snake)
	snakemake --dag | dot -Tpdf > protein_pipe.pdf

pdb2xyz: $(SRC)/pdb2xyz.cpp
		$(CXX) $(STD) $(fmath) $(gnu) $(OMP) -O3 -o $(OUT)/pdb2xyz $(SRC)/pdb2xyz.cpp

#viewer: $(SRC)/viewer.cpp
#		$(CXX) $(STD) $(fmath) $(gnu) $(OMP) -O3 -o $(OUT)/viewer $(SRC)/viewer.cpp

compare:
	@blender --python ./py/BlenderPDBRec.py -- -T ./protein/$(pdb).pdb.xyz -G ./protein/$(pdb).pdb.guess -s $(node_size) -l $(edge_thickness) 

paper: $(tex_dir)/$(paper_file) \
	   $(wildcard $(tex_dir)/img/**/*)
	cd $(tex_dir) && latexmk -synctex=1 -bibtex -interaction=nonstopmode -file-line-error -pdf $(basename $(paper_file)) -jobname=$(paper_out) && cd ..
	$(MAKE) clean

.PHONY: clean
clean: $(paper_out)
	$(remove) $(tex_dir)$(sep)$(paper_out).blg 2> $(empty)
	$(remove) $(tex_dir)$(sep)$(paper_out).log 2> $(empty)
	$(remove) $(tex_dir)$(sep)$(paper_out).out 2> $(empty)
	$(remove) $(tex_dir)$(sep)$(paper_out).fls 2> $(empty)
	$(remove) $(tex_dir)$(sep)$(paper_out).synctex.gz 2> $(empty)

.PHONY: cleanall
cleanall: $(paper_out) clean
	@$(remove) $(tex_dir)$(sep)$(paper_out).aux 2> $(empty)
	@$(remove) $(tex_dir)$(sep)$(paper_out).bbl 2> $(empty)
	@$(remove) $(tex_dir)$(sep)$(paper_out).fdb_latexmk 2> $(empty)
