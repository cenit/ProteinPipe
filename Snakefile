# snakemake
# snakemake -n # check workflow
# snakemake --dag | dot -Tpdf > protein_pipe.pdf

from __future__ import division
import numpy as np
import pandas as pd

# pipeline of protein reconstruction
import os # local, sep, files, dirs
import sys # exit
import platform # check system
import shutil # move file
from contextlib import suppress # makedir with no error 
import subprocess # capture output
import scipy # squareform, pdist, csgraph, eigsh
import functools # map with multiple args

# plots of results
import seaborn as sns # pretty/easy plots with pandas
import matplotlib.pylab as plt # plots

configfile: "config.yaml"

np.random.seed(int(config["seed"])) # random number generator
sns.set_context("paper", font_scale=1.5)
sep             =   os.sep
local           =   os.path.abspath(".")
if platform.system() == "Windows":
    extension   =   ".exe"
    build       =   "powershell .\\build.ps1"
    shell       =   ".ps1"
elif platform.system() == "Linux" or platform.system()[:6] == "CYGWIN":
    extension   =   ""
    build       =   'bash -c "./build.sh"'
    shell       =   ".sh"

pdb_extension   =   config["pdb_ext"]
atoms           =   list(config["atoms"]) # type of atoms to conserve
thr             =   int(config["thr"]) # threashold distance for contact map
n_population    =   int(config["n_population"]) # number of dna to use in GA
max_iter        =   int(config["max_iter"]) # max number of generation in GA
precision       =   float(config["precision"]) # minimum rmsd accepted
elit_rate       =   float(config["elit"]) # percentage of population to conserve in GA
mutation_rate   =   float(config["mutation"]) # probability of mutation

# folders
pdb_dir         =   config["folders"]["pdb_dir"] # directory of pdb files
protein_dir     =   config["folders"]["protein_dir"] # directory of protein (.xyz) files
guessed_dir     =   config["folders"]["guessed_dir"] # directory of guessed (eigenvectors) files
rec_dir         =   config["folders"]["rec_dir"] # directory of reconstructed protein with GA
cpp             =   config["folders"]["cpp"]
#scripts         =   config["folders"]["python"]

# thread rules
nth_pdb2xyz     =   int(config["NTH_PDB2XYZ"])
nth_guess       =   int(config["NTH_GUESS"])
nth_genetic     =   int(config["NTH_GENETIC"])

ELIT = n_population * elit_rate
HALF = n_population / 2
protein = [f.rsplit('.', 1)[0] for f in os.listdir( os.path.join(local, pdb_dir)) if os.path.isfile(os.path.join( os.path.join(local, pdb_dir), f))]

with(suppress(OSError)):
    os.makedirs(os.path.join(local, protein_dir))
    os.makedirs(os.path.join(local, guessed_dir))
    os.makedirs(os.path.join(local, rec_dir))

def kabsch(pt_true, pt_guessed):
    pt_true -= pt_true.mean(axis = 0)
    pt_guessed -= pt_guessed.mean(axis = 0)
    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD) 
    # of the covariance matrix.
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    V, S, W = np.linalg.svd( np.dot(pt_true.T, pt_guessed) ) # SVD of covariance matrix
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    return np.dot(pt_guessed, np.dot(V, W)) # return pt_guessed rotated by U ( = V * W, rotation matrix)

def random_population(cmap, N, scale = 1e-2):
    weights = [np.triu(np.random.normal(loc=0.0, scale=scale, size=(len(cmap), len(cmap))) * cmap) for i in N]
    return [weights + weights.T for i in N] # symmetric matrix

def get_cmap(protein, thr):
    return np.asarray(scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(protein.iloc[:,1:], metric="euclidean") > thr), dtype=np.int)

def get_lap(cmap):
    return scipy.sparse.csgraph.laplacian(cmap, normed=False)

def laplacian_coords(protein, thr):
    """
    Computing of laplacian (eigenvalues, eigenvectors) of contact map
    """
    vals, vecs = scipy.sparse.linalg.eigsh( get_lap(get_cmap(protein, thr)), k=4)
    return (vals[1:], vecs[:, 1:]) # remove first eigenvalue (null) and corresponding eigenvector

def mutate(cmap, type = "strong", scale=1e-2):
    if type == "strong":
        """
        mutation of an entire row and column
        """
        pos = int(np.random.uniform(low=0, high=len(cmap), size=1))
        new = np.random.normal(loc=0.0, scale = scale, size=len(cmap))
        cmap[pos] = new * cmap[pos]
        cmap[:, pos] = cmap[pos]
    elif type == "week":
        """
        mutation of only one weight
        """
        nn = np.nonzero(A)
        pos = int(np.random.uniform(low=0, high=len(nn[0]), size=1))
        couple = [x[pos] for x in np.nonzero(A)]
        cmap[pos[0], pos[1]] = np.random.normal(loc=0.0, scale = scale, size=1)
    return cmap

def crossover(cmap_a, cmap_b):
    """
    Cross over as submatrix of A extracted as random permutation and B
    """
    rn = np.arange(len(cmap_a))
    np.random.shuffle(rn)
    pos = int(np.random.uniform(low=0, high=len(cmap_a), size=1))
    W = cmap_b
    W[:pos, :pos] = cmap_a[rn[:pos]][:, rn[:pos]]
    return W

def new_generation(old_generation, elit_rate = ELIT, mutation_rate = mutation_rate, half = HALF, type = "strong"):
    # crossover
    if i < ELIT:
        new_gen = population[i]
    else:
        new_gen = crossover(population[idx[np.random.rand() * HALF]], population[idx[np.random.rand() * HALF]])
    # mutation
    if np.random.rand() < mutation_rate:
        new_gen = mutate(new_gen, type)
    return new_gen

def fitness(protein_a, protein_b):
    """
    Sum of square differences of coordinates
    """
    return np.sqrt( np.sum((protein_a - protein_b)**2, axis=1) / len(protein_a))

def protein_pipe(weights, real_coords, guess_cmap, thr=12):
    w = np.sqrt(weights)
    _, guess_coords = laplacian_coords( np.einsum('ij,jk,lk->il', w, guess_cmap, w.T), thr ) # product w*A*w^T
    return fitness(real_coords, kabsch(real_coords, guess_coords) ) 


rule all:
    input:
        expand(os.path.join(local, rec_dir, "{protein}.rec"), protein=protein)

rule build:
    input:
        os.path.join(local, cpp, "pdb2xyz.cpp"),
        os.path.join(local, "build" + shell),
    output:
        os.path.join(local, "bin", "pdb2xyz" + extension)
    benchmark:
        os.path.join("benchmark", "benchmark_build.dat")
    message:
        "Building softwares"
    shell:
        build

rule pdb2xyz:
    input:
        pdbfile  = os.path.join(local, pdb_dir, "{protein}." + pdb_extension),
        software = os.path.join(local, "bin", "pdb2xyz" + extension),
    output:
        out      = os.path.join(local, protein_dir, "{protein}.xyz"),
    benchmark:
        os.path.join("benchmark", "benchmark_pdb2xyz_{protein}.dat")
    threads:
        nth_pdb2xyz 
    message:
        "Conversion PDB2xyz for {wildcards.protein}"
    shell:
        os.path.join(local,
                     ' '.join(["{input.software}",
                               "{input.pdbfile}",
                               ' '.join(atoms)
                               ])
                    )

rule guess_protein:
    input:
        coord_file = os.path.join(local, protein_dir, "{protein}.xyz"),
    output:
        guess_file = os.path.join(local, guessed_dir, "{protein}.guess")
    benchmark:
        os.path.join("benchmark", "benchmark_reconstruction_{protein}.dat")
    threads:
        nth_guess
    message:
        "Simple reconstruction protein for {wildcards.protein}"
    run:
        coords = pd.read_csv(input.coord_file, sep=",", header=None, names=["atoms", "x", "y", "z"])
        vals, vecs = compute_guess(coords, thr)
        pd.DataFrame(vecs, index=vals).to_csv(output.guess_file, sep=",", header=True, index=False)

rule reconstructGA:
    input:
        coord_file = os.path.join(local, protein_dir, "{protein}.xyz"),
        guess_file = os.path.join(local, guessed_dir, "{protein}.guess"),
    output:
        rec_file  = os.path.join(local, rec_dir, "{protein}.rec")
    benchmark:
        os.path.join("benchmark", "benchmark_GA_{protein}.dat")
    threads:
        nth_genetic
    message:
        "GA reconstruction protein for {wildcards.protein}"
    run:
        protein_true  = pd.read_csv(input.coord_file, sep=",", header=None, usecols=["x", "y", "z"], names=["atoms", "x", "y", "z"])
        protein_guess = pd.read_csv(input.guess_file, sep=",", header=None, usecols=["x", "y", "z"], names=["atoms", "x", "y", "z"])
        cmap_guess = get_cmap(protein_guess, thr)
        population = random_population(cmap=get_cmap(protein_guess), N=n_population, scale=1e-2)

        for generation in range(max_iter):
            fit = list(map(functools.partial(protein_pipe, real_coords=protein_true, guess_cmap=cmap_guess, thr=thr), population))
            idx = np.argsort(fit)
            best = idx[0]
            if fit[best] < precision:
                break
            population = list(map(functools.partial(new_generation, elit_rate = ELIT, mutation_rate = mutation_rate, half = HALF, type = "strong"), population))
        population[best].to_csv(output.rec_file, sep=",", header=True, index=False)

