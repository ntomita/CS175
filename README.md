# CS175 Bioinfo Final Project

## Installation
### Requirements:

* Python 2.7
* BioPython (http://biopython.org/)

The library can be installed via pip

`pip2 install biopython`

## Usage
Currently, there are two classes,
* Blast
* PMF

Basic usage of these classes is presented in project.py file.
You can simply extend project.py to implement operations for optimization.

## Folders
fasta/ folder contains proteins in fasta format. Place proteins as you need.

xml/ folder is a place to store blasted result, works as a cache. Note: Blast via web api is very slow. Avoid duplicated requests.

etc/ folder, you can put whatever files related to this project.

## TODO:
Implement a class to
* perform optimization.
* evaluate a sequence in 3d structure.
