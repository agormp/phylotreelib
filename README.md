# Treelib: python module for manipulating phylogenetic trees

Using classes and methods in treelib.py it is possible to read treefiles in either
NEXUS or Newick format (each file can contain one or more trees), and to analyze 
and manipulate the trees in various ways. 

The main idea is to construct a treefile object from text (e.g., obtained by 
reading a NEXUS or Newick treefile). Treefile objects contain tree objects, and 
can be iterated over. Tree objects have a number of methods that can be used to 
analyze or alter the tree in question.

In the descriptions below it is assumed that treelib has been imported 
("import treelib"). The treelib.py library file can be placed anywhere in the user's 
PYTHONPATH, or simply in the same directory as the calling script.

