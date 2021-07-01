# Treelib: python module for manipulating and analyzing phylogenetic trees

Using classes and methods in treelib.py it is possible to read treefiles in either
NEXUS or Newick format (each file can contain one or more trees), and to analyze
and manipulate the trees in various ways.

The main idea is to construct a treefile object from text (e.g., obtained by
reading a NEXUS or Newick treefile). Treefile objects contain tree objects, and
can be iterated over. Tree objects have a number of methods that can be used to
analyze or alter the tree in question.

The treelib.py library file can be placed anywhere in the user's
PYTHONPATH, or simply in the same directory as the calling script.

## Constructing treefile objects

Typically, treefile objects are constructed by reading a textfile containing a
NEXUS or Newick format description of a tree. This is done as follows:

Newick files:
```python
  treefile = treelib.Newicktreefile(filename)
```

NEXUS files:
```python
  treefile = treelib.Nexustreefile(filename)
```

## Constructing tree objects

Treefile objects contain tree objects. Tree objects are retrieved from treefile
objects by iteration:

Doing something to all trees in a treefile object:
```python
  for tree in treefile:
      <do something with tree>
```

Getting just one tree from a treefile object:
  tree = next(treefile)

A tree object can also be constructed directly from a string (where the string is a Newick formatted tree):
  tree = treelib.Tree.from_string(mystring)

Tree objects consist of external nodes (leafs), which are identified by strings
(e.g. "Chimpanzee"), and internal nodes, which are identified by integers
(e.g., 5). Branches between nodes may have a label (string) and/or a branch
length (float) associated with them.

One functionality in treelib that is handy during debugging is that tree objects
can be printed in fairly clear format using simply: "print(tree)". The resulting
output is a child-list representation of the tree followed by an alphabetical
list of leafs, along these lines:

>>> print(tree)
|------------------------------------------------------------------|
|  Node  |              Child               |  Distance  |  Label  |
|------------------------------------------------------------------|
|     0  |  tetM_AF333235_Clostridium_diff  |  0.045404  |         |
|     0  |                         9830414  |  0.003285  |         |
|     0  |                         9830359  |  0.001099  |         |
|     0  |                               1  |  0.018480  |   1.00  |
|     1  |                              24  |  0.007180  |   1.00  |
|     1  |                               2  |  0.026523  |   1.00  |
|    24  |                         9830470  |  0.001096  |         |
|    24  |                         9830409  |  0.000563  |         |
|     2  |                               3  |  0.016300  |   0.77  |
|     2  |  tetO_AY190525_Campylobacter_je  |  0.565353  |         |
|     3  |                              17  |  0.007065  |   0.82  |
|     3  |                               4  |  0.006912  |   1.00  |
|    17  |                              18  |  0.004701  |   1.00  |
|    17  |                              19  |  0.001937  |   0.87  |
|     4  |                               5  |  0.003834  |   1.00  |
|     4  |                               6  |  0.007932  |   1.00  |
|    18  |  tetM_X90939_Streptococcus_pneu  |  0.000804  |         |
|    18  |   tetM_AF376746rc_Streptococcus  |  0.001385  |         |
|    19  |                              20  |  0.001780  |   0.98  |
|    19  |                              21  |  0.007687  |   0.97  |
|     5  |                           20028  |  0.001084  |         |
|     5  |                           18854  |  0.000562  |         |
|     6  |  tetM_U58986_Gardnerella_vagina  |  0.007218  |         |
|     6  |                               7  |  0.010529  |   1.00  |
|    20  |  tetM_AJ580978_Streptococcus_mi  |  0.004311  |         |
|    20  |  tetM_AJ580977_Streptococcus_mi  |  0.002311  |         |
|    21  |  tetM_L12242_Neisseria_gonorrho  |  0.001760  |         |
|    21  |                              22  |  0.011284  |   1.00  |
|     7  |  tetM_L12241_Neisseria_gonorrho  |  0.001992  |         |
|     7  |                               8  |  0.008244  |   1.00  |
|     8  |                               9  |  0.006822  |   1.00  |
|     8  |  tetM_U58985_Gardnerella_vagina  |  0.000684  |         |
|    22  |  tetM_X04388_Enterococcus_faeca  |  0.000610  |         |
|    22  |                              23  |  0.002168  |   1.00  |
|     9  |                              10  |  0.002313  |   1.00  |
|     9  |                              12  |  0.001110  |   0.83  |
|    23  |                           20074  |  0.000577  |         |
|    23  |  tetM_AB054984_Clostridium_sept  |  0.002784  |         |
|    10  |                           20032  |  0.000537  |         |
|    10  |                              11  |  0.002783  |   1.00  |
|    10  |  tetM_X92947_Enterococcus_faeca  |  0.001092  |         |
|    10  |                         9830090  |  0.000543  |         |
|    10  |                         9830498  |  0.000536  |         |
|    10  |  tetM_M85225_Enterococcus_faeca  |  0.000536  |         |
|    10  |  tetM_U09422_Enterococcus_faeca  |  0.000541  |         |
|    10  |                         9830479  |  0.000547  |         |
|    10  |                           18836  |  0.000542  |         |
|    10  |                           15109  |  0.000540  |         |
|    10  |                         9830457  |  0.000536  |         |
|    10  |                         9830491  |  0.000537  |         |
|    12  |                              16  |  0.001123  |   0.99  |
|    12  |                              13  |  0.005133  |   1.00  |
|    12  |                              14  |  0.002217  |   1.00  |
|    16  |  tetM_ABO39845_Erysipelothrix_r  |  0.000567  |         |
|    16  |  tetM_AF329848_Clostridium_perf  |  0.000547  |         |
|    11  |  tetM_X56353_Enterococcus_faeca  |  0.000559  |         |
|    11  |  tetM_AJ580976_Streptococcus_or  |  0.004525  |         |
|    13  |  tetM_U08812_Ureaplasma_urealyt  |  0.004481  |         |
|    13  |  tetM_X75073_Neisseria_meningit  |  0.001664  |         |
|    13  |  tetM_AOE14233_Streptococcus_ag  |  0.000558  |         |
|    13  |  tetM_AF440277_Lactobacillus_pl  |  0.001113  |         |
|    14  |   tetM_AF491293_Bacillus_cereus  |  0.000589  |         |
|    14  |                              15  |  0.022849  |   1.00  |
|    15  |  tetM_M21136_Staphylococcus_aur  |  0.001704  |         |
|    15  |  tetM_APO03359rc_Staphylococcus  |  0.000573  |         |
|------------------------------------------------------------------|

41 Leafs:
------------------------------
15109
18836
18854
20028
20032
20074
9830090
9830359
9830409
9830414
9830457
9830470
9830479
9830491
9830498
tetM_AB054984_Clostridium_sept
tetM_ABO39845_Erysipelothrix_r
tetM_AF329848_Clostridium_perf
tetM_AF333235_Clostridium_diff
tetM_AF376746rc_Streptococcus
tetM_AF440277_Lactobacillus_pl
tetM_AF491293_Bacillus_cereus
tetM_AJ580976_Streptococcus_or
tetM_AJ580977_Streptococcus_mi
tetM_AJ580978_Streptococcus_mi
tetM_AOE14233_Streptococcus_ag
tetM_APO03359rc_Staphylococcus
tetM_L12241_Neisseria_gonorrho
tetM_L12242_Neisseria_gonorrho
tetM_M21136_Staphylococcus_aur
tetM_M85225_Enterococcus_faeca
tetM_U08812_Ureaplasma_urealyt
tetM_U09422_Enterococcus_faeca
tetM_U58985_Gardnerella_vagina
tetM_U58986_Gardnerella_vagina
tetM_X04388_Enterococcus_faeca
tetM_X56353_Enterococcus_faeca
tetM_X75073_Neisseria_meningit
tetM_X90939_Streptococcus_pneu
tetM_X92947_Enterococcus_faeca
tetO_AY190525_Campylobacter_je
