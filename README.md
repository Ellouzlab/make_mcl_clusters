
# PURPOSE
Builds genesharing network for multiple purposes, including:
 - phage taxonomic classification
 - understanding plasmidic gene mobility
 - phage gene mobility
 - understanding gene flow between different taxa
 - Building an AMG database



# INSTALL
To install program, navigate into Mobile_clusters and use the command - 

```pip install .```

Please install these dependencies via mamba or conda, otherwise ensure the commands 'mcl', 'mmseqs' and 'diamond' are in your path:
Requires: 
 - diamond >= 2.0.0
   - ```mamba install diamond=2.0.0```
 - mmseqs2 >= 15.6f452
   - ```mamba install mmseqs=15.6f452```
 - mcl >= 22.282
   - ```mamba install mcl=14.137```
  
NOTE: This program may have trouble installing if not your base env.

