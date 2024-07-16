
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
 - mmseqs2 >= 15.6
   - ```mamba install mmseqs2=15.6```
 - mcl >= 22.282
   - ```mamba install mcl=22.282```
  
NOTE: This program may have trouble installing if you are not in your base env. This will be fixed in the future.

