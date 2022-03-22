# PDB-to-DistanceMap

## **Task**
Convert PDB Co-ordinates to Distance Map

###### **Basic Concept Needs to Know:**
We can represent a tertiary (3D) protein structure as a contact map or distance matrix, which are both 2D-based representations of a 3D structure.

###### **Code Description:**

-> This code will download PDB files from RCSB website automatically and then will read these PDBs and using PDBParser, 
will generate 3D ( tertiary) protein structures.

-> In the next step, it will calculate distanc map/ matrix for each PDB so that, each 3D protein structure could be represent as 2D distance map. 
We can crop each distance map to any size (example: 64 by 64, 20 by 20 etc.) if needed.  


###### **Requirement:**

Install BioPython: pip install biopython

