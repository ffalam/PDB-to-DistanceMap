
'''
@Author : Fardina F Alam


@Task:
-> Convert PDB Co-ordinates to Distance Map

@Basic Concept need to Know:

We can represent a tertiary (3D) protein structure as a contact map or distance matrix, which are both 2D-based representations of a 3D structure.

@Code Description: 

-> This code will download PDB files from RCSB website automatically and then will read these PDBs and using PDBParser, 
will generate 3D ( tetiary) protein structures.

-> In the next step, it will calculate distanc map/ matrics for each PDB so that, each 3D protein structure could be represent as 2D distance map. 
We can crop each distance map to any size (example: 64 by 64, 20 by 20 etc.) if needed.  


@Requirement:

Install BioPython: pip install biopython

'''

import os 
import scipy
import random
import Bio.PDB
import numpy as np
from pathlib import Path
from scipy import spatial
import urllib.request
from Bio.PDB.PDBParser import PDBParser
import sys



# Here, I have defined 3 functions for calculating distance map from given coordinates of PDB. We can use any one of them

def calc_distance_map(coordinates) :
    """
    Calculate distance map

    :param : coordinates of a PDB file
    :return: equivalent distance map
    """

    distance_mat = np.sum((coordlist[:,np.newaxis,:] - coordlist[np.newaxis,:,:])**2, axis=-1)**(1./2)
    #print (distance_mat.shape)
    return distance_mat


# Alternative way to calculate distance map using norm
def calc_distance_map_using_norm(coordinates) :
    """
    Calculate distance map

    :param : coordinates of a PDB file
    :return: equivalent distance map
    """
    map_arr=[]
    for i in range (len(coordlist)):
        res_arr=[]
        a=[]
        for j in range(len(coordlist)):
            #print(i,j)
            res_arr.append(np.linalg.norm(coordlist[i]-coordlist[j])) # Calculating Distance Map using np.linalg.norm
        map_arr.append(res_arr)
    return map_arr



# Alternative way to calculate distance map using norm
def calc_distance_map_using_cdist(coordinates) :
      
    """
    Calculate distance map

    :param : coordinates of a PDB file
    :return: equivalent distance map
    """  
    map_arr=scipy.spatial.distance.cdist(coordinates,coordinates,'Euclidean')
    #print(map_arr.shape)
    return map_arr

 

# Here, I have defined the function to download PDB files from RCSB website automatically

def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    count = 0
    pdbfn = pdbcode + ".pdb"
    #print("pdbfn",pdbfn)
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        count +=1
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None


#--------------------------------------------Download PDB files------------------------------------------------------------------------------------#


datadir = "Spring2022/"
file_count = 0

print("Downloading PDB files")
with open('pdb-list.txt') as inputfile:
    for line in inputfile:
         j=line.rstrip('\n') 
         pdb_id=j
         print("pdb_id",pdb_id)
         download_pdb(pdb_id,datadir)



#--------------------------------------------Extracting Co-ordinates from PDB------------------------------------------------------------------------------------#

print("Extracting Co-ordinates from PDB")
with open('pdb-list.txt') as inputfile:
    for line in inputfile:                                             # Read each PDB name from "PDB-List"
        print("-------------------------------------------")
        pdb_id=line.rstrip('\n')
        print("pdb_id",pdb_id)
        length=len(pdb_id)
        structure_id = pdb_id

        p1 = PDBParser(PERMISSIVE=1)                                   # PDBPARSER parses every PDB file 
        pdbfn = structure_id + ".pdb"
        
        if os.path.isfile(datadir+pdbfn):
            print("File exists") 
            file_count += 1    
            structure= p1.get_structure(structure_id,datadir+pdbfn)    # Used PDBParser to extract pdb structure directly from url: Example : # p1.get_structure("3NIR", '3NIR.pdb') 
            structure1 = structure[0]                                  # Considering only first structure of the PDB as there may be more than one structures: User choice     

            # Extracting C-Alpha (CA) co-ordinates only from PDB: user choice

            chain_id='A'                                               # Considering only Chain A : user choice
            coordlist=[]
            for chain in structure1:
                                print("chain",chain.get_id())
                                if(chain.get_id()==chain_id):
                                    for residue in chain:
                                        for atom in residue:
                                            if(atom.get_id()=="CA"):
                                                coordlist.append(atom.get_coord())

            np.savetxt(datadir+pdb_id+"_coord.txt",coordlist)          # Saving CA Coordinates of PDBS 
            


#--------------------------------------------Generate Distance Maps------------------------------------------------------------------------------------#

            distance_map_arr=np.array(calc_distance_map_using_cdist(coordlist))
            #distance_map_crop= distance_map_arr[:20,:20]                  # If needed, we can crop distance map to any size : user choice
            #print("distance_map_crop",distance_map_crop.shape)
            
            print("STRUCTURE with Distance map shape",structure_id, distance_map_arr.shape)
            np.savetxt(datadir+pdb_id+"_distancemap.txt",distance_map_arr). # Saving Distance Map
            print("-------------------------------------------")
        else:
            print("PDB File does not exist! IOError has occured")
            continue


print("Total PDB Structure Saved:", file_count)