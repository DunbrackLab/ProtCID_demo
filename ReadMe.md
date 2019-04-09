# Protein Common Interfaces Database

## Overview
ProtCID contains comprehensive, PDB-wide structural information on the interactions of proteins and individual protein domains with other molecules, including four types of interactions: chain interfaces, Pfam domain interfaces, Pfam-peptide interfaces and Pfam-ligand/nucleic acids interactions. A common interaction here indicates chain-chain or Pfam domain-domain interfaces that occur in different crystal forms or Pfam-peptide or Pfam-ligand interactions that occur in multiple homologous proteins.


## Repo Contents
-	[**Code**](httpshttps://github.com/DunbrackLab/ProtCID_demo/ProtCid_demo): folder containing source code to use ProtCID dll libraries.
-	[**installer**](httpshttps://github.com/DunbrackLab/ProtCID_demo/ProtCid_demo_setup): a demo program is used to generate interfaces from a user input of PDB crystal structures, cluster the interfaces and output coordinates of each cluster. 

## Web Site
Interfaces and clusters can be searched and downloaded from [ProtCID web site](http://dunbrack2.fccc.edu/protcid)

## Demo
This is a Windows console program to generate chain interfaces on a list of PDB entries, cluster interfaces, output interface files in PDB format, output result text files (e.g. similarity Q scores), compile coordinates of each cluster including PyMol scripts to visualize each cluster. 
This demo program uses ProtCID libraries, and can be installed by the Windows installer.    

## System Requirements
All source code of ProtCID including the demo program is written in C# in Visual studio 2013. It requires Windows operating system and .NET Framework 4.5. The installer will automatically install .NET Framework. 

## Installation Guide
It is very easy to install the demo program. Just download [protcid_demo_setup.msi](https://github.com/DunbrackLab/ProtCID_demo/ProtCid_demo_setup) 
or [setup.exe](https://github.com/DunbrackLab/ProtCID_demo/ProtCid_demo_setup), double click the installer, follow the steps, but to change the installation directory where the program can read and write to it.  

## Instructions for Use

### Synopsis
Protcid_demo –infile ls-pdb.txt –datadir datadir [options]
ProtCid_demo generates interfaces for each PDB entry in “ls-pdb.txt”, clusters interfaces and stores all files in “datadir”. 
```
-infile           a text file containing a list of PDBs, one PDB per line 
```  
```
 -datadir       the path where result data are to be saved
```

### Options
    -alnfile     a text file containing a multiple sequence alignment for the input PDBs. The program will use this file to map residues when calculating similarity Q score of two interfaces from different entries. If this file is not provided, the author residue numbers are used when calculating Q scores. The alignment can be clustal omega format, or a simple text file, one line for each PDB sequence or chain, gaps must be filled by ‘-’, like 
    
```
clustal omega format
CLUSTAL O(1.2.4) multiple sequence alignment
        
1gwnC      MGSSHHHHHHSSGLVPRGSHMDPNQNVKCKIVVVGDSQCGKTALLHVFAKDCFPENYVPT	60
5p21A      -------------------------MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPT	35                                                                                  

1gwnC      VFENYTASFEIDTQRIELSLWDTSGSPYYDNVRPLSYPDSDAVLICFDISRPETLDSVLK	120
5p21A      IEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIH-	94        
```         
OR        
```
1gwnC MGSSHHHHHHSSGLVPRGSHMDPNQNVKCKIVVVGDSQCGKTALLHVFAKDCFPENYVPTVFENYTASFEIDTQRIELSLWDTSGSPYYDNVRPLSYPDS DAVLICFDISRPETLDSVLK
5p21A -------------------------MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTG EGFLCVFAINNTKSFEDIH- 
```

    -groupname       to create a folder so all results are saved into a specific folder with the name, and also named clusters. For instance, -groupname ras. All results are saved into a folder under “datadir”, and clusters are named by ras_cluster ID.tar.gz, e.g. ras_1.tar.gz for the first cluster of user group “ras”. 

### How to run 
Change the current directory to the directory where ProtCID_demo.exe is located. You don't have to change the directory, but provide full path to ProtCID_demo.exe. 
- **1.**	ProtCID_demo –infile demo_data\ls-pdb_ST1A1.txt –datadir demo_data 
This will generate and cluster interfaces from PDB entries in ls-pdb_ST1A1.txt, save all files and coordindates to demo_data folder. The PDB entries contain same protein ST1A1_HUMAN. When calculating Q scores, author residue numbers are used. 
- **2.**	ProtCID_demo –infile demo_data\ls-pdb_ST1A1.txt –datadir demo_data –groupname sulf
This will generate and cluster interfaces from PDB entries in ls-pdb_ST1A1.txt, same all files and coordinates to demo_data\sulf folder and named cluster coordinates files in sulf_cluster ID.tar.gz (e.g. sulf_1.tar.gz).
- **3.**	ProtCID_demo –infile demo_data\ls-pdb_RAS.txt –datadir demo_data –alnfile demo_data\ RasMonomers_clustalO.aln –groupname ras
This will generate and cluster interfaces from PDB entries in ls-pdb_RAS.txt, calculate similarity Q scores by updating all residue numbers to alignment positions (1:N of the alignment). The sequence alignment can be generated by Clustal omega or  other sequence or structure alignment programs like PsiBlast or FATCAT. 

## Results
Output of program include
- **PDB files** (PDB xml files downloaded from rcsb web site, and converted into ProtCID specific xml files. <folder: xml>)
- **Text files**
EntryInfo.txt, EntityChainInfo.txt, EntityUnpDbRef.txt (parsed from PDB xml files)
CrystInterfaces.txt (definition of interfaces computed from crystals)
SameEntryInterfaceCompInfo.txt, DiffEntryInterfaceCompInfo.txt (similarity Q scores between interfaces of each entry, between interfaces of two entries)
InterfaceClusters.txt (Clusters of interfaces)
- **Interface files** (interface coordinate files in PDB format)
- **Cluster coordinates files**  (coordinates of each cluster including Pymol scripts to visualize a cluster )


