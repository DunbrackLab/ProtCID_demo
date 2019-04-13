# Protein Common Interfaces Database

## Overview
ProtCID contains comprehensive, PDB-wide structural information on the interactions of proteins and individual protein domains with other molecules, including four types of interactions: interfaces between full-length protein chains, Pfam domain/domain interfaces, Pfam domain/peptide interfaces and Pfam domain/ligand interactions (including nucleic acids as ligands). A common interactionn occurs in different crystal forms or in non-crystallographic structures across a family of homologous proteins.

Its main goal is to identify and cluster homodimeric and heterodimeric interfaces observed in multiple crystal forms of homologous proteins, and interactions of peptide and ligands in homologous proteins. Such interfaces and interactions, especially of non-identical proteins or protein complexes, have been associated with biologically relevant interactions. For more details about the algorithm and benchmarking, please refer to our paper "Statistical Analysis of Interface Similarity in Crystals of Homologous Proteins." and the [ProtCID web site](http://dunbrack2.fccc.edu/ProtCiD).

We provide source code and dynamically linked libraries for ProtCID. However, we do not suggest rebuilding the ProtCID databases on the entire PDB database since it contains hundreds of GBs of databases and millions of interface files, cluster files, and text files. We provide a web site http://dunbrack2.fccc.edu/protcid, so users can query on our database.

The demo program is used to generate interfaces from a user input of PDB crystal structures, cluster the interfaces and output coordinates of each cluster.

## Repo Contents
-	[**Code**](httpshttps://github.com/DunbrackLab/ProtCID_demo/ProtCid_demo): source code to use ProtCID dll libraries to generate and cluster interfaces, and output results.
-	[**installer**](httpshttps://github.com/DunbrackLab/ProtCID_demo/ProtCid_demo_setup): A Windows installer is to install ProtCid_demo, .NET Framework 4.5 and ProtCID dll libraries. 

## Web Site
http://dunbrack2.fccc.edu/protcid
### Interfaces and clusters can be searched and downloaded from [ProtCID web site](http://dunbrack2.fccc.edu/protcid)
The ProtCID web site contains two services of Windows Communication Foundation (WCF). One is to query the ProtCID database, and the other one is to assign Pfams to input sequences. For details about ProtCID and how to use it, please refer the [HELP pages](http://dunbrack2.fccc.edu/ProtCiD/Help/Help.aspx) on ProtCID web site. 

## Demo
The Demo program is a Windows console program for generating chain interfaces on a list of PDB entries, clustering those interfaces, outputting the interface files in PDB format, outputting result text files (e.g. similarity Q scores), and grouping together coordinate files for interfaces in each cluster and producing PyMol scripts for visualization. 
This demo program uses ProtCID libraries, and can be installed by the Windows installer (ProtCID_demo_setup.msi).    

## System Requirements
All source code of ProtCID including the demo program is written in C# in Visual studio 2013. It requires Windows operating system and .NET Framework 4.5. The installer will automatically install .NET Framework and all required dll libraries. 

## Installation Guide
It is very easy to install the demo program. Just download [protcid_demo_setup.msi](https://github.com/DunbrackLab/ProtCID_demo/tree/master/ProtCid_demo_setup/Release) 
or [setup.exe](https://github.com/DunbrackLab/ProtCID_demo/tree/master/ProtCid_demo_setup/Release), double click the installer, and follow the indicated steps. Please change the installation directory where the program can read and write to it, since ProtCID uses [third-party software](https://github.com/DunbrackLab/ProtCID_demo/tree/master/ProtCid_demo_setup/Release/tools) which might need read and write to the directory.  

## Instructions for Use

### Synopsis
Protcid_demo –infile ls-pdb.txt –datadir datadir [options]
The ProtCid_demo generates interfaces for each PDB entry in “ls-pdb.txt”, clusters interfaces, and stores all files in “datadir”. A cluster contains at least two entries. 
```
-infile        a text file containing a list of PDB IDs (e.g. 5ldg or 5LDG), one PDB per line 
```  
```
-datadir       the path where result data are to be saved
```

### Options
```
-alnfile     a text file containing a multiple sequence alignment for the input PDBs. 
```
The program will use this file to map residues when calculating similarity Q score of two interfaces from different entries. If this file is not provided, the residue numbers in the PDB files (usually named author residue numbers in PDB xml files) are used when calculating Q scores (e.g., this would work for structures of the same protein with identical numbering). The alignment can be [clustal omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) format, or a simple text file, one line for each PDB sequence or chain, gaps must be filled by ‘-’, e.g., 
    
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
```
-groupname   a folder name and filename base for the clusters 
```
For instance, -groupname ras. All results are saved into a folder named "ras" under “datadir”, and clusters are named by ras_cluster ID.tar.gz, e.g. ras_1.tar.gz for the first cluster of user group “ras”. 

### How to run 
Change the current directory to the directory where ProtCID_demo.exe is located. It is not required to change the directory, but please provide the full path to all parameters except -groupname.  
```
ProtCID_demo –infile C:\Users\Qifang\ProtCid_demo_setup\demo_data\ls-pdb_ST2A1.txt –datadir C:\Users\Qifang\ProtCid_demo_setup\demo_data 
```
This will generate and cluster interfaces from PDB entries in ls-pdb_ST2A1.txt, save all files, and coordindates to demo_data folder. The PDB entries contain the same protein, ST2A1_HUMAN. When calculating Q scores, residue numbers in the PDB file (author residue numbers) are used. The cluster 1 in folder "clusterCoordinates" contains the coordindates of the biological dimer of human sulfotransferases. 
```
ProtCID_demo –infile C:\Users\Qifang\ProtCid_demo_setup\demo_data\ls-pdb_ST2A1.txt –datadir C:\Users\Qifang\ProtCid_demo_setup\demo_data –groupname ST2A1
```
This will generate and cluster interfaces from PDB entries in ls-pdb_ST2A1.txt, and save all files and coordinates to demo_data\ST2A1 folder and name cluster coordinates files in ST2A1_cluster ID.tar.gz (e.g. ST2A1_1.tar.gz).
```
ProtCID_demo –infile C:\Users\Qifang\ProtCid_demo_setup\demo_data\ls-pdb_RAS.txt –datadir C:\Users\Qifang\ProtCid_demo_setup\demo_data –alnfile C:\Users\Qifang\ProtCid_demo_setup\demo_data\RasMonomers_clustalO.aln –groupname ras
```
This will generate and cluster interfaces from PDB entries in ls-pdb_RAS.txt, calculate similarity Q scores by updating all residue numbers to alignment positions (1:N of the alignment). The sequence alignment can be generated by [Clustal omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) or  other sequence or structure alignment programs like [PSI-BLAST]( https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) or [FATCAT](http://fatcat.sanfordburnham.org/). 

## Results
Output of the program includes
- **PDB files** (PDB xml files downloaded from [RCSB web site](https://www.rcsb.org/), and converted into ProtCID specific xml files. <folder: xml>)
- **Text files** <folder: textfiles>
   - EntryInfo.txt, EntityChainInfo.txt, EntityUnpDbRef.txt (parsed from PDB xml files)
   - CrystInterfaces.txt (definition of interfaces computed from crystals)
   - SameEntryInterfaceCompInfo.txt, DiffEntryInterfaceCompInfo.txt (similarity Q scores between interfaces of each entry, between interfaces of two entries)
    - InterfaceClusters.txt (Clusters of interfaces)
- **Interface files** (interface coordinate files in PDB format) <folder: CrystInterfaces>
- **Cluster coordinates files**  (coordinates of each cluster including Pymol scripts for visualization <folder: clusterCoordinates>


