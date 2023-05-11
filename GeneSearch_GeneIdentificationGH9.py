# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:59:59 2022

@author: Fjodor
"""

# Packages
from Bio import SeqIO
import pandas as pd

# D:\4_Uni\MastersThesis\Master's_Thesis (FromDocuments)\Python\LocalBlastTest
# Find all celluloses

# Manual test:
# VmyS280g6323.t1
# TRINITY_GG_21169_c26_g1_i1


"""
Before running:
    
Set working directory
Download files from Vaccinium.org, both Fasta and CSV are required
    
Set filename (without extentions, make sure both the table and the fasta have the same (base) name) # GH9
Set the Workingname extention (e.g. VmPL in the case of PactateLyase) # VmGH9


What has to be done; 
- make sure unaligned contigs are treated properly
- Add protein lengths, PL, Full length (yes or no)

"""

# Xyloglucan_endotransglucosylase
# PactateLyase

Filename = "GH9"

def LocalnBlast (QueryFasta, Database = "TestDBName"):
    

    """
    Script which will run a local nBlast on a specified sequence. 
     
    Command line command for making a database out of the Fasta file:
    makeblastdb -in transcripts.fsa -input_type fasta -parse_seqids -blastdb_version 5 -title "TestDB" -dbtype nucl -out TestDBName
    !! Make sure the working directory is set correctly

    Parameters
    ----------
    QueryFasta : TYPE
        DESCRIPTION.
    Database : TYPE, optional
        DESCRIPTION. The default is "TestDBName".

    Returns
    -------
    Alignment with the best e value

    """
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
    
    blastx_cline = NcbiblastnCommandline(query = QueryFasta, db= Database, evalue = 1e-20, outfmt = 5, out= "my_blast.xml")
    stdout, stderr = blastx_cline()
    
    #Reading the result
    Good_results = []
    for record in NCBIXML.parse(open("my_blast.xml")):
        if record.alignments : #skip queries with no matches
            print("QUERY: %s" % record.query[:60])
            for align in record.alignments:            
                for hsp in align.hsps:
                   # print("MATCH: ", align.title[:60], f"E_value: {hsp.expect}")
                    Good_results.append(align)
    try :
        print(Good_results[0].title)
        return Good_results[0]
    except :
        print("No match was found")
        return "No match was found"
    

FastaFilename = Filename+".fasta"
parsedobj = SeqIO.parse(FastaFilename, "fasta")
parsedlist = []
for item in parsedobj:
    parsedlist.append(item)


Datadictionary = {}
Trinitynames = []
StructuralNames = []

for gene in parsedlist:
    name = gene.name
    namepart = name.split("_")
    lengthofgene = len(gene.seq)
    StructuralNames.append(namepart[0])
    Datadictionary[namepart[0]] = lengthofgene
    with open("fastatemp.fsa", "w") as fasta: 
        fasta.write(f">{gene.name}\n{gene.seq}")
    QueryFasta = "fastatemp.fsa"
    Trinityname = LocalnBlast(QueryFasta)
    try :
        print(Trinityname.title)
        Trinitynames.append(Trinityname.title.split()[0])
    except :
        print(f"No match was found for {gene.name}")
        Trinitynames.append(f"No Trinityname found, {gene.name}")

#Reading the CSV 
TableFilename = Filename+".csv"
DataTable = pd.read_csv(TableFilename)

# Splitting location into seperate Chromosome and Location columns
DataTable2 = DataTable["Location"].str.split(':', expand = True)

# Adding Universal identifier
DataTable2["Universal_identifier"] = StructuralNames

DataTable2.columns = ["Chromosome", "Location", "Universal_identifier"]

DataTable2["Length"] = DataTable2["Universal_identifier"].map(Datadictionary)


DataTable2["TrinityName"] = Trinitynames

#Sorting the table
DataTable2["StartBase"] = DataTable2["Location"].str.split(".").str[0].astype(int)
DataTable2["SortChromosome"] = DataTable2["Chromosome"].str.extract('(\d+)', expand=False).astype(int)
DataTable2 = DataTable2.sort_values(["SortChromosome", "StartBase"])
DataTable2 = DataTable2.drop(["StartBase", "SortChromosome"], axis = 1)

# Working Name:
# Xyloglucan: VmXTH
# PactateLyase: VmPL
WorkingName = "VmGH9-"
WorkingNameSeries = [WorkingName+"%.2d" % (i+1) for i in range (len(parsedlist))]


DataTable2["Study_Name"] = WorkingNameSeries

DataTable2 = DataTable2[["Universal_identifier", "Study_Name", "TrinityName", "Chromosome", "Location", "Length"]]

print(DataTable2)


DataTable2.to_excel("NameTable.xlsx")












