#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 21:18:50 2020

@author: nicolebowers
Input: VCF
Output: CSV
"""
#import packages
import pandas as pd 

#import package to help convert the VCF file into a dataframe https://github.com/dceoy/pdbio
from pdbio.vcfdataframe import VcfDataFrame

vcf_path = '/Users/nicolebowers/Challenge_data.vcf'

vcfdf = VcfDataFrame(path=vcf_path)

#name variable now that it is a dataframe 
df=vcfdf.df

#Setting up new columns for each item of information that will be found within the from INFO column
#The numbers correlate to the Tempus Bioinformatics Tehcnical Challenge guidelines 

#1. to annotate type of variation:
#using TYPE from INFO column 

df["TYPE"]= ""

#2. to annotate depth of sequence coverage at the site of variation
#using DP Read Depth 
#Description:"Total read depth at the locus"

df["DP"]= 0.0

#3. to annotate number of reads supporting variant 
#Using AO
#Description: Alternate allele observation count

df["AO"]= 0.0

#4. To calculate the precentage of reads supporting the variant vs. those supporting the reference reads 
#Using AO
#Alternate allele observation count 

#Using RO
#Alternate Reference allele observation count

df["RO"]= 0.0
df["read_percentage"] = ""
df["SNP Type"]=""

#This is a loop to extract information from the paired keys and values from the INFO column 
#Then puts the outputs into their own columns 
for item in range(len(df)):
    info = df["INFO"][item].split(";") #; separate the keys and values 
   #Slicing strings to remove the keys, keeping values
    df["AO"][item]= info[5][3:] 
    df["DP"][item]= info[7][3:]
    df["RO"][item]= info[28][3:]
    df["TYPE"][item]= info[40][5:]

#This loop handles the Alt Read count for reads with insetions which contain a comma (canot be int or float type)
for item in range(len(df)):
        if type(df["AO"][item]) == str:
            ao = df["AO"][item].split(",")
            ao = [int(a) for a in ao]
            ro = df["RO"][item]
            o = []
#Find the percentage of alt read by dividing the number of alt reads by total reads (alt reads + reference reads)      
            for a in ao:
                o.append(a/(a+ro))
  
#Make all items in the read_percentage column into the same type: list (rather than a mix of floats and list items )
            
            df["read_percentage"][item] = o
        else: 
            df["read_percentage"][item]= [df["AO"][item]/(df["AO"][item] + df["RO"][item])]

#identify the type of SNP present 
#Possible Transversions
for item in range(len(df)):
    if df["REF"][item] == "A" and df["ALT"][item] == "C" or \
        df["REF"][item] == "C" and df["ALT"][item] == "A" or \
        df["REF"][item] == "G" and df["ALT"][item] == "T" or \
        df["REF"][item] == "T" and df["ALT"][item] == "G" or \
        df["REF"][item] == "T" and df["ALT"][item] == "A" or \
        df["REF"][item] == "A" and df["ALT"][item] == "T" or \
        df["REF"][item] == "C" and df["ALT"][item] == "G" or \
        df["REF"][item] == "G" and df["ALT"][item] == "C":
        df["SNP Type"][item] = "Transversion"
#Posible Transitions 
    elif df["REF"][item] == "C" and df["ALT"][item] == "T" or \
        df["REF"][item] == "T" and df["ALT"][item] == "C" or \
        df["REF"][item] == "G" and df["ALT"][item] == "A" or \
        df["REF"][item] == "A" and df["ALT"][item] == "G":
        df["SNP Type"][item] = "Transition"



'''
For item:
5. Allele frequency of variant from ExAC API (API documentation is available here:
http://exac.hms.harvard.edu/).

I could not figure out how to read the VCF file into the ExAC API tool. So instead, I annotated the 
ExAC frequencies directly into the file using the online version of ANNOVAR (http://wannovar.usc.edu/)
This tool allowed me to upoad the VCF directly. I merged the ExAC frequences into the output dataframe.
'''
#read in output from ANNOVAR
ANNOVAR = pd.read_csv("/Users/nicolebowers/query.output.genome_summary.csv")

#Merge two Dataframes on index of both the dataframes
mergedDf = df.merge(ANNOVAR["ExAC_Freq"], left_index=True, right_index=True)

#Export csv        
mergedDf.to_csv(r'/Users/nicolebowers/Annotated_Tempus_Challenge.csv', index = False)
'''
See first tab of Annotated_Tempus_Challenge.csv for instructions as to how to read this document.
See tab two for the output data.

Thank you Tempus Bioinformatics Team for reviewing my code.
Nicole 
'''

