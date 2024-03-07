#! /usr/bin/env python

import argparse
import os,sys
from decimal import Decimal
from substrate import get_substrate_data
from chebi_parser import ce_categorization,role_categorization


'''
A program developed to aid in performing cross-transcriptome studies using the Transporter
Classfication Database (TCDB). This program takes in multiple analysis files, previously
processed by evaluateAlignments.py. It then creates a combined table that contains TCDB
proteins and the most similar proteins (determined through the analysis) from each of
the proteomes.
'''

def readLines(file):
    """
    A function to open a file and read the lines

    Parameters
    ----------
    file
        the path to the file to be read

    Returns
    -------
    lines
        the list containing the file line by line
    """
    lines = open(file,'r').readlines()

    if lines[0][0] == '#':
        return lines[1:]

    return lines


def mapTCID(lines,tcids,accessions,tms,genome,repeats):
    """
    For each analysis file, updates the dictionaries that contain TCIDs, Accessions, TMS data,
    proteome specifc data, and repeat protein hits.
    """
    tcidMap = {}

    for line in lines:

        fields = line.rstrip().split('\t')

        acc = fields[2]
        tms[acc] = fields[3]
        tcid = fields[4]
        query = [fields[0]]+fields[5:]
        tms[query[0]] = fields[1]

        tcids.add(tcid)

        #Overall Data
        if tcid not in accessions:
            accessions[tcid] = []

        if acc not in accessions[tcid]:
            accessions[tcid].append(acc)

        #Genome-specifc data
        if tcid not in tcidMap:
            tcidMap[tcid] = {}

        if acc in tcidMap[tcid]:

            if '{}:{}-{}'.format(genome,tcid,acc) not in repeats:
                repeats.append('{}:{}-{}'.format(genome,tcid,acc))

            if float(query[1]) < float(tcidMap[tcid][acc][1]):
                tcidMap[tcid][acc] = query

        else:

            tcidMap[tcid][acc] = query

    return tcidMap,tcids,accessions,tms,repeats


def getTCID(line):
    """
    Helper function to grab TCID from an input file
    """
    fields = line.split('\t')

    return fields[2]


def getSubstrate(tcid,substrate_data,primary_ce,secondary_ce,primary_role,secondary_role):
    """
    Given a specifc TCID, search for the substrate and attempt to classify it using the provided
    primary and secondary chemical categorizations and role categorizations
    """
    ce = []
    role = []

    if tcid in substrate_data:

        for substrate in substrate_data[tcid]:

            id,name = substrate

            cat = ce_categorization(id,primary_ce=primary_ce,secondary_ce=secondary_ce)

            if cat != (None,None):
                ce.append('{}({})-{}({})'.format(cat[0],cat[1],name,id))
            else:
                ce.append('{}({})'.format(name,id))

            r_cat = role_categorization(id,primary_role=primary_role,secondary_role=secondary_role)

            if r_cat != (None,None):
                role.append('{}({})-{}({})'.format(r_cat[0],r_cat[1],name,id))
            else:
                role.append('{}({})'.format(name,id))

        return ', '.join(ce),', '.join(role)

    return 'none','none'


def getProtein(line):
    """
    Helper Function to get Protein Accession
    """
    fields = line.split('\t')

    return fields[0]


def addEntry(tcidMap,tcid,protein):
    """
    Adds a TCID and the prt
    """
    if tcid not in tcidMap:

        tcidMap[tcid] = []

    tcidMap[tcid].append(protein)

    return tcidMap


def getGenomes(directory):
    """
    helper Function to get genome name and file paths
    """
    genomes = []
    genomeFiles = {}

    for file in os.listdir(directory):

        filePath = '{}/{}'.format(directory,file)
        fileName = file.replace('.tsv','')

        genomes.append(fileName)

        genomeFiles[fileName] = filePath

    return genomes,genomeFiles


def printTable(genomes,tcids,tcidMaps,accessions,tms,substrate_data,primary_ce,secondary_ce,primary_role,
               secondary_role,repeats,outdir,output):
    """
    Uses all the previous calculations and outputs the final table.
    """
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    outputFile = open('{}/{}'.format(outdir,output),'w')
    repeatFile = open('{}/repeats.tsv'.format(outdir),'w')

    queryFields = '\t'.join(['query\tq_tms\tevalue\tpident\tqcov\tscov']*len(genomes))

    outputFile.write('#TCID\tAcc\tCE\tRole\thit_tms_no\t{}\t{}\n'.format('\t'.join(genomes),queryFields))

    for tcid in tcids:

        for acc in accessions[tcid]:

            hits = []
            pos = []

            for genome in genomes:

                if tcid in tcidMaps[genome]:

                    if acc in tcidMaps[genome][tcid]:

                        query = tcidMaps[genome][tcid][acc][0]
                        tms_no = tms[query]
                        values = tcidMaps[genome][tcid][acc][1:]

                        hits.append('\t'.join([query]+[tms_no]+values))
                        pos.append('+')

                    else:

                        hits.append('none\tnone\tnone\tnone\tnone\tnone')
                        pos.append('-')
                else:

                    hits.append('none\tnone\tnone\tnone\tnone\tnone')
                    pos.append('-')

            ce,role = getSubstrate(tcid,substrate_data,primary_ce,secondary_ce,primary_role,secondary_role)

            outputFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(tcid,acc,ce,role,tms[acc],
                                                                   '\t'.join(pos),'\t'.join(hits)))

    outputFile.close()

    for repeat in repeats:
        repeatFile.write('{}\n'.format(repeat))

    repeatFile.close()


def parse_arguments():
    """
    Argument Parser for CLI
    """
    parser = argparse.ArgumentParser(description="A command line tool to process proteome analysis tables")

    parser.add_argument('-i', '--input_directory', action='store',
                        help='The path to the directory containing all the tsv files containing analysis.')
    parser.add_argument('-of', '--outfile', action='store',
                        help='The name of the final file containing the master table.')
    parser.add_argument('-od', '--outdir', action='store',
                        help='The path to the directory where all the analysis files will be output.')

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    return args.input_directory, args.outdir, args.outfile

if __name__ == "__main__":


    #Initialize tcids
    tcids = set()

    #initialize maps
    tcidMaps = {}
    substrates = {}
    accessions = {}
    tms = {}
    repeats = []

    directory, outdir, outfile = parse_arguments()
    genomes,genomeFiles = getGenomes(directory)

    #get substrate information
    primary_ce = set(['CHEBI:33696','CHEBI:33838','CHEBI:36976','CHEBI:18282','CHEBI:18059','CHEBI:22563',
                    'CHEBI:36916','CHEBI:33575','CHEBI:16670','CHEBI:35381','CHEBI:50699','CHEBI:18154','CHEBI:72813'
                    'CHEBI:88061','CHEBI:10545','CHEBI:25106'])

    secondary_ce = set(['CHEBI:24870','CHEBI:25367','CHEBI:83821'])

    primary_role = set(['CHEBI:33281','CHEBI:26672','CHEBI:31432','CHEBI:33229'])

    secondary_role = set(['CHEBI:23888','CHEBI:25212','CHEBI:23357'])

    substrate_data = get_substrate_data('https://tcdb.org/cgi-bin/substrates/getSubstrates.py')


    for genome in genomes:

        lines = readLines(genomeFiles[genome])

        tcidMaps[genome],tcids,accessions,tms,repeats = mapTCID(lines,tcids,accessions,tms,genome,repeats)

    tcids = sorted(list(tcids),key=lambda x: x.split('.'))

    printTable(genomes,tcids,tcidMaps,accessions,tms,substrate_data,primary_ce,secondary_ce,primary_role,secondary_role,repeats,outdir,outfile)
