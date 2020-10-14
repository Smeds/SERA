#!/usr/bin/python2.7

import argparse
import re
import sys
from os import walk
import pprint
from collections import OrderedDict

# Parse commandline
parser = argparse.ArgumentParser(description = "")
parser.add_argument('-f', '--file',
                    help = 'File(s) that will be processed (ClinicalVariant/ClinicalPosition)',
                    action='append',
                    default = [],
                    required = True)
parser.add_argument('-c', '--codons',
                    help = 'File with codons that should be extracted (Gene name\tcondon nr',
                    type = str,
                    required = False)

args = parser.parse_args()

#Variable that will store csv information
data = {}

def columnsToDict(headerMap, columns):
    """
    Creates a dict with information from a row in a csv file. It will use a header dict to
     extract information from the correct column.

    :param headerMap a dict mapping column name to colum position, ex {"#Gene": 0, "Found": 1}
    :param columns: row from csv file
    :return: a dict with data

    >>> pprint.pprint(columnsToDict({'#Gene': 0, \
                        'Cosmic_id': 1, \
                        'Report': 2, \
                        'Found': 3, \
                        'Min_read_depth30': 4, \
                        'Min_read_depth100': 5, \
                        'Min_read_depth300': 6, \
                        'Min_read_depth500': 7, \
                        'Total_read_depth': 8, \
                        'Chr': 9, \
                        'Start': 10, \
                        'End': 11, \
                        'Reference': 12, \
                        'Variant': 13, \
                        'Variant_read_depth': 14, \
                        'Reference_read_depth': 15, \
                        'Variant_allele_ratio': 16, \
                        'CDS_change': 17, \
                        'AA_change': 18, \
                        'Reference_plus_amplicons': 19, \
                        'Reference_minus_amplicons': 20, \
                        'Variant_plus_amplicons': 21, \
                        'Variant_minus_amplicons': 22, \
                        'Transcript': 23, \
                        'Protein': 24, \
                        'Ref_amplicons': 25, \
                        'Var_amplicons': 26}, \
                        ['NRAS','-','2','no','ok','ok','ok','ok','5768','NC_000001.10','115252202','115252203','GG', \
                         'AA','-','-','-', 'c.437_438CC>TT','p.A146V','3_3','3_3','-','-', 'NM_002524.4', \
                         'NP_002515.1','chr1:115252142-115252239:-:31','-']))
    {'#Gene': 'NRAS',
     'AA_change': 'p.A146V',
     'CDS_change': 'c.437_438CC>TT',
     'Chr': 'NC_000001.10',
     'Cosmic_id': '-',
     'End': '115252203',
     'Found': 'no',
     'Min_read_depth100': 'ok',
     'Min_read_depth30': 'ok',
     'Min_read_depth300': 'ok',
     'Min_read_depth500': 'ok',
     'Protein': 'NP_002515.1',
     'Ref_amplicons': 'chr1:115252142-115252239:-:31',
     'Reference': 'GG',
     'Reference_minus_amplicons': '3_3',
     'Reference_plus_amplicons': '3_3',
     'Reference_read_depth': '-',
     'Report': '2',
     'Start': '115252202',
     'Total_read_depth': '5768',
     'Transcript': 'NM_002524.4',
     'Var_amplicons': '-',
     'Variant': 'AA',
     'Variant_allele_ratio': '-',
     'Variant_minus_amplicons': '-',
     'Variant_plus_amplicons': '-',
     'Variant_read_depth': '-'}
    """

    return {'#Gene': columns[headerMap['#Gene']],
            'Cosmic_id': columns[headerMap['Cosmic_id']],
            'Report': columns[headerMap['Report']],
            'Found': columns[headerMap['Found']],
            'Min_read_depth30': columns[headerMap['Min_read_depth30']],
            'Min_read_depth100': columns[headerMap['Min_read_depth100']],
            'Min_read_depth300': columns[headerMap['Min_read_depth300']],
            'Min_read_depth500': columns[headerMap['Min_read_depth500']],
            'Total_read_depth': columns[headerMap['Total_read_depth']],
            'Chr': columns[headerMap['Chr']],
            'Start': columns[headerMap['Start']],
            'End': columns[headerMap['End']],
            'Reference': columns[headerMap['Reference']],
            'Variant': columns[headerMap['Variant']],
            'Variant_read_depth': columns[headerMap['Variant_read_depth']],
            'Reference_read_depth': columns[headerMap['Reference_read_depth']],
            'Variant_allele_ratio': columns[headerMap['Variant_allele_ratio']],
            'CDS_change': columns[headerMap['CDS_change']],
            'AA_change': columns[headerMap['AA_change']],
            'Reference_plus_amplicons': columns[headerMap['Reference_plus_amplicons']],
            'Reference_minus_amplicons': columns[headerMap['Reference_minus_amplicons']],
            'Variant_plus_amplicons': columns[headerMap['Variant_plus_amplicons']],
            'Variant_minus_amplicons': columns[headerMap['Variant_minus_amplicons']],
            'Transcript': columns[headerMap['Transcript']],
            'Protein': columns[headerMap['Protein']],
            'Ref_amplicons': columns[headerMap['Ref_amplicons']],
            'Var_amplicons': columns[headerMap['Var_amplicons']]}

codons = {}
if args.codons is not None:
    with open(args.codons, 'r') as c:
        for line in c:
            line = line.rstrip()
            columns = line.split('\t')
            if columns[0] in codons:
                codons[columns[0]]['AA_change'].append(columns[1])
                codons[columns[0]]['CDS_change'].append(columns[2])
            else:
                codons[columns[0]] = {'AA_change': [columns[1]], 'CDS_change': [columns[2]]}

#Process each provided file
for f in args.file:
    #Extract sample name
    sample = re.search('\/*([A-Za-z0-9]+)\.h\.sapiens\.clinicalPositions\.txt$',f).group(1)
    data[sample] = {}
    #Process the provided file
    with open(f, 'r') as f:
        headerMap = {}
        for line in f:
            line = line.rstrip()
            columns = line.split('\t')
            #Header line found
            if line.startswith('#Gene') or line.startswith('#Sample'):
                #Create a header map
                for i, c in enumerate(columns):
                    if c.startswith('Gene'):
                        c = '#Gene'
                    headerMap[c] = i
            else:
                #Should be a data row
                aa_change = ""
                #Extract AA change codon and type
                if columns[headerMap['AA_change']] == 'p.?' or columns[headerMap['AA_change']] == '-':
                    aa_change = "Unknown"
                else:
                    aa_change = re.search('^[a-zA-Z.*]*([0-9]+|[0-9]+_[A-Z*]*[0-9]+[delins]*)[0-9A-Za-z*>]*$',columns[headerMap['AA_change']]).group(1)
                #Extract the CDS change
                cds_change = ""
                if columns[headerMap['CDS_change']] == '-':
                    cds_change = "Unknown"
                else:
                    cds_change = re.search('^c\.(.+)*$',columns[headerMap['CDS_change']]).group(1)

                if (not codons or #Check if filter list have been provided
                       (columns[headerMap['#Gene']] in codons and #Check that gene is in filter list
                            #Check that the AA change should be extracted or
                            (aa_change in codons[columns[headerMap['#Gene']]]['AA_change'] or
                             #that the CDS change should be extracted
                             cds_change in codons[columns[headerMap['#Gene']]]['CDS_change']))):
                    key = aa_change
                    #If no AA_change information was found the cds change should be used
                    if aa_change == "Unknown" and cds_change != "Unknown":
                        key = cds_change
                    #Add the data
                    if columns[headerMap['#Gene']] in data[sample]:
                        if (aa_change + cds_change) in data[sample][columns[headerMap['#Gene']]]:
                            data[sample][columns[headerMap['#Gene']]][key].append(columnsToDict(headerMap,columns))
                        else:
                            data[sample][columns[headerMap['#Gene']]][key] = [columnsToDict(headerMap,columns)]
                    else:
                        data[sample][columns[headerMap['#Gene']]] = {key: [columnsToDict(headerMap,columns)]}

patterm = re.compile('^.*(del|ins).*$')
printData = OrderedDict()

#Merge overlapping codon changes
for sample in data:
    for gene in data[sample]:
        for codon in data[sample][gene]:
            row_info = ""
            depth = 0
            counter = 0
            #Calculate average depth if multiple change were found for a codon
            if len(data[sample][gene][codon]) > 1:
                for row in data[sample][gene][codon]:
                    #Skip if no depth information were found
                    if not row['Total_read_depth'] == "-" :
                        row_info = gene + "\t" + codon + "\t"
                        counter = counter + 1
                        depth = depth + int(row['Total_read_depth'])
            else:
                #Skip if no depth information were found
                if not data[sample][gene][codon][0]['Total_read_depth'] == "-" :
                    row_info = data[sample][gene][codon][0]['#Gene'] + "\t" + codon + "\t"
                    depth = int(data[sample][gene][codon][0]['Total_read_depth'])
                    counter = 1
            #Save the information
            if len(row_info) > 0:
                if row_info in  printData:
                    printData[row_info][sample] = str(depth/counter)
                else:
                    printData[row_info] = OrderedDict()
                    printData[row_info][sample] = str(depth / counter)


#Print header information
sys.stdout.write("#Gere\tcodon\t")
for sample in data:
    sys.stdout.write("\t" + sample)
sys.stdout.write("\n")

#Print data for each sample and codon
for codon in printData:
    sys.stdout.write(codon)
    for sample in data:
        if sample in printData[codon]:
            sys.stdout.write("\t" + printData[codon][sample])
        else:# If no data were found for sample print a '-'
            sys.stdout.write("\t-")
    sys.stdout.write("\n")

#Perform unit testing
if __name__ == '__main__':
    import doctest
    doctest.testmod()
