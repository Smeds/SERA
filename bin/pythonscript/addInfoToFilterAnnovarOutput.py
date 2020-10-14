#!/usr/bin/python2.7

import argparse

# Parse commandline
parser = argparse.ArgumentParser(description = "")
parser.add_argument('-fa', '--filteredannovar',
                    help = 'Filtered annovar output file',
                    type = str,
                    required = True)
parser.add_argument('-i', '--info',
                    help = 'File with information that should be added, a csv file. It must contain a column named Projektnr',
                    type = str,
                    required = True)

def generate_header_dict(line, required_column_name = 'Projektnr'):
    if not required_column_name in line:
        raise Exception('No valid header found, column with name', required_column_name, 'is missing.')
    else:
        return dict((v, i) for i, v in enumerate(line.rstrip().split('\t')))

args = parser.parse_args()

#
# Start by importing the information that will be added to filter annovar file.
#

# Variable used to store data
data = {}
headerString = ""
if args.info is not None:
    with open(args.info, 'r') as info:
        # populate the header dict
        headerString = info.readline().rstrip();
        header = generate_header_dict(headerString, 'Projektnr')
        # Process the information found in the file.
        for line in info:
            line = line.rstrip()
            data[line.split('\t')[header['Projektnr']]] = line

#
# Process the filter annovar file.
#

with open(args.filteredannovar, 'r') as fa:
    # Read header line
    line = info.readline()
    # Create a dict with name -> position information
    header = generate_header_dict(line, 'Sample')
    # Remove the sample nr from the filter annovar file, it already exists in the headerString variable
    columns = line.rstrip.split('\t')
    del columns[header['Sample']]
    print headerString + "\t".join(columns)
    # Process the data
    for line in fa:
        columns = line.split('\t')
        sample = columns[header['Sample']]
        if sample in data:
            del columns[header['Sample']]
            print data[sample] + "\t" + "\t".join(columns)
        else:
            raise Exception('Sample ', columns[header['Sample']] , " not found in information file")

# Perform unit testing
if __name__ == '__main__':
    import doctest
    doctest.testmod()
