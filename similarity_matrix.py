'''
+++++++++++++++ Imports ++++++++++++++++++
'''
import os
import argparse
import parasail
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
 
description = ''' 

███████╗██╗███╗   ███╗██╗██╗      █████╗ ██████╗ ██╗████████╗██╗   ██╗    ███╗   ███╗ █████╗ ████████╗██████╗ ██╗██╗  ██╗
██╔════╝██║████╗ ████║██║██║     ██╔══██╗██╔══██╗██║╚══██╔══╝╚██╗ ██╔╝    ████╗ ████║██╔══██╗╚══██╔══╝██╔══██╗██║╚██╗██╔╝
███████╗██║██╔████╔██║██║██║     ███████║██████╔╝██║   ██║    ╚████╔╝     ██╔████╔██║███████║   ██║   ██████╔╝██║ ╚███╔╝ 
╚════██║██║██║╚██╔╝██║██║██║     ██╔══██║██╔══██╗██║   ██║     ╚██╔╝      ██║╚██╔╝██║██╔══██║   ██║   ██╔══██╗██║ ██╔██╗ 
███████║██║██║ ╚═╝ ██║██║███████╗██║  ██║██║  ██║██║   ██║      ██║       ██║ ╚═╝ ██║██║  ██║   ██║   ██║  ██║██║██╔╝ ██╗
╚══════╝╚═╝╚═╝     ╚═╝╚═╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝   ╚═╝      ╚═╝       ╚═╝     ╚═╝╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝  ╚═╝


    Hello, bioinformatician! :3

    This program calculates the protein similarity matrix. To do this, a directory with protein sequences in fasta format should be provided as input to the program, using the flag -i. Then, i) a pairwise alignment is performed between all sequence files; ii) similarity is calculated; iii) the values are stored in a matrix, and finally, iv) the program returns a CSV file with the similarity matrix values and a graphical representation of the matrix, a heatmap in PNG format. This analysis can be used to quickly visualize clusters of proteins with high similarity. 
'''

'''
+++++++++++++++ Menu ++++++++++++++++++
'''

# Function to create the initialize menu using argparse
def argument_parsing():
    parser = argparse.ArgumentParser(prog='Similarity Matrix',
                                     description=description)


    parser.add_argument('-i', '--input-directory',
                        dest='input_directory',
                        type=str,
                        metavar='',
                        help='Defines the directory path with the sequence files to be analysed')

    parser.add_argument('-v', '--verbose',
                        dest='verbosity',
                        action='count',
                        default=0,
                        help='Controls verbosity')
    initial_args = parser.parse_args()

    return initial_args

# Parsing the initial arguments
initial_args = argument_parsing()

print(description)

try:
    # Set the working directory
    work_dir = os.path.join(initial_args.input_directory)

    # Get a list of all the FASTA files in the directory
    if initial_args.verbosity > 0:
        print("Getting all the FASTA files in the provided directory")
    fasta_files = [f for f in os.listdir(work_dir) if f.endswith(".fasta")]

    # Set the alignment parameters 
    matrix = parasail.blosum62
    gap_open = 10
    gap_extend = 1

    # Initialize a dataframe to store the results
    names =  [name.split('_')[0] for name in fasta_files]
    df = pd.DataFrame(index=names, columns=names)

    # Loop over the FASTA files performing pairwise alignment
    for i in range(len(fasta_files)):
        for j in range(len(fasta_files)):
            # Read in the sequences from the FASTA files
            with open(os.path.join(work_dir, fasta_files[i]), "r") as f1:
                seq1 = f1.read().strip()
            with open(os.path.join(work_dir, fasta_files[j]), "r") as f2:
                seq2 = f2.read().strip()
            # Perform the pairwise alignment
            if initial_args.verbosity > 0:
                print( "Performing parwise alignment between" + " " + str(fasta_files[i]) + " " + "and" + " " + str(fasta_files[j]))
            alignment = parasail.sg_stats_striped_16(seq1, seq2, gap_open, gap_extend, matrix)
            # Calculate the identity score in percentage
            percent_id = round(((alignment.matches)/alignment.length*100), 2)
            if initial_args.verbosity > 0:
                print("Alignment identity score percentage: " + str(percent_id))
            # Add the results to the dataframe
            df.loc[fasta_files[i].split("_")[0],fasta_files[j].split("_")[0]]= percent_id

    # Export matrix in a csv file
    df.to_csv('similarityMatrix.csv', sep=',')

    '''
    +++++++++++++++ Heatmap Plot ++++++++++++++++++
    '''

    #Set canva theme
    sns.set_theme()

    # Convert the palette to vectors that will be drawn on the side of the matrix
    df = df.fillna(0)
    proteins = df.columns.get_level_values(0)

    # Draw the full plot
    g = sns.clustermap(df.corr(), center=0, cmap = sns.diverging_palette(230, 0, 90, 60, as_cmap=True),
                    dendrogram_ratio=(.1, .2),
                    cbar_pos=(.02, .32, .03, .2),
                    linewidths=.75, figsize=(12, 13))

    g.ax_row_dendrogram.remove()

    g.savefig("similarityMatrix.png")

except:
    print("Provide a directory with protein sequences in fasta format to perform the analysis :3")