from Bio.Blast import NCBIXML
import argparse
import time

#### ARGS DEFINITIONS
from bin.tax import TaxMatch,TaxTree

parser = argparse.ArgumentParser(description='Augmenting BLAST results.')
parser.add_argument("--input", "-i",
                    help=("The path of the BLAST XML results."))
parser.add_argument("--out", "-o",
                    help=("The path for the processed results"))
parser.add_argument("--verbose", "-v",
                    help=("Print complete output in the terminal."),
                    action="store_true")
parser.add_argument("--ascii", "-at",
                    help=("Print an ASCII representation of the processed results in terminal."
                         "Will require use of taxonomy database matching and may therefore take a long time for large result sets."),
                    action="store_true")
parser.add_argument("--save_ascii", "-sa",
                    help=("Save the ASCII representation of the taxonomy tree to a specified file or enter 'default' to append to the end of existing output file."))
parser.add_argument("--most", "-m", nargs='?',
                    help=("List the most commonly occurring matches. Format: [Rank][Number of results] e.g: -m 'Genus 5 Family 5 Class 5'"
                         "Accepts: Species Genus Family Order Class Phylum Kingdom Superkingdom or All"
                         "Will require use of taxonomy database matching and may therefore take a long time for large result sets."),
                    action="store",
                    default=None)
parser.add_argument("--least", "-l", nargs='?',
                    help=("List the least commonly occurring matches. Format: [Rank][Number of results] e.g: -l 'Genus 5 Family 5 Class 5'"
                         "Accepts: Species Genus Family Order Class Phylum Kingdom Superkingdom or All"
                         "Will require use of taxonomy database matching and may therefore take a long time for large result sets."),
                    action="store",
                    default = None)
parser.add_argument("--visualise","-vt",
                    help="Show a GUI representation of the taxonomy tree after the results have been processed.",
                    action="store_true")

###METHODS

"""Takes in the filepath of the BLAST XML file based on args provided, processes it and, if successful" returns an XML parser. """

def open_BLAST_XML_file(args):
    if args.input:
        results_handle = open(args.input)
    else:
        print('Error: Input not specified')
        exit()

    blastRecords = NCBIXML.parse(results_handle)
    return blastRecords

"""Incredibly basic tree filepath substitution. """

def check_tree_args(args):
    if args.save_ascii == "default":
        if args.out is None:
            print("--save_ascii error: No default output defined.")
            exit()
        else:
            args.save_ascii = args.out

"""Basic check for the correctness of indexing arguments. """

def check_indexing_args(argument):

    # All of this is incredibly prone to breaking and abuse,
    # but frankly it's a whole other project writing a 'perfect' args parser.

    if argument is None:
        return argument

    ranks = ['Species', 'Genus', 'Family', 'Order',
             'Class', 'Phylum', 'Kingdom', 'Superkingdom']
    list = argument.split()

    if len(list) % 2 is not 0:
        print("Least/Most Common error: incorrect number of arguments.")
        exit()

    if list[0] == "All":
        list = ['Species', list[1], 'Genus', list[1], 'Family', list[1], 'Order', list[1],
             'Class', list[1], 'Phylum', list[1], 'Kingdom', list[1], 'Superkingdom', list[1]]

    correct_list = []

    while len(list) > 0:
        rank = list.pop(0)
        number = list.pop(0)
        if rank not in ranks:
            print("Least/Most Common Error: incorrect rank", rank)
            exit()
        else:
            correct_list.append(rank)
        if not number.isdigit():
            print("Least/Most Common Error: incorrect value", number)
            exit()
        else:
            correct_list.append(number)

    return " ".join(correct_list)

"""Calculates the percentage match of a given hit. """

def calculate_match_percentage(hsp):
    result = (hsp.identities / hsp.align_length) * 100
    return "{:.3f}".format(result)

"""Calculates the amount of mismatches in a hit. """

def calculate_mismatches(hsp):
    result = hsp.align_length - hsp.identities - hsp.gaps
    return result

"""Calculates the total number of open gaps between the query and the subject """

def calculate_gap_opens(hsp):
    if hsp.gaps == 0:
        return 0
    else:
        return count_gaps_open(hsp.query) + count_gaps_open(hsp.sbjct)

"""Counts the amount of open gaps within an individual sequence """

def count_gaps_open(sequence):
    last = ""
    gaps_open = 0
    for letter in sequence:
        if letter == "-":
            if last != "-":
                gaps_open += 1
        last = letter
    return gaps_open


""" Returns the total alignment coverage based on the amount of identities matched, relative to the size of a query.
"""

def calculate_total_alignment_coverage(blast_record, hsp):
    coverage = (hsp.identities / blast_record.query_length) * 100
    return "{:.3f}".format(coverage)

"""A wrapper function for counting individual gaps. Takes in a given HSP and its respective subject/query sequence """

def get_gap_details(hsp, sequence):
    gaps, starts = count_individual_gaps(sequence)
    reverse_flag = (hsp.query_start - hsp.query_end) > 0
    gap_result = show_individual_gaps(gaps, starts, hsp.query_start, reverse_flag)

    return gap_result


"""Iterates over an individual sequence, counting gap lengths and their respective positions.

    Returns a list of gap sizes, followed by a list of their starting positions. 
"""

def count_individual_gaps(sequence):
    last = ""
    index = 0
    gaps = []
    start = []
    for letter in sequence:
        index += 1
        if letter == "-":
            if last != "-":
                gaps.append(1)
                start.append(index)
            else:
                gaps[-1] += 1
        last = letter

    return gaps, start


"""Formats the information returned from count_individual_gaps above into a single string,
    in the format of [gap length]:[gap location]/[gap length]:[gap location]/..

    Arguments:
    gap_sizes -- list of single [int] values representing the length of an individual gap
    gap_starts -- list of single [int] values representing the position of an individual gap
    start_position -- the [int] starting point of the subject sequence matched
    reverse -- a boolean value indicating whether the subject sequence was matched as in a typical
               left-to-right fashion or in reverse (e.g start: 50, end: 0) 
"""

def show_individual_gaps(gap_sizes, gap_starts, start_position, reverse):
    g = []
    if len(gap_sizes) >= 1:
        for gaps_size, gap_start in zip(gap_sizes, gap_starts):

            g.append(str(gaps_size))
            g.append(':')

            if (reverse):
                position = start_position - gap_start
            else:
                position = start_position + gap_start

            g.append(str(position))
            g.append("/")

        g.pop()
        g = ''.join(g)


    #  elif len(gap_sizes) == 1:
    #      g.append(str(gap_sizes[0]))
    #     g.append(':')
    #      if (reverse):
    #          g.append(str(start_position - gap_starts[0]))
    #      else:
    #          g.append(str(start_position + gap_starts[0]))
    #      g = ''.join(g)
    else:
        g = 'None'

    return g

"""A wrapper function for printing out query coverage and gap information. Takes the current BLAST record and HSP being processed """

def get_gaps_and_coverage(blast_record, hsp):
    if hsp.gaps != 0:
        query_gaps = get_gap_details(hsp,hsp.query)
        subject_gaps = get_gap_details(hsp,hsp.sbjct)

        output = ('Individual gaps: Q:', query_gaps, " S:", subject_gaps, 'Query-Subject Coverage:',
                  calculate_total_alignment_coverage(blast_record, hsp))
        output = ' '.join(output)
        output += '\n'
    else:
        output = ('Query-Subject Coverage:', calculate_total_alignment_coverage(blast_record, hsp))
        output = ' '.join(output)
        output += '\n'

    return output

"""A wrapper function for compiling the tabular BLAST output. Takes in current blast_record, alignment and HSP """

def compile_tab_output(blast_record, alignment, hsp):
    tab_output = (blast_record.query,
                  alignment.title,
                  calculate_match_percentage(hsp),
                  hsp.align_length,
                  calculate_mismatches(hsp),
                  calculate_gap_opens(hsp),
                  hsp.query_start,
                  hsp.query_end,
                  hsp.sbjct_start,
                  hsp.sbjct_end,
                  hsp.expect,
                  round(hsp.bits))

    tab_output = '\t'.join(map(str, tab_output))
    tab_output += '\n'

    return tab_output

"""The main BLAST XML iterator function. """

def blast_parse(blast_records):

    blast_total_processed = 0
    blast_results_matched = 0
    blast_individual_hits = 0

    match_names = []

    if args.out: file = open(args.out, 'w')

    start = time.time()

    for blast_record in blastXMLRecords:
        blast_total_processed += 1
        for alignment in blast_record.alignments:
                blast_results_matched += 1
                for hsp in alignment.hsps:
                    blast_individual_hits += 1

                    tab = compile_tab_output(blast_record, alignment, hsp)
                    print(tab)

                    stats = get_gaps_and_coverage(blast_record, hsp)
                    print(stats)
                    match_names.append(alignment.title)

                    if args.out:
                        file.write(tab)
                        file.write(stats)


    end = time.time()
    print("BLAST parsing finished in ", end - start)
    print("BLAST results processed:", blast_total_processed)
    print("BLAST results with hits:", blast_results_matched)
    print("BLAST total hits:", blast_individual_hits)

    return match_names

################
# PROGRAM MAIN #
################


# PARSING CHECKS

start = time.time()

# Initial safety checks.
args = parser.parse_args()
check_tree_args(args)
args.most = check_indexing_args(args.most)
args.least = check_indexing_args(args.least)

blastXMLRecords = open_BLAST_XML_file(args) #BLAST Output + extra stats
match_names = blast_parse(blastXMLRecords)

if(args.most or args.least or args.visualise or args.ascii or args.save_ascii): #last second fix, do refactor

    species_indexes = TaxMatch.match_results_to_taxonomy(match_names) # Matching for most/least & tree work.
    TaxTree.main(species_indexes, args)
# Sample tree used in the unit tests.
#species_indexes = ["9606","63221","741158","11676","7157","842271","1350382","1289305","1427098","1042"]
# There may come a time that tree functionality will be expanded to the point of needing to generate new pickles.
# Use the above to generate them.
#



end = time.time()
print("Finished in", end - start)

#Personal reference for Tabular Output
#Query ID, Subject ID, Percentage of Identical Matches, Alignment Length, Number Of Mismatches, Number Of Gap Openings,
#Start Of Alignment Query, End Of Alignment Query, Start Of Alignment Subject, End Of Alignment Subject, Expect Value, Bit Score
# also known as:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore