#/usr/bin/python


#	***************************************************
# |   This script was written by Alice Feurtey in 2018	|
# |   It was last modified on the 28/01/19				|
#	***************************************************


import os
import sys
import argparse
import itertools
from collections import namedtuple, defaultdict
from Parsing_functions import *
from Gathering_functions import *
from Input_functions import *
from difflib import SequenceMatcher

#   <<>><<>><<>>
# |   WARNINGs	|
#   <<>><<>><<>>


# I have modified some of the software to parse their outputs more conveniently.
#
# On line 47 of the file deeploc, I have change the opening of the file to "a" instead of
# "w". On line 48, I have added a "#" at the beginning of the header line.
#
# In the file targetP and phobius, I have changed the output to add a # in front of
# every line which is not a per-gene result.


# Default paths for wallace cluster. Change to match yours.
paths = {"targetP": "/home/feurtey/Documents/Software/targetp-1.1/targetp",
         "signalP": "/home/feurtey/Documents/Software/signalp-4.1/signalp",
         "phobius": "/home/feurtey/Documents/Software/phobius/phobius.pl",
         "DeepLoc": "deeploc",
         "TMHMM": "/home/feurtey/Documents/Software/tmhmm-2.0c/bin/tmhmm",
         "EffectorP": "/home/feurtey/Documents/Software/EffectorP_2.0/Scripts/EffectorP.py",
         "ApoplastP": "/home/feurtey/Documents/Software/ApoplastP_1.0.1/Scripts/ApoplastP.py",
         "LOCALIZER": "/home/feurtey/Documents/Software/LOCALIZER_1.0.4/Scripts/LOCALIZER.py",
         "antismash": "/home/feurtey/bin/antismash-4.2.0/run_antismash.py"}


# Default options for Alice's project. Change to match yours, if necessary.
options = {"targetP": "",
           "signalP": " -t euk",
           "phobius": " -short",
           "DeepLoc": "",
           "TMHMM": "",
           "EffectorP":  "",
           "ApoplastP": "",
           "LOCALIZER": " -e -S 20",
           "antismash": "--taxon fungi --input-type nucl  --minlength 200"}


# These are the command lines. It is very probable that you do not have to touch them. If you do, make sure to keep the number and order of elements. This is important, since the input name will be inserted between the first and second element and the output name between the second and third.
commands = {"targetP": [paths["targetP"] + " " + options["targetP"] + " -N ", " >> ", ".targetP.out"],
            "signalP": [paths["signalP"] + " " + options["signalP"] + " ", " >> ", ".signalP.out"],
            "phobius": [paths["phobius"] + " " + options["phobius"] + " ", " >> ", ".phobius.out"],
            "TMHMM": [paths["TMHMM"] + " " + options["TMHMM"] + " ", " >> ", ".tmhmm.out"],
            "DeepLoc": [paths["DeepLoc"] + " -f ", " -o ", ".deeploc.out", " " + options["DeepLoc"]],
            "EffectorP": ["python " + paths["EffectorP"] + " -i ", " -o ", ".EffectorP.out", " " + options["EffectorP"]],
            "ApoplastP": ["python " + paths["ApoplastP"] + " -i ", " -o ", ".ApoplastP.out", " " + options["ApoplastP"]],
            "LOCALIZER": ["python " + paths["LOCALIZER"] + " -i ", " -o ", ".LOCALIZER.out", " " + options["LOCALIZER"]],
            "antismash": ["python " + paths["antismash"] + " --gff3 ", " --outputfolder ", ".antismash.out", " " + options["antismash"]]}

all_secretion_software = {"targetP", "signalP", "phobius", "DeepLoc"}
all_TM_software = {"phobius", "DeepLoc", "TMHMM"}
all_location_software = {"targetP", "DeepLoc"}
all_effector_software = {"EffectorP", "ApoplastP", "LOCALIZER"}
all_cluster_software = {"antismash"}


#   <<>><<>><<>><<>><<>><<>>
# |   Functions & classes	|
#   <<>><<>><<>><<>><<>><<>>


def fasta_iter_and_count(fasta_name):
    """ Given a fasta file, yields tuple of count, header, sequence"""
    fh = open(fasta_name)
    faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == ">"))
    count = 0
    for header in faiter:
        count += 1
        header = header.next().strip()
        seq = "".join([s.strip() for s in faiter.next()])
        yield count, header, seq  # modified


def check_length_fasta_header_and_find_common_substr(fasta_name, max_length=0):
    too_long_names = []
    for count, header, seq in fasta_iter_and_count(fasta_name):
        # The plus one below accounts for the '>' which is included in the header.
        if len(header) >= (max_length + 1):
            if len(too_long_names) <= 25:
                print "BEWARE! This fasta header is longer than targetP likes:", header, " This will cause this protein to not be recognized later on."
            if len(too_long_names) == 25:
                print "I have found more than 25 headers which are deemed too long. I will stop printing them on the screen."
            too_long_names.append(header[1:])

    if len(too_long_names) > 1:
        substr0 = too_long_names[0]
        for substr1 in too_long_names[1:]:
            match = SequenceMatcher(None, substr0, substr1).find_longest_match(
                0, len(substr0), 0, len(substr1))
            if match.size < 3:
                break
            else:
                substr0 = substr0[match.a:(match.a + match.size)]
        print "\nI have found", len(too_long_names), "headers which are longer than liked by targetP."
        return True, substr0
    else:
        return False, ""


def check_if_substring_in_seq(fasta_name, substring):
    is_in = False
    for count, header, seq in fasta_iter_and_count(fasta_name):
        if substring in seq:
            is_in = True
            break
    return is_in


def run_software(fasta_name, output_prefix, list_software, run=True):
    for soft in list_software:
        command = commands[soft][:]
        command.insert(1, fasta_name)
        command.insert(3, output_prefix)
        to_write = "".join(command)
        if run:
            print "\nRunning", soft
            try:
                os.system(to_write)
            except:
                print "This command failed", to_write
        else:
            print "".join(["   - ", soft, ": ", to_write])


def run_antismash(gff_name, genome_fasta_name, output_prefix, run=True):
    """ Currently, this function is used to run antismash only, which is currently the only one with two input files. At the """
    input_string = " ".join([gff_name, genome_fasta_name])
    run_software(input_string, output_prefix, ["antismash"], run)


def effector_prediction(fasta_name, output_prefix, list_software):

    if "EffectorP" in list_software:
        run_software(fasta_name, output_prefix, ["EffectorP"])

    result_effectors = parsing_effectorP(output_prefix)
    effector_list = []
    for protein in result_effectors:
        if result_effectors[protein] == "Effector":
            effector_list.append(protein)
    subset_fasta(effector_list, fasta_name, output_prefix + ".effectors.fasta")

    list_software_loc = [x for x in list_software if x != "EffectorP"]
    run_software(output_prefix + ".effectors.fasta",
                 output_prefix, list_software_loc)


#   <<>><<>><<>><<>><<>><<>><<>>
# |  Inputs, info and warnings  |
#   <<>><<>><<>><<>><<>><<>><<>>


#   -----------------------------
#  |   Inputs from command line  |
#   -----------------------------
parser = argparse.ArgumentParser(
    description=("WARNING: I have modified some of the software to parse their outputs more conveniently.On line 47 of the file deeploc, I have change the opening of the file to 'a' instead of	'w'. On line 48, I have added a '#' at the beginning of the header line. In the file targetP, I have changed the output to add a '#' in front of	every line which is not a per-gene result."))
group_inputs = parser.add_argument_group(
    'Input arguments', 'These are the input files. You must provide either 1. a fasta file containing protein sequences, 2. a fasta file containing the genome AND an annotation file (GFF or GTF) containing the annotation for the exons of each gene or 3. all of the above. If you do not provide a protein sequence fasta, I will generate my own. ')
group_inputs.add_argument(
    "-p", "--protein_fasta_input", type=str, help="This is a protein fasta file on which the predictions will be done. This is one option for the input files. This will disable antismash which can only run from genomic sequences and annotation file in this pipeline. The corresponding column in the final table will be empty unless you provide the antismash results (the file gene_clusters.txt specifically) in a folder named output_prefix.antismash.out ")
group_inputs.add_argument(
    "-g", "--genomic_fasta_input", type=str, help="This is the genome in a multi-fasta format. (Necessary with -a)")
group_inputs.add_argument(
    "-a", "--annotation_input", type=str, help="This is a gene annotation file. This should be in a GFF or GTF format. (Necessary with -g)")
group_inputs.add_argument("--which_type", type=str,
                          help="This is the type of the annotation to create the protein sequences from. The default is CDS. (Only with -g and -a)", default="CDS")
group_inputs.add_argument("--which_attribute", type=str,
                          help="This is the name of the attributes to create the protein sequences from. The default is Parent. (Only with -g and -a)", default="Parent")

group_outputs = parser.add_argument_group('Output arguments', '')
group_outputs.add_argument("--check_command_lines_only", action="store_true")
group_outputs.add_argument("--check_protein_names_only", action="store_true",
                           help="I have found that targetP does not like long protein names. It will create a table with the names shortened and thus impossible to recognized and distinguish from one another. This is to check that your proteins all have names compatible with targetP.")
group_outputs.add_argument(
    "-o", "--output", help="This is a prefix for all of the output files. The default is to follow the format of the fasta input (either the genome or the protein one).")
group_outputs.add_argument("-t", "--temp_directory", help=" In case the number of protein you are interested is too high, I will split the file and work chunk by chunk. I need a temporary directory for this. The default is to simply write the temp files to the current directory.", default="./")

group_software = parser.add_argument_group('Software and steps arguments', '')
group_software.add_argument(
    "--software", help="The supported software so far are: targetP, signalP, phobius, DeepLoc, TMHMM, EffectorP, ApoplastP and LOCALIZER. Giving this option a random string, such as 'None' will lead this script to only gather the results from the outputs of the different software. The output need to follow the specific format of this script, though.", nargs="+")
group_software.add_argument("--secretion_threshold", default="0.75",
                            help="Proportion of software identifying a secreted protein (including the threshold value)", type=float)
group_software.add_argument("--TM_threshold", default="0.75",
                            help="Proportion of software identifying a transmembrane protein (including the threshold value)", type=float)
group_ex = parser.add_mutually_exclusive_group()
group_ex.add_argument("--effector_only", action="store_true")
group_ex.add_argument("--secretion_only", action="store_true")
group_ex.add_argument("--location_only", action="store_true")
group_ex.add_argument("--TM_only", action="store_true")
group_ex.add_argument("--clusters_only", action="store_true")
A = parser.parse_args()


#   ----------------------------------------------
#  |   Inputs checks, creation and output names   |
#   ----------------------------------------------

print "\n Input/output checks and information: "
print " -------------------------------------"
fasta_extensions = ["fasta", "fa"]

# Output prefix and input presence checks
if A.output:
    output_name = A.output

if A.protein_fasta_input:
    if any([A.protein_fasta_input.endswith(end) for end in fasta_extensions]):
        if not A.output:
            output_name = A.protein_fasta_input.rsplit(".", 1)[0]
        protein_fasta_input = A.protein_fasta_input
    else:
        sys.exit("Your fasta file extension is neither fasta, nor fa. Please make sure it is indeed a fasta file. If it is, either add an output prefix in the command line using the option '-o' or change the extension of your input file.")

if A.genomic_fasta_input and A.annotation_input:
    if any([A.genomic_fasta_input.endswith(end) for end in fasta_extensions]):
        if not A.output:
            output_name = A.genomic_fasta_input.rsplit(".", 1)[0]
        print "  - Your input genome sequence file is ", A.genomic_fasta_input
        print "  - Your input annotation file is ", A.annotation_input

        # If no protein sequence fasta is provided, I will create my own.
        if not A.protein_fasta_input:
            protein_fasta_input = output_name + ".prot.fasta"
            from_annot_to_prot_fasta(protein_fasta_input, A.genomic_fasta_input,
                                     A.annotation_input, A.which_type, A.which_attribute)
            os.system(" ".join(["sed 's/\*//g' -i", protein_fasta_input]))
    else:
        sys.exit("Your fasta file extension is neither fasta, nor fa. Please make sure it is indeed a fasta file. If it is, either add an output prefix in the command line using the option '-o' or change the extension of your input file.")
elif A.genomic_fasta_input:
    sys.exit(
        "If you choose to input a genome fasta, you need to give me an annotation file too.")
elif A.annotation_input:
    sys.exit(
        "If you choose to input an annotation file, you need to give me a genome file too.")


if not any([A.annotation_input, A.genomic_fasta_input, A.protein_fasta_input]):
    print "I do need some form of input..."
    parser.print_help()
    sys.exit()


# TargetP does not like names which are too long, so I am checking for this here and will change names by removing any common string in these names.
names_too_long, longest_common_substr = check_length_fasta_header_and_find_common_substr(
    protein_fasta_input, 20)
if names_too_long:
    if A.check_protein_names_only:
        print "\nYou have chosen to check the size of your protein names only. Wise move."
        print "For your convenience, I have looked for the longest substring shared by these names. It might help you change the names."
        print "Longest substring:", longest_common_substr, "\n"
        sys.exit()

    else:
        new_protein_fasta_input = protein_fasta_input.rsplit(
            ".", 1)[0] + ".shorter_names.fasta"
        correspondance_table = protein_fasta_input.rsplit(
            ".", 1)[0] + ".correspondance_to_shorter_names.tab"
        print "Because the protein names in your input are too long for TargetP, I will change them. You can find the correspondance between the original and new names in the table", correspondance_table
        corr_out = open(correspondance_table, "w")
        corr_out.write("Original name\tShorter name\n")
        with open(new_protein_fasta_input, "w") as new_fasta_out:
            for count, header, sequence in fasta_iter_and_count(protein_fasta_input):
                shorter_header = "".join([">seq_", str(count), "\n"])
                new_fasta_out.write(shorter_header)
                new_fasta_out.write(sequence + "\n")
                corr_out.write(header + "\t" + shorter_header)
        corr_out.close()
    protein_fasta_input = new_protein_fasta_input

else:
    if A.check_protein_names_only:
        print "\nYou have chosen to check the size of your protein names only. Wise move."
        print "I have not detected anything wrong with the size of the protein names."
        sys.exit()

        # print "  - Your input protein sequence file is ", protein_fasta_input, "& it contains", nb_sequences, "sequences."
# Let's try to open the protein file and count the number of sequences
try:
    nb_sequences = 0
    with open(protein_fasta_input) as fasta_file:
        for x in fasta_file:
            if x.startswith(">"):
                nb_sequences += 1
    print "  - Your input protein sequence file is ", protein_fasta_input, "& it contains", nb_sequences, "sequences."
except IOError:
    sys.exit(
        " ".join(["Oops! The file", protein_fasta_input, "could not be opened.\n"]))
except:
    print "Oops! Unexpected error:", sys.exc_info()[0]
    raise


# If the number of sequences is too large for some of the software, inform the user that their file will be split up in chunks.
if nb_sequences > 1000:
    temp_file = A.temp_directory + "temp.fasta"
    temp = open(temp_file, "w")
    print "  - Because your input file is too large, I will use temporary files. The temporary fasta file is", temp_file

print "  - Your output files will start with", output_name
print ""

#   ---------------------------
#  |   Software for each step  |
#   ---------------------------

secretion = True
cellular_location = True
TM_detection = True
effectors = True
clusters = True

if A.secretion_only:
    cellular_location = False
    TM_detection = False
    effectors = False
    clusters = False
elif A.effector_only:
    cellular_location = False
    TM_detection = False
    secretion = False
    clusters = False
elif A.location_only:
    TM_detection = False
    effectors = False
    secretion = False
    clusters = False
elif A.TM_only:
    cellular_location = False
    effectors = False
    secretion = False
    clusters = False
elif not A.genomic_fasta_input:
    clusters = False
    if A.clusters_only:
        sys.exit("Sorry! If you want to get results from the secondary metabolite gene clusters (from antismash), you need to input a genome fasta and the corresponding annotation in gff or gtf format, NOT a protein fasta file.\n")
elif A.clusters_only:
    cellular_location = False
    effectors = False
    secretion = False
    TM_detection = False


# The software used for each steps are:
if A.software:
    all_possible_software = set(A.software)
else:
    all_possible_software = set(commands.keys())

secretion_software = set()
TM_software = set()
location_software = set()
effector_software = set()
cluster_software = set()

if effectors:
    effector_software = all_possible_software.intersection(
        all_effector_software)
if clusters:
    cluster_software = all_possible_software.intersection(all_cluster_software)
if TM_detection:
    TM_software = all_possible_software.intersection(all_TM_software)
if cellular_location:
    location_software = all_possible_software.intersection(
        all_location_software)
if secretion:
    secretion_software = all_possible_software.intersection(
        all_secretion_software)
general_software = TM_software | location_software | secretion_software

print "\n Command lines: "
print " ---------------"
run_software(protein_fasta_input, output_name, effector_software |
             general_software, run=False)
if clusters:
    run_antismash(A.genomic_fasta_input,
                  A.annotation_input, output_name, False)
print ""
if A.check_command_lines_only:
    sys.exit()


#   <<>><<>><<>>
# |  Script body |
#   <<>><<>><<>>

input_fasta = fasta_iter_and_count(protein_fasta_input)
if any([secretion, TM_detection, cellular_location]) and len(general_software) != 0:
    print "\nI will now run:\n", "\n  - " + "\n  - ".join(general_software)
    print "\n"
    if nb_sequences > 1000:
        for count, header, seq in input_fasta:
            temp.write("\n".join([header, seq]) + "\n")
            if count % 1000 == 0 or count == nb_sequences:
                temp.close()
                print "Running until sequence", count, "over", nb_sequences
                run_software(
                    temp_file, output_name, general_software)
                os.system("wc -l " + temp_file)
                if count != nb_sequences:
                    temp = open(temp_file, "w")
    else:
        run_software(protein_fasta_input, output_name, general_software)

if secretion:
    results_from_secretion(protein_fasta_input, output_name,
                           output_name, threshold=A.secretion_threshold)

if TM_detection:
    results_from_TM(protein_fasta_input, output_name,
                    output_name, threshold=A.TM_threshold)

if cellular_location:
    results_from_localisation(protein_fasta_input, output_name, output_name)


if effectors and len(effector_software) != 0:
    print "I will now run:\n", "\n  - ".join(effector_software)
    print "\n"
    if secretion and TM_detection:
        new_fasta = combining_TM_and_secreted_fasta(output_name + ".secreted.fasta",
                                                    output_name + ".TM.fasta", output_name)
        effector_prediction(new_fasta, output_name, effector_software)
    else:
        effector_prediction(protein_fasta_input,
                            output_name, effector_software)
        results_from_effectors(fasta_name, output_name,
                               output_name, True)

if clusters and len(cluster_software) != 0:
    print "I will now run:\n", "\n  - ".join(cluster_software)
    print "\n"
    run_antismash(A.annotation_input, A.genomic_fasta_input,
                  output_name, run=True)

if sum([secretion, TM_detection, effectors, clusters, cellular_location]) > 1:
    all_results(protein_fasta_input, output_name, output_name)
