#/usr/bin/python
from Parsing_functions import *
from collections import defaultdict
import itertools


#  ----------------
# | Tool functions |
#  ----------------

def fasta_iter(fasta_name):
    """ Given a fasta file, yields tuple of header, sequence. """
    with open(fasta_name) as fh:
        faiter = (x[1] for x in itertools.groupby(
            fh, lambda line: line[0] == ">"))
        for header in faiter:
            header = header.next().strip()
            seq = "".join([s.strip() for s in faiter.next()])
            yield header, seq


def subset_fasta(list_proteins, fasta_input, fasta_output_name):
    """ Takes a fasta file and creates a new one with only the sequences for which the names were given in the list. """
    kept_records = []
    with open(fasta_output_name, "w") as out:
        for header, seq in fasta_iter(fasta_input):
            if header[1:] in list_proteins:
                out.write("\n".join([header, seq]) + "\n")


def liste_from_fasta(fasta_input, striped=True):
    """ Takes a fasta file and creates a list of headers. """
    headers = []
    with open(fasta_input) as input_file:
        for line in input_file:
            if line.startswith(">"):
                if striped:
                    headers.append(line.strip().strip(">"))
                else:
                    headers.append(line.strip())
    return headers


def write_tab_output_from_dict(output_name, list_header, dico, list_keys):
    """ Takes a dictionary and a list of keys. Write a separated tab file from the dictionary keys and values. """
    with open(output_name, "w") as out:
        out.write("\t".join(list_header) + "\n")
        for clef in list_keys:
            to_write = [clef] + dico[clef]
            out.write("\t".join(to_write) + "\n")


def parse_my_output(output_name, dico):
    """This function gets the consensus value for the output of the different steps in the pipeline. If the file can not be open, NA is writtent instead. """
    try:
        input_file = open(output_name)
    except:
        for key in dico:
            dico[key].append("NA")
    else:
        with open(output_name) as input_file:
            for line in input_file:
                if line.strip() == "":
                    pass
                elif line.startswith("#"):
                    pass
                else:
                    sp = filter(None, line.strip().split())
                    dico[sp[0]].append(sp[-1])
        input_file.close()


#  ---------------------
# | Gathering functions |
#  ---------------------

def results_from_secretion(fasta_input, input_prefix, output_prefix, threshold=0.75):
    """ Uses original fasta input to create list of proteins. Uses input prefix to find and read the files created by the software used for secretion detection. Gather the results from these in one tab separated file (output_prefix.secreted.fasta) and add the consensus value (secreted or not) based on the value of threshold. """
    parsing = {'targetP': parsing_targetP,
               'signalP': parsing_signalP,
               'phobius': parsing_phobius,
               'deeploc': parsing_deeploc}

    proteins = liste_from_fasta(fasta_input)

    dict_results = defaultdict(list)
    for software in sorted(parsing.keys()):
        soft_dict = parsing[software](input_prefix)
        for protein in proteins:
            temp = "NA"
            if soft_dict:
                temp = soft_dict.get(protein, "NA")
                if temp != "NA":
                    result = temp.secretion
                else:
                    result = temp
            else:
                result = "NA"
            dict_results[protein].append(result)

    # Based on a predefined threshold for the consensus,
    # creating a list of secreted proteins
    kept_proteins = []
    for protein in proteins:
        count_secreted = float(
            sum(i == "Secreted" for i in dict_results[protein]))
        count_data = float(sum(i != "NA" for i in dict_results[protein]))

        if float(count_secreted / count_data) >= threshold:
            dict_results[protein].append("Secreted")
            kept_proteins.append(protein)
        else:
            dict_results[protein].append("Not_secreted")

    # Writing fasta with only the secreted proteins
    subset_fasta(kept_proteins, fasta_input, output_prefix + ".secreted.fasta")

    # Writing the table with the results from the secretion prediction step
    header = ["#Protein_ID"] + sorted(parsing.keys()) + ["Consensus"]
    write_tab_output_from_dict(output_prefix + ".secreted.tab",
                               header, dict_results, proteins)

    # Printing information about this step on the screen
    print "\n\n"
    print "  ************************************"
    print "   Detection of secretion is complete"
    print "  ************************************"
    print "Among", len(proteins), "proteins,", len(kept_proteins), "were predicted to be secreted."
    print "You can find the table summarizing the results in the file", output_prefix + ".secreted.tab"
    print "You can find the sequences of the kept proteins in the file", output_prefix + ".secreted.fasta"
    print ""


def results_from_TM(fasta_input, input_prefix, output_prefix, threshold=0.5):
    """ Uses original fasta input to create list of proteins. Uses input prefix to find and read the files created by the software used for TM detection. Gather the results from these in one tab separated file (output_prefix.TM.fasta) and add the consensus value (secreted or not) based on the value of threshold. """
    parsing = {'phobius': parsing_phobius,
               'deeploc': parsing_deeploc,
               'tmhmm': parsing_tmhmm}

    proteins = liste_from_fasta(fasta_input)

    dict_results = defaultdict(list)
    for software in sorted(parsing.keys()):
        soft_dict = parsing[software](input_prefix)
        for protein in proteins:
            temp = "NA"
            if soft_dict:
                temp = soft_dict.get(protein, "NA")
                if temp != "NA":
                    result = temp.TM
                else:
                    result = temp
            else:
                result = "NA"
            dict_results[protein].append(result)

    # Based on a predefined threshold for the consensus,
    # creating a list of secreted proteins
    kept_proteins = []
    for protein in proteins:
        count_secreted = float(
            sum(i == "TM" for i in dict_results[protein]))
        count_data = float(sum(i != "NA" for i in dict_results[protein]))

        if float(count_secreted / count_data) >= threshold:
            dict_results[protein].append("TM")
            kept_proteins.append(protein)
        else:
            dict_results[protein].append("Not_TM")

    # Writing fasta with only the TM proteins
    subset_fasta(kept_proteins, fasta_input, output_prefix + ".TM.fasta")

    # Writing the table with the results from the TM prediction step
    header = ["#Protein_ID"] + sorted(parsing.keys()) + ["Consensus"]
    write_tab_output_from_dict(output_prefix + ".TM.tab",
                               header, dict_results, proteins)

    # Printing information about this step on the screen
    print "\n\n"
    print "  ************************************************"
    print "   Detection of transmembrane domains is complete"
    print "  ************************************************"
    print "Among", len(proteins), "proteins,", len(kept_proteins), "were predicted to be transmembrane."
    print "You can find the table summarizing the results in the file", output_prefix + ".TM.tab"
    print "You can find the sequences of the kept proteins in the file", output_prefix + ".TM.fasta"
    print "\n"


def results_from_localisation(fasta_input, input_prefix, output_prefix):
    """ Uses original fasta input to create list of proteins. Uses input prefix to find and read the files created by the software used for cell location detection. Gather the results from these in one tab separated file (output_prefix.cell_location.fasta). The consensus value here is a little bit more complicated as it has to take more complex information into account. It also uses the TM consensus. """
    parsing = {'targetP': parsing_targetP,
               'deeploc': parsing_deeploc}

    proteins = liste_from_fasta(fasta_input)

    dict_results = defaultdict(list)
    for software in sorted(parsing.keys()):
        soft_dict = parsing[software](input_prefix)
        for protein in proteins:
            temp = "NA"
            if soft_dict:
                temp = soft_dict.get(protein, "NA")
                if temp != "NA":
                    result = temp.localisation
                else:
                    result = temp
            else:
                result = "NA"
            dict_results[protein].append(result)

    consistent_with_secretion = [
        "Cell_membrane", "Endoplasmic_reticulum", "Lysosome/Vacuole", "Extracellular"]

    dict_TM = {}
    for protein in proteins:
        dict_TM[protein] = []
    parse_my_output(input_prefix + ".TM.tab", dict_TM)

    for protein in proteins:
        R = dict_results[protein]
        corrected_result = R[0]
        if "Plastid" in R or "Chloroplast" in R:
            if "Chloroplast" in R and "Plastid" in R:
                corrected_result = "Chloroplast"
            elif "Plastid" in R and "Other" in R:
                corrected_result = "Plastid"
            else:
                corrected_result = "NA"
        elif "Mitochondrion" in R:
            if len(set(R)) > 1:
                corrected_result = "NA"
        elif "Secreted" in R:
            if dict_TM[protein] == ["TM"]:
                corrected_result = "NA"
            elif R[0] not in consistent_with_secretion:
                corrected_result = "NA"
        elif "Cell_membrane" in R:
            if dict_TM[protein] == ["Not_TM"]:
                corrected_result = "NA"
        dict_results[protein].append(corrected_result)

    # Writing the table with the results from the TM prediction step
    header = ["#Protein_ID"] + sorted(parsing.keys()) + ["Consensus"]
    write_tab_output_from_dict(output_prefix + ".cell_location.tab",
                               header, dict_results, proteins)

    # Printing information about this step on the screen
    print "\n\n"
    print "  ************************************************"
    print "   Detection of cellular location is complete"
    print "  ************************************************"
    print "You can find the table summarizing the results in the file", output_prefix + ".cell_location.tab"
    print "\n"


def combining_TM_and_secreted_fasta(secreted_fasta, TM_fasta, output_prefix):
    """ To detect effectors, we need to find the protein which have a secretion signal but are not transmembrane. This function takes the results from the secretion and TM prediction and creates a fasta with only these proteins. It also returns the name of this fasta file. """
    TM_proteins = liste_from_fasta(TM_fasta, False)
    output_name = output_prefix + ".secreted_notTM.fasta"
    with open(output_name, "w") as out:
        secreted_proteins = fasta_iter(secreted_fasta)
        for header, seq in secreted_proteins:
            if header not in TM_proteins:
                out.write(header + "\n")
                out.write(seq + "\n")
    return output_name


def results_from_effectors(fasta_input, input_prefix, output_prefix, tab_output=True):
    """ Gathers the results from the the software effectorP, apoplastP and LOCALIZER. If tab_output == True, writes a table with the prediction from the software. Returns a dictionary with the same results. """
    parsing = {'effectorP': parsing_effectorP,
               'apoplastP': parsing_apoplastP,
               'LOCALIZER': parsing_LOCALIZER, }
    columns = ['effectorP', 'apoplastP', 'LOCALIZER']
    proteins = liste_from_fasta(fasta_input)

    dict_results = defaultdict(list)
    effector_nb = 0
    for software in columns:
        soft_dict = parsing[software](input_prefix)
        for protein in proteins:
            if soft_dict:
                result = soft_dict.get(protein, "-")
                if result == "Effector":
                    effector_nb += 1
            else:
                result = "NA"
            dict_results[protein].append(result)

    # Outputing the results from the effector prediction step
    print "\n\n"
    print "  *************************************"
    print "   Detection of effectors is complete"
    print "  *************************************"
    print "Among", len(proteins), "proteins,", effector_nb, "were predicted to be effectors."
    if tab_output:
        write_tab_output_from_dict(output_prefix + ".effectors.tab",
                                   ["#Protein_ID"] + columns, dict_results, proteins)
        print "You can find the table summarizing the results in the file", output_prefix + ".effectors.tab\n"

    print "\n"
    return dict_results


def all_results(fasta_input, input_prefix, output_prefix):
    """ This will gather the results from all of the previous steps in one large tab delimited file. """
    dict_all = {}
    proteins = liste_from_fasta(fasta_input)
    for protein in proteins:
        dict_all[protein] = []

    # Adding information that has already been gathered in previous steps
    secreted_file_name = input_prefix + ".secreted.tab"
    parse_my_output(secreted_file_name, dict_all)

    TM_file_name = input_prefix + ".TM.tab"
    parse_my_output(TM_file_name, dict_all)

    cell_loc_file_name = input_prefix + ".cell_location.tab"
    parse_my_output(cell_loc_file_name, dict_all)

    # Adding information that has not yet been gathered
    dict_cluster = parsing_antismash(input_prefix + ".antismash.out")
    dict_effectors = results_from_effectors(
        fasta_input, input_prefix, output_prefix)
    for protein in proteins:
        dict_all[protein].append(dict_cluster.get(protein, "-"))
        dict_all[protein].extend(dict_effectors.get(protein, ["-", "-", "-"]))

    header = ["#Protein_ID", "Secretion",
              "Transmembrane", "Cell_location", "Cluster(nb, function)", "Effector", "Apoplast", "LOCALIZER"]
    write_tab_output_from_dict(output_prefix + ".all_predictions.tab",
                               header, dict_all, proteins)

    print "\n"
    print "  *************************************"
    print "     The whole pipeline is complete"
    print "  *************************************"
    print "You can find the table summarizing the results in the file", output_prefix + ".all_predictions.tab\n\n"
