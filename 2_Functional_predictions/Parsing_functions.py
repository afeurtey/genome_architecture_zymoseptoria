#/usr/bin/python
from collections import namedtuple
import sys


def check_file_opens(file_to_check):
    # Checking file
    try:
        temp = open(file_to_check)
        temp.close()
    except IOError:
        print "File", file_to_check, "could not be opened."
        return None
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise


def parsing_signalP(input_suffix):
    """Takes the output file of signalP and returns a dictionary  with a named tuple as value and the genes as keys."""
    input_name = input_suffix + ".signalP.out"
    check_file_opens(input_name)

    Results = namedtuple("Results", "secretion")
    signalP_dict = {}
    with open(input_name) as input_file:
        for line in input_file:
            if line.strip() == "":
                pass
            elif line.startswith("#"):
                pass
            else:
                sp = line.strip().split()
                if sp[9] == "N":
                    result = Results("Not_secreted")
                elif sp[9] == "Y":
                    result = Results("Secreted")
                else:
                    result = Results("NA")
                signalP_dict[sp[0]] = result
    return signalP_dict


def parsing_targetP(input_suffix):
    """Takes the output file of targetP WITH ALL THE NON RESULT LINES STARTING WITH #. And returns a dictionary  with a named tuple as value and the genes as keys."""
    input_name = input_suffix + ".targetP.out"
    check_file_opens(input_name)

    Results = namedtuple("Results", "secretion, localisation")
    targetP_dict = {}
    with open(input_name) as input_file:
        for line in input_file:
            if line.strip() == "":
                pass
            elif line.startswith("#"):
                pass
            else:
                sp = line.strip().split()
                if sp[5] == "M":
                    result = Results("Not_secreted", "Mitochondrion")
                elif sp[5] == "S":
                    result = Results("Secreted", "Secreted")
                elif sp[5] == "C":
                    result = Results("Not_secreted", "Chloroplast")
                elif sp[5] == "_":
                    result = Results("Not_secreted", "Other")
                else:
                    result = Results("NA", "Na")

                targetP_dict[sp[0]] = result
    return targetP_dict


def parsing_phobius(input_suffix):
    """Takes the output file of phobius WITH ALL THE NON RESULT LINES STARTING WITH #. And returns a dictionary  with a named tuple as value and the genes as keys. """
    input_name = input_suffix + ".phobius.out"
    check_file_opens(input_name)

    Results = namedtuple("Results", "secretion, TM, TM_domains")
    phobius_dict = {}
    with open(input_name) as input_file:
        for line in input_file:
            if line.strip() == "":
                pass
            elif line.startswith("#"):
                pass
            else:
                sp = line.strip().split()
                if sp[2] == "0":
                    is_secreted = "Not_secreted"
                elif sp[2] == "Y":
                    is_secreted = "Secreted"
                else:
                    result = Results("NA", "NA", "NA")
                if int(sp[1]) == 0:
                    is_TM = "Not_TM"
                else:
                    is_TM = "TM"
                result = Results(is_secreted, is_TM, str(sp[1]))
                phobius_dict[sp[0]] = result
    return phobius_dict


def parsing_deeploc(input_suffix):
    input_name = input_suffix + ".deeploc.out.txt"
    check_file_opens(input_name)

    Results = namedtuple("Results", "secretion, TM, localisation")
    deeploc_dict = {}
    with open(input_name) as input_file:
        for line in input_file:
            if line.strip() == "":
                pass
            elif line.startswith("#"):
                pass
            else:
                sp = line.strip().split()
                if sp[1] == "Extracellular":
                    is_secreted = "Secreted"
                else:
                    is_secreted = "Not_secreted"

                if sp[1] == "Cell_membrane":
                    is_TM = "TM"
                else:
                    is_TM = "Not_TM"

                result = Results(is_secreted, is_TM, sp[1].strip())
                deeploc_dict[sp[0]] = result
    return deeploc_dict


def parsing_tmhmm(input_suffix):
    """Takes the SHORT output file of tmhmm. And returns a dictionary  with a named tuple as value and the genes as keys. """
    input_name = input_suffix + ".tmhmm.out"
    check_file_opens(input_name)

    Results = namedtuple("Results", "TM, TM_domains")
    tmhmm_dict = {}
    with open(input_name) as input_file:
        for line in input_file:
            if line.strip() == "":
                pass
            elif line.startswith("#"):
                pass
            else:
                sp = line.strip().split()
                nb_TM = int(sp[4].split("=")[1])
                if nb_TM == 0:
                    is_TM = "Not_TM"
                else:
                    is_TM = "TM"
                result = Results(is_TM, str(nb_TM))
                tmhmm_dict[sp[0]] = result
    return tmhmm_dict


def parsing_effectorP_apoplastP(input_name):
    """ Common function to get the results from EffectorP and ApoplastP. Returns a simple dictionary with the name of the protein as key and the result as value. """

    check_file_opens(input_name)

    P_dict = {}
    with open(input_name) as input_file:
        for line in input_file:
            if line.strip() == "":
                pass
            elif line.startswith("#"):
                pass
            elif line.strip() == "-----------------":
                pass
            else:
                sp = filter(None, line.strip().split())
                if sp[1] == "Unlikely":
                    P_dict[sp[0]] = "Unlikely effector"
                else:
                    P_dict[sp[0]] = sp[1]
    return P_dict


def parsing_effectorP(input_prefix):
    input_name = input_prefix + ".EffectorP.out"
    return parsing_effectorP_apoplastP(input_name)


def parsing_apoplastP(input_prefix):
    input_name = input_prefix + ".ApoplastP.out"
    return parsing_effectorP_apoplastP(input_name)


def parsing_LOCALIZER(input_prefix):
    input_name = input_prefix + ".LOCALIZER.out/Results.txt"
    check_file_opens(input_name)

    Loc_dict = {}
    with open(input_name) as input_file:
        for line in input_file:
            if line.strip() == "":
                pass
            elif line.startswith("#"):
                pass
            elif line.startswith("Identifier"):
                pass
            else:
                sp = line.strip().split("\t")
                where = []
                for localization, prediction in zip(["Chloroplast", "Mitochondrion", "Nucleus"], sp[1:]):
                    if prediction == '-':
                        pass
                    else:
                        if prediction.split()[0] == "Y":
                            where.append(localization)

                if len(where) == 0:
                    Loc_dict[sp[0].strip()] = "-"
                else:
                    Loc_dict[sp[0].strip()] = "/".join(where)
    return Loc_dict


def parsing_antismash(directory_path):
    """Parses the simplest output of antismash (v4.2.0). And returns a dictionary  with the gene names as keys and, as value, the cluster number (as numbered in the output I have seen of antismash v.4.2.0) and the putative function of the cluster. """
    input_name = directory_path + "/geneclusters.txt"
    check_file_opens(input_name)
    antism_dict = {}
    nb_genes = 0
    with open(input_name) as input_file:
        for cluster_nb, line in enumerate(input_file):
            if line.strip() == "":
                pass
            else:
                sp = line.strip().split()
                genes = sp[3].split(";")
                function = sp[2]
                for gene in genes:
                    nb_genes += 1
                    antism_dict[gene] = ";".join(
                        ['cluster' + str(cluster_nb).zfill(3), function])
    # Outputing the results from the effector prediction step
    print "\n\n"
    print "  ******************************"
    print "   Detection of MGC is complete"
    print "  ******************************"
    print "We have found", str(cluster_nb + 1), "metabolite gene clusters containing", str(nb_genes),  "genes."
    return antism_dict
