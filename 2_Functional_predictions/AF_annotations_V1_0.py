

#    ******************************************************
# |    This class and associated methods were               |
# |           written by Alice Feurtey in fall of 2017      |
# |        It was last modified on the 18/04/18.            |
#    ******************************************************
import sys

accepted_format_list = ["GFF3", "GFF", "GTF", "BED"]


class GFF_feature:
    """ GFF_feature class represents a gff file line."""

    def __init__(self, ligne):
        """ Creates a GFF_feature with seqid, source, type, start, end, score, strand, phase, attributes and optionally supplementary from one ligne in a GFF format."""
        split1 = ligne.strip().split("\t")
        self.seqid = split1[0]
        self.source = split1[1]
        self.type = split1[2]
        self.start = split1[3]
        self.end = split1[4]
        self.score = split1[5]
        self.strand = split1[6]
        self.phase = split1[7]
        self.attributes = {}
        split2 = filter(None, split1[8].strip().split(";"))
        for i in split2:
            split3 = i.strip().split("=")
            if len(split3) == 1:
                self.attributes["Name"] = split3[0]
            else:
                self.attributes[split3[0]] = split3[1]
        if len(split1) > 9:
            self.supplementary = "\t".join(split1[9:])

    @classmethod
    def from_GTF_line(cls, ligne):
        """ Creates a GFF_feature with seqid, source, type, start, end, score, strand, phase, attributes and optionally supplementary from one ligne in a GTF format."""
        split1 = ligne.strip().split("\t")
        split1[8] = split1[8].replace(" ", "=")
        new_ligne = "\t".join(split1)
        return cls(new_ligne)

    @classmethod
    def from_BED_line(cls, ligne):
        """ Creates a GFF_feature with seqid, source, type, start, end, score, strand, phase, attributes and optionally supplementary from one ligne in a BED format. The bed format being on base 0, I add 1 to the start coordinate when converting to GFF_feature."""
        split1 = ligne.strip().split("\t")
        try:
            name = "Name=" + split1[3]
        except:
            name = "Name=" + split1[0] + "_" + split1[1]

        new_ligne = "\t".join([split1[0], "BED_file", "unk", str(int(split1[1]) + 1),
                               split1[2], ".", "+", ".", name])
        return cls(new_ligne)

    @classmethod
    def from_scratch(cls, chromosome, start, end, name, strand="+"):
        """ From the chromosome, start, end and name, creates a GFF_feature."""
        ligne = "\t".join([chromosome, "from_scratch", "unk",
                           start, end, ".", strand, ".", "Name=" + name])
        return cls(ligne)

    @classmethod
    def from_biopython_maf_record(cls, biopython_maf_record):
        """ From the a record based on the BioPython creates a GFF_feature."""
        chromosome = biopython_maf_record.id.split(".")[1]
        sample = biopython_maf_record.id.split(".")[0]
        strand = biopython_maf_record.annotations["strand"]
        if strand == 1:
            start = int(biopython_maf_record.annotations["start"]) + 1
            end = start + int(biopython_maf_record.annotations['size'])
            name = "_".join([sample, chromosome, str(start)])
            ligne = "\t".join([chromosome, "from_scratch", "unk",
                               str(start), str(end), ".", "+", ".", "Name=" + name])
            return cls(ligne)
        else:
            print "Alice was too lazy to implement this when the strand is minus because of possible +1 or -1."
            print biopython_maf_record.id, biopython_maf_record.annotations["strand"]

    def write_BED_line(self):
        """ Returns a string in GTF format from a GFF feature. Untested as yet. The bed format being on base 0, I remove 1 to the GFF feature start coordinate."""
        attributes_as_str = ";".join(
            [self.attributes[clef] for clef in sorted(self.attributes.keys())])
        written_GFF_line = "\t".join([self.seqid, str(int(self.start) - 1),
                                      self.end, attributes_as_str])
        return written_GFF_line

    def write_GFF_line(self):
        """ Returns a string in GFF format from a GFF feature. """
        attributes_as_str = ";".join(
            [str(clef + "=" + self.attributes[clef]) for clef in sorted(self.attributes.keys())])
        written_GFF_line = "\t".join([self.seqid, self.source, self.type, self.start,
                                      self.end, self.score, self.strand, self.phase, attributes_as_str])
        return written_GFF_line

    def write_GTF_line(self):
        """ Returns a string in GTF format from a GFF feature. Untested as yet."""
        attributes_as_str = ";".join(
            [clef + " " + self.attributes[clef] for clef in sorted(self.attributes.keys())])
        written_GFF_line = "\t".join([self.seqid, self.source, self.type, self.start,
                                      self.end, self.score, self.strand, self.phase, attributes_as_str])
        return written_GFF_line


def annot_parse(annotation_input_name, annotation_input_format="GFF"):
    """ Takes an annotation file name and its format which can be GFF, GTF or BED (with GFF as the default) and returns a tuple with a list of header lines and a list of features of classe GFF_feature """
    import sys
    headers = []
    features = []
    with open(annotation_input_name, "r") as infile:
        for ligne in infile:
            if ligne.startswith("#"):
                headers.append(ligne.strip())
            elif ligne.strip() == "":
                pass
            else:
                if annotation_input_format.upper() == "GFF" or annotation_input_format.upper() == "GFF3":
                    temp = GFF_feature(ligne)
                elif annotation_input_format.upper() == "GTF":
                    temp = GFF_feature.from_GTF_line(ligne)
                elif annotation_input_format.upper() == "BED":
                    temp = GFF_feature.from_BED_line(ligne)
                else:
                    sys.exit(
                        "As of now, this function only supports BED, GFF and GTF.")

                features.append(temp)
    return (headers, features)


def annot_parse_and_sort_per_chr(annotation_input_name, annotation_input_format="GFF"):
    """ Takes an annotation file name and its format (either GFF, the default, GTF or BED) and returns a tuple with a list of header lines and a dictionary (defaultdict) of features of classe GFF_feature """
    from collections import defaultdict
    import sys
    headers = []
    features = defaultdict(list)
    with open(annotation_input_name, "r") as infile:
        for ligne in infile:
            if ligne.startswith("#"):
                headers.append(ligne.strip())
            elif ligne.strip() == "":
                pass
            else:
                if annotation_input_format.upper() == "GFF" or annotation_input_format.upper() == "GFF3":
                    temp = GFF_feature(ligne)
                elif annotation_input_format.upper() == "GTF":
                    temp = GFF_feature.from_GTF_line(ligne)
                elif annotation_input_format.upper() == "BED":
                    temp = GFF_feature.from_BED_line(ligne)
                else:
                    sys.exit(
                        "As of now, this function only supports BED, GFF and GTF.")
                features[temp.seqid].append(temp)
    return (headers, features)


def other_parse(file_name, file_format):
    headers = []
    features = []
    with open(file_name, "r") as infile:
        for ligne in infile:
            if ligne.startswith("#"):
                headers.append(ligne.strip())
            else:
                ligne = ligne.strip().split("\t")
                temp = GFF_feature.from_scratch(ligne[file_format - 1])
                features.append(temp)
    return (headers, features)


def write_in_annot_format(gff_feature, writing_format):
    """ Takes a GFF_feature object and a writing format. Then calls the functions for  """
    writing_format = writing_format.upper()
    if writing_format not in accepted_format_list:
        print "Your writing format is incorrect since I only support GFF, GTF and BED at present. I will use bed by default."
        writing_format = "BED"
    if writing_format == "GFF" or writing_format == "GFF3":
        return gff_feature.write_GFF_line()
    elif writing_format == "GTF":
        return gff_feature.write_GTF_line()
    else:
        return gff_feature.write_BED_line()


def check_file_format(annotation_input_name, annotation_input_format=None):
    """This function takes a input file and, optionally, the format of this input file as bed, gtf or gff. It then checks that the format name (ie it does not check that the format in the file itself is correct) is indeed one of the accepted one. If not, it will look for the format as the extension of the file and check again against the list of accepted formats. It will return the file name and the file format."""
    file_format_found = False
    accepted_formats = accepted_format_list
    while not file_format_found:
        if annotation_input_format:
            annotation_input_format = annotation_input_format.upper()
            if annotation_input_format not in accepted_formats:

                sys.exit("ERROR. You have specified the format of your annotation file named " + annotation_input_file + " as" + annotation_input_format +
                         "which is not one of the accepted format. Check that your format is indeed one of the accepted one (listed in the help). I will try to look at the extension instead.")
                annotation_input_format = annotation_input_file.rsplit(
                    ".")[-1].upper()
            else:
                file_format_found = True

        else:
            annotation_input_format = annotation_input_name.rsplit(
                ".")[-1].upper()
            if annotation_input_format not in accepted_formats:
                sys.exit("ERROR. You have not specified the format of your annotation file " + annotation_input_name + ". Your file name ends with" + annotation_input_format +
                         "which is not one of the accepted format. Check that your format is indeed one of the accepted ones: " + " ".join(accepted_formats) + ". If it is, you need to either change your file name or specify the format explicitly.")
            else:
                file_format_found = True

    return annotation_input_name, annotation_input_format


def reverse_annotations(file_name, file_format, output_file_name="based_on_input", output_file_format="based_on_input"):
    """ This script transforms a gff annotation into its reverse annotations: the output annotations are  the 'empty' spaces between the input annotations. WARNING. Because there is no way to know where the chromosomes end, the last annotation has 'INF' instead of an end coordinate. To get rid of them, just do grep -v 'INF' or grep -v '_and_end'. """
    import sys
    if output_file_format == "based_on_input":
        output_file_format = file_format
    if output_file_name == "based_on_input":
        output_file_name = file_name.rsplit(
            ".", 1)[0] + ".reverse." + output_file_format.lower()

    remembered_chrom = "Whatever"
    remembered_end = None
    out = open(output_file_name, "w")
    with open(file_name, "r") as infile:
        for indice, ligne in enumerate(infile):
            if ligne.startswith("#"):
                out.write(ligne)
            elif ligne.strip() == "":
                pass
            else:
                if file_format.upper() == "GFF" or file_format.upper() == "GFF3":
                    temp = GFF_feature(ligne)
                elif file_format.upper() == "GTF":
                    temp = GFF_feature.from_GTF_line(ligne)
                elif file_format.upper() == "BED":
                    temp = GFF_feature.from_BED_line(ligne)
                else:
                    sys.exit(
                        "As of now, this function only supports BED, GFF and GTF.")
                if "Name" in temp.attributes:
                    name = temp.attributes["Name"]
                elif "ID"in temp.attributes:
                    name = temp.attributes["ID"]
                else:
                    name = str("annotation" + str(indice))

                if remembered_chrom == temp.seqid:
                    annot_name = "Between_" + name + "_and_" + remembered_name
                    to_write = GFF_feature.from_scratch(temp.seqid, remembered_end, str(
                        int(temp.start) - 1), annot_name, strand=".")
                else:
                    if remembered_end:
                        annot_name = "Between_" + remembered_name + "_and_" + remembered_chrom + "_end"
                        to_write = GFF_feature.from_scratch(
                            remembered_chrom, remembered_end, "INF", annot_name, strand=".")
                        out.write(write_in_annot_format(
                            to_write, output_file_format) + "\n")
                    annot_name = "Between_" + temp.seqid + "_start_and_" + name
                    to_write = GFF_feature.from_scratch(temp.seqid, "1", str(
                        int(temp.start) - 1), annot_name, strand=".")

                remembered_end = str(int(temp.end) + 1)
                remembered_chrom = temp.seqid
                remembered_name = name
                out.write(write_in_annot_format(
                    to_write, output_file_format) + "\n")

        annot_name = "Between_" + remembered_name + "_and_" + remembered_chrom + "_end"
        to_write = GFF_feature.from_scratch(
            remembered_chrom, remembered_end, "INF", annot_name, strand=".")
        out.write(write_in_annot_format(to_write, output_file_format) + "\n")
    out.close()


def compare_2_annotations(gff_feature1, gff_feature2, full_comparison=True):
    """ Takes two GFF_features. Does not check for strands.
    If full_comparison, returns a sentence comparing the annotations. Else, returns True if overlap and False if no overlap. """
    if gff_feature1.seqid != gff_feature2.seqid:
        if full_comparison:
            return "on different chromosomes"
        else:
            return False
    else:
        start1 = int(gff_feature1.start)
        end1 = int(gff_feature1.end)
        start2 = int(gff_feature2.start)
        end2 = int(gff_feature2.end)
        if start1 <= start2:
            if end1 <= start2:
                if full_comparison:
                    return "Annotation 2 is after annotation 1"
                else:
                    return False
            elif end2 <= end1:
                if full_comparison:
                    return "Annotation 2 is included in annotation 1"
                else:
                    return True
            else:
                if full_comparison:
                    return "Annotation 2 begins in annotation 1"
                else:
                    return True
        else:
            if end2 <= start1:
                if full_comparison:
                    return "Annotation 2 is before annotation 1"
                else:
                    return False
            elif end2 >= end1:
                if full_comparison:
                    return "Annotation 2 includes annotation 1"
                else:
                    return True
            else:
                if full_comparison:
                    return "Annotation 2 ends in annotation 1"
                else:
                    return True
