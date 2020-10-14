#/usr/bin/python

from AF_annotations_V1_0 import *
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_rna
from collections import namedtuple
import argparse
import sys


def from_annot_to_prot_fasta(output, genomic_fasta, annotation_file, which_annot="CDS", which_gene_id="Parent", translated=True):
    """ This function takes a genome multi-fasta and an annotation file in GTF or GFF format. It then extracts the protein sequences."""
    print "\n... Wait. Creating a protein sequences from the input annotation and genome files."
    temp, annotation_format = check_file_format(annotation_file)
    if annotation_format == "BED":
        sys.exit("Only file in GTF and GFF format are allowed here.")
    headers, annotations = annot_parse(
        annotation_file, annotation_format)
    Short_annot = namedtuple("Short_annot", 'seqid, start, end, strand, name')

    # We start by reading the whole annotation file and gathering the exons for each gene.
    all_genes = []
    current_exons = []
    current_gene = ""
    for annot in annotations:
        if annot.type == which_annot:
            if annot.attributes[which_gene_id] != current_gene:
                all_genes.append(current_exons)
                current_exons = []
                current_gene = annot.attributes[which_gene_id]
            current_exons.append(Short_annot(annot.seqid, int(annot.start), annot.end, annot.strand,
                                             annot.attributes[which_gene_id]))
    all_genes.append(current_exons)

    # Now that we have the exons to extract, we get the sequences (either DNA or translated).
    genome = SeqIO.to_dict(SeqIO.parse(genomic_fasta, "fasta"))
    # print genome.keys()
    with open(output, "w") as out:
        for gene in all_genes:
            sequence = Seq("", generic_rna)
            if len(gene) > 0:
                for exon in sorted(gene):
                    sequence += genome[exon.seqid].seq[int(
                        exon.start) - 1: int(exon.end)]
                out.write(">" + gene[0].name + "\n")
                if gene[0].strand == "-":
                    sequence = sequence.reverse_complement()
                if translated:
                    out.write(str(sequence.translate()) + "\n")
                else:
                    out.write(str(sequence) + "\n")
    print "... Done\n"
