#!/usr/bin/env python
import argparse as ap
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def prepare_argparser():
        description = "Output Reverse Complement of Fasta."
        epilog = "For command line options of each command, type %(prog)s COMMAND -h"
        argparser = ap.ArgumentParser(description=description, epilog = epilog)
        argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "input FASTA file")
        argparser.add_argument("-o","--output",dest = "outfile", type = str,required = True, help = "output file, default is stdout")
        return(argparser)

def make_rc_record(record):
    """Returns a new SeqRecord with the reverse complement sequence."""
    return SeqRecord(seq = record.seq.reverse_complement(), \
                 id = "rc_" + record.id, \
                 description = "reverse complement")

def main():
	argparser = prepare_argparser()
	args = argparser.parse_args()
	
	try:
        	infile = open(args.infile,"rb")
        except IOError,message:
                print >> sys.stderr, "cannot open Input file",message
                sys.exit(1)

	records = map(make_rc_record, SeqIO.parse(args.infile, "fasta"))
	SeqIO.write(records, args.outfile, "fasta")


if __name__=="__main__":
        main()

