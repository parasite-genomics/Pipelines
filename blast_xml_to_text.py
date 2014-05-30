__author__ = 'squiresrb'

def parse_blast_result(blast_file_path):
    from Bio.Blast import NCBIXML
    output_file = "%s.alignments.txt" % blast_file_path
    output_handle = open(output_file, 'w')
    result_handle = open(blast_file_path)
    blast_records = NCBIXML.parse(result_handle)
    #E_VALUE_THRESH = 0.04
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            output_handle.write("Query ID: %s\n\n" % blast_record.query)
            for hsp in alignment.hsps:
    #            if hsp.expect < E_VALUE_THRESH:
                 #if (hsp.query.find('-') != -1) or (hsp.sbjct.find('-') != -1):
                 output_handle.write('sequence: %s\n' % alignment.title)
                 output_handle.write('length: %s\n' % alignment.length)
                 output_handle.write('e value: %s\n' % hsp.expect)
                 output_handle.write(hsp.query + '\n')
                 output_handle.write(hsp.match + '\n')
                 output_handle.write(hsp.sbjct + '\n\n')

def main(args):
#    format_bed_file(args.bed_file)
    parse_blast_result(args.blast_file)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='blast_xml_to_text.py',
                                     usage="%(prog)s -in (BLAST XML file)\n",
                                     description='Reads BLAST XML file and saves to text file.',
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=15),
                                     add_help=False)
    required_group = parser.add_argument_group('Required')
    required_group.add_argument('-in', '--blast_file',  required=True, help='The BLAST file.') #type=argparse.FileType('r'),

    #options_group = parser.add_argument_group('Options')
    #options_group.add_argument('-sheet', '--excel_sheet', type=int, help='What sheet do you want to use as input?')
    #options_group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    #options_group.add_argument("-path", default=".", help=argparse.SUPPRESS)

    args = parser.parse_args()

    main(args)
