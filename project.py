from blast import Blast
import os


def test():
    """Testing methods
    """
    fasta_file = os.path.join("fasta", "DROME_HH_Q02936.fasta")
    xml_file = os.path.join("xml", "blast_{}.xml".format(os.path.splitext(os.path.basename(fasta_file))[0]))
    blast = Blast()

    if not os.path.exists(xml_file):
        blast.blast(fasta_file_path=fasta_file, xml_file_path=xml_file)
    else:
        blast.load_xml(xml_file)
        blast.load_sequences(fasta_file)

    print "# of records(protein): {}".format(blast.records_size)
    print "Original Sequence: {}".format(blast.sequences[0])
    # showing first alignment of first protein
    blast.pretty_print(record_index=0, alignment_index=0)
    print ""

    # Now construct a character-distribution for each position
    pmf = blast.build_positional_distributions(record_index=0)
    for i in range(10):
        print "at {}: {}".format(i, pmf.get_distribution(i))  # showing first 10 distributions


if __name__ == '__main__':
    test()
