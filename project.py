from blast import Blast
import os


def test():
    """Testing methods
    """
    fasta_file = os.path.join("fasta", "DROME_HH_Q02936.fasta")
    xml_file = os.path.join("xml", "blast_result.xml")
    blast = Blast()

    if not os.path.exists(xml_file):
        blast.blast(fasta_file_path=fasta_file, xml_file_path=xml_file)
    blast.open_xml(xml_file)

    print "# of records(protein): {}".format(blast.records_size)

    # showing first alignment of first protein
    blast.pretty_print(record_index=0, alignment_index=0)
    print ""

    # Now construct a character-distribution for each position
    pmf = blast.build_positional_distributions(record_index=0)
    for i in range(10):
        print "at {}: {}".format(i, pmf.get_distribution(i))  # showing first 10 distributions


if __name__ == '__main__':
    test()
