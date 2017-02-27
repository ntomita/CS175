from blast import Blast
import os
from mutation import probable_mutations

flexibility_dict = {'A':-0.605,'C':-0.693,'D':-0.279,'E':-0.160,'F':-0.719
                      ,'G':-0.537,'H':-.662,'I':-0.682,'K':-0.043,'L':-0.631
                      ,'M':-0.626,'N':-0.381,'P':-0.271,'Q':-0.369,'R':-0.448
                      ,'S':-0.423,'T':-0.525,'V':-0.669,'W':-0.727,'Y':-0.721}

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
    original_string = blast.sequences[0]
    print "our original: ",original_string
    # showing first alignment of first protein
    blast.pretty_print(record_index=0, alignment_index=0)
    print ""

    # Now construct a character-distribution for each position
    pmf = blast.build_positional_distributions(record_index=0)
    for i in range(10):
        print "at {}: {}".format(i, pmf.get_distribution(i))  # showing first 10 distributions

    (mutated, score) = probable_mutations(original_string, pmf, flexibility_dict, 4, method='min')
    print mutated


if __name__ == '__main__':
    test()
