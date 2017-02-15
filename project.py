from Bio.Blast import NCBIWWW, NCBIXML
import os
from collections import Counter


class PMF():
    """Probability Mass Function
    Having bins with the length of n
    each bin is a counter dictionary, representing the number of
    occurence of a symbol at a position

    EX:
        >>> p = PMF(3)
        >>> p.increment(0, 'A')
        >>> p.get_distribution(0)
        {'A': 1.0}
        >>> p.increment(0, 'B')
        >>> p.get_distribution(0)
        {'A': 0.5, 'B': 0.5}

    """
    def __init__(self, length):
        self.bins = [Counter() for i in range(length)]

    def get_distribution(self, pos):
        return {k: v*(1.0/sum(self.bins[pos].values())) for k, v in self.bins[pos].items()}

    def increment(self, pos, symbol):
        self.bins[pos][symbol] += 1


def blast(fasta_file_path, xml_file_path="", hitlist_size=50):
    """Blast a given protein in a fasta file and returns handler
    If file path is given, it saves data in xml format.
    Args:
        fasta_file_path (string): a path for a protein in fasta format
        xml_file_path (string): a path to save data retrieved from NCBI
        hitlist_size (int): maximum number of different alignment sequence
    Returns:
        (Blast records handler)
    """
    fasta_string = open(fasta_file_path).read()
    result_handle = NCBIWWW.qblast(
        program='blastp',
        database='nr',
        sequence=fasta_string,
        hitlist_size=hitlist_size)
    if xml_file_path != "":
        save_file = open(xml_file_path, 'w')
        save_file.write(result_handle.read())
        save_file.close()
    return result_handle


def open_xml(xml_file_path):
    """Returns blast records stored in a xml file
    Args:
        xml_file_path (string) where a file containing a result of blast
    Returns:
        (list of Blast record)
    """
    result_handle = open(xml_file_path)
    blast_records = NCBIXML.parse(result_handle)
    return blast_records


def pretty_print(alignment, max_length=75, threshold=0.04):
    """Show the contents of an alignment object in a nice fashion
    Args:
        alignment (alignment of Blast record)
        max_length (int): maximum length of amino acid sequences to show
    """
    for hsp in alignment.hsps:  # HPS: high scoring pair
        if hsp.expect < threshold:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query[0:min(max_length, len(hsp.query))] + '...')
            print(hsp.match[0:min(max_length, len(hsp.match))] + '...')
            print(hsp.sbjct[0:min(max_length, len(hsp.sbjct))] + '...')


def build_positional_distributions(blast_record, threshold=0.04):
    """Count possible mutation at each position and returns distributions
    Args:
        blast_record (Bio.Blast.Record.Blast)
    Returns:
        (PMF) containing a list of dictionaries that represents distributions
            of symbols at each position of protein sequence
    """
    query_length = blast_record.query_length
    pmf = PMF(query_length)

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect > threshold:
                continue
            pos = hsp.query_start - 1

            for i in range(len(hsp.query)):
                if hsp.query[i] == '-':
                    continue  # gap in query: ignore
                elif hsp.query[i] == hsp.match[i]:
                    pos += 1
                    continue  # identity: ignore
                elif hsp.sbjct[i] == '-':
                    pos += 1
                    continue  # gap in sbjct: ignore
                else:
                    pmf.increment(pos, hsp.sbjct[i])
                    pos += 1
    return pmf


def test():
    """Testing methods
    """
    fasta_file = os.path.join("fasta", "DROME_HH_Q02936.fasta")
    xml_file = os.path.join("xml", "blast_result.xml")
    if not os.path.exists(xml_file):
        blast(fasta_file_path=fasta_file, xml_file_path=xml_file)

    # extract only one because the fast file have only one protein
    blast_record = list(open_xml(xml_file))[0]
    for alignment in blast_record.alignments:
        pretty_print(alignment)
        break  # just showing one alignment as an example

    pmf = build_positional_distributions(blast_record)
    for i in range(10):
        print(pmf.get_distribution(i))  # showing first 10 distributions


if __name__ == '__main__':
    test()
