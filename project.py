from Bio.Blast import NCBIWWW, NCBIXML
import os


def blast(fasta_file_path, xml_file_path=""):
    """Blast a given protein in a fasta file and returns handler
    If file path is given, it saves data in xml format.
    Args:
        fasta_file_path (string): a path for a protein in fasta format
        xml_file_path (string): a path to save data retrieved from NCBI
    Returns:
        (Blast records handler)
    """
    fasta_string = open(fasta_file_path).read()
    result_handle = NCBIWWW.qblast(
        program='blastp',
        database='nr',
        sequence=fasta_string,
        hitlist_size=50)
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
    """Show the contents of alignment object in a nice fashion
    Args:
        alignment (alignment of Blast record)
        max_length (int): maximum length of amino acid sequences to show
    """
    for hsp in alignment.hsps:
        if hsp.expect < threshold:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query[0:min(max_length, len(hsp.query))] + '...')
            print(hsp.match[0:min(max_length, len(hsp.match))] + '...')
            print(hsp.sbjct[0:min(max_length, len(hsp.sbjct))] + '...')


def main():
    fasta_file = "DROME_HH_Q02936.fasta"
    xml_file = "blast_result.xml"
    if not os.path.exists(xml_file):
        blast(fasta_file_path=fasta_file, xml_file_path=xml_file)

    # extract only one because the fast file have only one protein
    blast_records = list(open_xml(xml_file))[0]
    for alignment in blast_records.alignments:
        pretty_print(alignment)



if __name__ == '__main__':
    main()