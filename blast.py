from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from pmf import PMF


class Blast():

    def __init__(self):
        self.records = None  # placeholder
        self.sequences = None  # placeholder

    @property
    def records_size(self):
        return len(self.records) if self.records is not None else 0

    def alignments_size(self, record_index):
        try:
            return len(self.records[record_index].alignments)
        except:
            return 0

    def get_record(self, record_index):
        return self.records[record_index]

    def get_alignment(self, record_index, alignment_index):
        return self.records[record_index].alignments[alignment_index]

    def blast(self, fasta_file_path, xml_file_path="", hitlist_size=50,
              program='blastp', database='nr'):
        """Blast a given protein in a fasta file and returns handler
        If file path is given, it saves data in xml format.
        Args:
            fasta_file_path (string): a path for a protein in fasta format
            xml_file_path (string): a path to save data retrieved from NCBI
            hitlist_size (int): maximum number of different alignment sequence
            program (string): program blastn, blastp, blastx, tblastn, or tblastx (lower case)
            database Which database to search against (e.g. "nr").
        Returns:
            (list of Blast record)
        """
        fasta_string = open(fasta_file_path).read()
        result_handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=fasta_string,
            hitlist_size=hitlist_size)
        if xml_file_path != "":
            save_file = open(xml_file_path, 'w')
            save_file.write(result_handle.read())
            save_file.close()
        self.records = list(self.load_xml(xml_file_path))
        self.load_sequences(fasta_file_path)
        return self.records

    def load_xml(self, xml_file_path):
        """Returns blast records stored in a xml file
        Args:
            xml_file_path (string) where a file containing a result of blast
        Returns:
            (list of Blast record)
        """
        result_handle = open(xml_file_path)
        blast_records = NCBIXML.parse(result_handle)
        self.records = list(blast_records)
        return self.records

    def load_sequences(self, fasta_file_path):
        self.sequences = [str(s.seq) for s in SeqIO.parse(fasta_file_path, "fasta")]
        return self.sequences

    def pretty_print(self, record_index=0, alignment_index=0,
                     alignment=None, max_length=75, threshold=0.04):
        """Show the contents of an alignment object in a nice fashion
        It access to a stored blast object using given indices. If alignment
            parameter is given, it prints the contents in the object instead.
            In this case, indices given are ignored

        Args:
            record_index (int)
            alignment_index (int)
            alignment (alignment of Blast record)
            max_length (int): maximum length of amino acid sequences to show
            threshold (float): show pair whose expect value is under threshold
        """
        if alignment is None:
            print self.records[record_index].__class__
            alignment = self.records[record_index].alignments[alignment_index]

        for hsp in alignment.hsps:  # HPS: high scoring pair
            if hsp.expect < threshold:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print(hsp.query[0:min(max_length, len(hsp.query))] + '...')
                print(hsp.match[0:min(max_length, len(hsp.match))] + '...')
                print(hsp.sbjct[0:min(max_length, len(hsp.sbjct))] + '...')

    def build_positional_distributions(self, record_index=0, blast_record=None,
                                       threshold=0.04, min_length_ratio=0.8):
        """Count possible mutation at each position and returns distributions
        Construct pmf distributions of a record at given index. If blast_record
            object is given, construct distributions based on the object instead.

        Args:
            record_index (int)
            blast_record (Bio.Blast.Record.Blast)
            threshold (float): used to compare with hsp.expect
                (but not sure how this can be effective)
            min_length_ratio (float): value between 0 and 1, defining
                the minimum length compared to the length of original string
                match results should have. (we don't need a tiny sequence)
        Returns:
            (PMF) containing a list of dictionaries that represents distributions
                of symbols at each position of protein sequence
        """
        if blast_record is None:
            blast_record = self.records[record_index]

        query_length = blast_record.query_length
        pmf = PMF(query_length)

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect > threshold:
                    continue
                if (hsp.query_end - hsp.query_start + 0.0) / query_length < min_length_ratio:
                    continue
                pos = hsp.query_start - 1

                for i in range(len(hsp.query)):
                    if hsp.query[i] == '-':
                        # gap in query: ignore
                        continue
                    elif hsp.query[i] == hsp.match[i]:
                        # identity
                        pmf.increment(pos, hsp.sbjct[i])
                        pos += 1
                        continue
                    elif hsp.sbjct[i] == '-':
                        # gap in sbjct: ignore
                        pos += 1
                        continue
                    else:
                        # different symbol
                        pmf.increment(pos, hsp.sbjct[i])
                        pos += 1
        self.pmf = pmf
        return pmf
