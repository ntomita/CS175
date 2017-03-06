from blast import Blast
import os
import operator

# A function to generate the score dictionary.
def generate_dictionary(score_file):
    if not os.path.exists(score_file):
        print "Please provide an existing score file in the format of character:score per line."
        return

    f = open(score_file, "r")
    score_dictionary = {}
    for line in f:
        if line.strip():
            character = line.split(':')[0]
            score = line.split(':')[1]
        #print score
            score_dictionary[character] = float(score)
    f.close()

    return score_dictionary

 # A generic function to calculate score for a given sequence and dictionary based on a feature.
def calculate_score(sequence, dictionary):
    # Initialize the score to be 0
    score = 0

    # Sum the score based on the dictionary
    for c in sequence:
        if c in dictionary:
            score += dictionary.get(c)
    return score

def test():
    fasta_file = os.path.join("fasta", "DROME_HH_Q02936.fasta")
    xml_file = os.path.join("xml", "blast_{}.xml".format(os.path.splitext(os.path.basename(fasta_file))[0]))
    blast = Blast()

    if not os.path.exists(xml_file):
        blast.blast(fasta_file_path=fasta_file, xml_file_path=xml_file)
    else:
        blast.load_xml(xml_file)
        blast.load_sequences(fasta_file)

    print "# of records(protein): {}".format(blast.records_size)
    original_string = blast.sequences[0]
    print "our original: ",original_string

    stability_dictionary = generate_dictionary("./data/stabilityScoreFile.txt")
    calculate_score(original_string, stability_dictionary)


if __name__ == '__main__':
    test()
