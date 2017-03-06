from blast import Blast
import os
from mutation import probable_mutations
from generate_blueprint import generate_blueprint
from remodel import remodel, relax
from score import generate_dictionary
from score import calculate_score


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

    # Load the dictionary.
    stability_dict = generate_dictionary(os.path.join("data", "stabilityScoreFile.txt"))

    (mutated, score) = probable_mutations(original_string, pmf, stability_dict, 2, method='min')
    print mutated
    print score


    # The original score
    original_score = calculate_score(original_string, stability_dict)
    print original_score

def extract_fasta_header(fasta_file_path):
    return open(fasta_file_path).read().split("\n")[0]

def save_in_fasta(file_path, header, sequence):
    mutated_fasta_file = open(file_path, "w")
    mutated_fasta_file.write(header + "\n" + sequence)
    mutated_fasta_file.close()

def produce(target, dict_path, num_mutate=2, method='min', output_path='rosetta_out'):
    """
        1. Blast a sequence in fasta file
    """
    fasta_file = os.path.join("fasta", "{}.fasta".format(target))
    xml_file = os.path.join("xml", "blast_{}.xml".format(os.path.splitext(os.path.basename(fasta_file))[0]))
    blast = Blast()
    if not os.path.exists(xml_file):
        print "Blast using NCBI..."
        blast.blast(fasta_file_path=fasta_file, xml_file_path=xml_file)
    else:
        blast.load_xml(xml_file)
        blast.load_sequences(fasta_file)
    original_string = blast.sequences[0]
    """
        2. Generate a possible mutation dictionary based on blasted results
    """
    pmf = blast.build_positional_distributions(record_index=0)
    stability_dict = generate_dictionary(os.path.join("data", dict_path))

    """
        3. Apply mutation and store it in fasta file
    """
    print "Mutating a sequence..."
    (mutated, score) = probable_mutations(original_string, pmf, stability_dict, num_mutate, method=method)
    print "Score: {}".format(score)

    save_in_fasta(os.path.join(output_path, target, "{}_mutated.fasta".format(target)),
                  extract_fasta_header(fasta_file),
                  mutated)
    """
        4. Generate a blueprint file from origianl and mutated sequences
    """
    blueprint_path = os.path.join("blueprint", "{}.remodel".format(target))
    generate_blueprint(original_string, mutated, blueprint_path)

    """
        4.5 Relax original pdb file
    """
    rosetta_path = open("rosetta_path.config").read().strip('\n')
    pdb_path = os.path.join("pdb", "{}.pdb".format(target))
    print "Start relaxing..."
    relax(rosetta_path, pdb_path, outfolder=os.path.join(output_path, target))

    """
        5. Run rosetta remodel to generate a mutated pdb file
    """
    # Now relaxed pdb file generated is stored at root named with 0001
    pdb_path = os.path.join(output_path, target, "{}_0001.pdb".format(target))
    print "Start remodeling... (Will take about an hour or more)"
    remodel(rosetta_path, pdb_path, blueprint_path, outfolder=os.path.join(output_path, target))

    """
        5.5 Relax the remodeled pdb
    """
    pdb_path = os.path.join(output_path, target, "{}_0001_0001.pdb".format(target))
    print "Start relaxing..."
    relax(rosetta_path, pdb_path, outfolder=os.path.join(output_path, target))


if __name__ == '__main__':
    #target = "1c20"
    #target = "1a7s"
    #target = "1b9o"
    #target = "1bkr"
    #target = "1ctq"
    #target = "1jsf"
    target = "1mfm"
    produce(target=target, dict_path="stabilityScoreFile.txt")

