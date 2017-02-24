from blast import Blast
import os
import operator
def probable_mutations(original,pmf,n):
    """Based on hsp sequences and Location parameters values of Amino acids,
    calculate the probable mutations depending on a threshold value

    Args:
        original - the original sequence
        pmf - the PMF distribution
        n -  number of deltas given by the user

        Returns:
        mut_string - the mutated string
    """

    #print "input inside function: ", original
    #print "pmf inside fiunction: ",pmf.get_distribution(0)
    delta_tuples = []     
    location_param_dict = {'A':-0.605,'C':-0.693,'D':-0.279,'E':-0.160,'F':-0.719
                        ,'G':-0.537,'H':-.662,'I':-0.682,'K':-0.043,'L':-0.631
                        ,'M':-0.626,'N':-0.381,'P':-0.271,'Q':-0.369,'R':-0.448
                        ,'S':-0.423,'T':-0.525,'V':-0.669,'W':-0.727,'Y':-0.721} 
    # Step 1: traverse through the original string to carry out the mutations
    for itr in xrange(0, len(original)):
        #Step 2: for every character of the original string, look at the corresponding
        #dictionary in the PMF
        pmf_dict = pmf.get_distribution(itr)
        str_char = original[itr]
        curr_tuple = (0,0,0)
        # Step 3: traverse through all key value pairs in the pmf dictionary
        for key in pmf_dict:
            # Step 4: caclulate Delta of the scores of the key and string's character
            delta = abs(location_param_dict[key]-location_param_dict[str_char])
            if curr_tuple[1] <= delta:
                curr_tuple = (key, delta, itr)

        # Step 5: add the tuples into the corresponding index of the list of 
        #deltas for the original string
        delta_tuples.append(curr_tuple)


    #delta_tuples.sort(key=lambda x: x[1])
    #import pdb; pdb.set_trace()
    sorted_delta_tuples = sorted(delta_tuples, key=lambda x: x[1], reverse=True)
    
    """Mutate the original string by taking highest n values of
     deltas from sorted list"""
    list_orig = list(original)
    
    for i in xrange(n):
        curr_t = sorted_delta_tuples[i]
        #print original[curr_t[2]
        list_orig[curr_t[2]] = curr_t[0]

    mutated_string ="".join(list_orig)

    print original
    print mutated_string

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

    probable_mutations(original_string,pmf,4)

    

if __name__ == '__main__':
    test()
