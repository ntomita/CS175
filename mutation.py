import operator
from score import calculate_score
import math


def weight(value, regularize=True):
    return math.log(math.exp(1) + value) if regularize else value


def probable_mutations(original, pmf, location_param_dict, n=2, method='min', regularize=True):
    """Based on hsp sequences and Location parameters values of Amino acids,
    calculate the probable mutations depending on a threshold value

    Args:
        original - the original sequence
        pmf - the PMF distribution
        n -  number of mutations given by the user
        method (string): optimization method:
            either minimization or maximization of objective function
            acceptable input: 'min', 'max'

    Returns:
        (string, score): the mutated string and corresponding score
    """

    delta_tuples = []
    # Step 1: traverse through the original string to carry out the mutations
    for itr in xrange(0, len(original)):
        # Step 2: for every character of the original string, look at the corresponding
        # dictionary in the PMF
        pmf_dict = pmf.get_distribution(itr)
        original_key = original[itr]
        if method == 'min':
            curr_tuple = (None, float('Inf'), None)
        else:
            curr_tuple = (None, -float('Inf'), None)
        # Step 3: traverse through all key value pairs in the pmf dictionary
        for key in pmf_dict:
            # Step 4: caclulate Delta of the scores of the key and string's character

            # Different cases to handle:
            #  a. the original character is in the location_param_dict, the mutated character is not in the location_param_dict
            #  b. the original character is not in the location_param_dict, the mutated character is in the location_param_dict
            #  c. both are not in the dictionary
            #  d. both are in the dictionary
            if original_key not in location_param_dict and key not in location_param_dict:
                continue
            original_score = location_param_dict[original_key] * weight(pmf_dict[original_key])
            if original_key in location_param_dict and key in location_param_dict:
                mutated_score = location_param_dict[key] * weight(pmf_dict[key])
                delta = mutated_score - original_score
                # Taking conservation of each residue into consideration
                # pmf_dict[key] serves as weight of each residue

                if method == 'min' and delta < 0 and curr_tuple[1] > delta:
                    curr_tuple = (key, delta, itr)
                elif method == 'max' and delta > 0 and curr_tuple[1] < delta:
                    curr_tuple = (key, delta, itr)

        # Step 5: add the tuples into the corresponding index of the list of
        # deltas for the original string
        delta_tuples.append(curr_tuple)
    if method == 'min':
        sorted_delta_tuples = sorted(delta_tuples, key=lambda x: x[1], reverse=False)
    else:
        sorted_delta_tuples = sorted(delta_tuples, key=lambda x: x[1], reverse=True)

    """Mutate the original string by taking highest n values of
     deltas from sorted list"""
    list_orig = list(original)

    for i in xrange(n):
        curr_t = sorted_delta_tuples[i]
        list_orig[curr_t[2]] = curr_t[0]

    print "List of Mutation Candidate:"
    for i in range(30):
        curr_t = sorted_delta_tuples[i]
        print "to {} at {} with positional delta score {}".format(curr_t[0], curr_t[2]+1, curr_t[1])

    mutated_string ="".join(list_orig)
    score = calculate_score(mutated_string, location_param_dict)
    return (mutated_string, score)
