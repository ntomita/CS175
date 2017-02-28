from FASTAReader import FASTAFile
import os


def generate_blueprint(original_fasta_seq, mutated_fasta_seq, output_name):
    # Shorten variable name
    o_seq = original_fasta_seq
    m_seq = mutated_fasta_seq

    out_file = open(output_name, 'w')
    out_seq = ""
    for i in range(len(o_seq)):
        if o_seq[i] == m_seq[i]:
            out_seq = out_seq + "{} {} .".format(i+1, o_seq[i])
        else:
            out_seq = out_seq + "{} {} . PIKAA {}".format(i+1, o_seq[i], m_seq[i])
        if i != len(o_seq)-1:
            out_seq = out_seq + "\n"
    out_file.write(out_seq)
    out_file.close()


def test():
    o_fasta = os.path.join("fasta", "1c20.fasta")
    m_fasta = os.path.join("fasta", "1c20_m.fasta")
    o_seq = FASTAFile(o_fasta)[0].seq
    m_seq = FASTAFile(m_fasta)[0].seq

    output_name = os.path.join("blueprint", "1c20.remodel")
    generate_blueprint(o_seq, m_seq, output_name)

if __name__ == '__main__':
    test()
