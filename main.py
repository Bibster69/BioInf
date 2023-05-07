import DNA_helper_functions as DNA_help



if __name__ == '__main__':
    dummy_seq = DNA_help.random_DNA_seq(100)
    print(DNA_help.return_specified_codon_usage(dummy_seq, 'G'))