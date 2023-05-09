import DNA_helper_functions as DNA_help
import structures
import BioSequences



if __name__ == '__main__':
    test_seq = BioSequences.BioSequence(seq = 'ACTG', type = 'DNA', label = 'testLabel')
    print(test_seq)
    print(test_seq.type)
    test_seq.generate_random_sequence()
    print(test_seq)
    test_seq.trainscribe()
    print(test_seq)