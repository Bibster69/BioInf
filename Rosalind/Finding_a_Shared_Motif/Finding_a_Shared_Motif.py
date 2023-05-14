import fasta_txt_reader


def find_largest_common_substring(path):
    sequences = sorted(list(fasta_txt_reader.read(path).values()), key = len)
    #print(sequences)
    shortest_seq = sorted(sequences, key = len)[0]
    sequences = sequences[1:]
    #print(shortest_seq)
    #print(sequences)
    shortest_seq_substrings = []
    for i in range(0, len(shortest_seq)):
        for k in range(0, len(shortest_seq) - i):
            shortest_seq_substrings.append(shortest_seq[k:i + k + 1])

    shortest_seq_substrings = sorted(list(set(shortest_seq_substrings)), key = len)
    largest_substring = ''
    for substr in shortest_seq_substrings:
        found_flag = False
        for seq in sequences:
            if substr in seq:
                found_flag = True
            else:
                found_flag = False
                break
        if found_flag:
            largest_substring = substr

    return largest_substring




if __name__ == '__main__':
    print(find_largest_common_substring('test.txt'))

    # largest_substring = ''
    # curr_substring = ''
    # for n in shortest_seq:
    #     found_flag = False
    #     curr_substring += n
    #     print(curr_substring)
    #     for seq in sequences:
    #         if curr_substring in seq:
    #             found_flag = True
    #         else:
    #             found_flag = False
    #             curr_substring = ''
    #             break
    #     if found_flag:
    #         largest_substring = curr_substring
