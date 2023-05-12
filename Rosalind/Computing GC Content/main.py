import DNA_helper_functions as DNA



if __name__ == '__main__':
    file_content = open("rosalind_gc.txt", 'r')
    seq_dictionary = {}
    for line in file_content:
        if '>' in line:
            current_seq = line.replace('\n', '')
            seq_dictionary[current_seq] = ""
        else:
            seq_dictionary[current_seq] += line.replace('\n', '')


    largest_seq_name = ''
    largest_seq_val = 0.0
    for k, v in seq_dictionary.items():
        gc_content = DNA.get_gc_content(v)
        if gc_content > largest_seq_val:
            largest_seq_name = k
            largest_seq_val = gc_content

    print(largest_seq_name[1:] + '\n' + str(largest_seq_val))