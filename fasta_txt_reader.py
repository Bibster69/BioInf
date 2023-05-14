

def read(path):
    file_content = open(path, 'r')
    seq_dictionary = {}
    for line in file_content:
        if '>' in line:
            current_seq = line.replace('\n', '')
            seq_dictionary[current_seq] = ""
        else:
            seq_dictionary[current_seq] += line.replace('\n', '')

    return seq_dictionary