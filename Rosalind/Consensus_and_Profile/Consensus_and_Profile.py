
import fasta_txt_reader

if __name__ == '__main__':
    seq_dict = fasta_txt_reader.read('test.txt')
    char_num = None
    for v in seq_dict.values():
        char_num = len(v)
        break

    A_counts = [0] * char_num
    C_counts = [0] * char_num
    G_counts = [0] * char_num
    T_counts = [0] * char_num


    for k, v in seq_dict.items():
        index = 0
        for char in list(v):
            if char == 'A':
                A_counts[index] += 1
            elif char == 'C':
                C_counts[index] += 1
            elif char == 'G':
                G_counts[index] += 1
            elif char == 'T':
                T_counts[index] += 1
            index += 1


    print(A_counts)
    print(C_counts)
    print(G_counts)
    print(T_counts)

    result = ''
    for i in range(len(A_counts)):
        counts = {'A': A_counts[i],
                  'C': C_counts[i],
                  'G': G_counts[i],
                  'T': T_counts[i]}

        result += max(counts, key=counts.get)

    print(result)
    print('A: ' + str(A_counts).replace(',', '').replace('[', ''))
    print('C: ' + str(C_counts).replace(',', '').replace('[', ''))
    print('G: ' + str(G_counts).replace(',', '').replace('[', ''))
    print('T: ' + str(T_counts).replace(',', '').replace('[', ''))






