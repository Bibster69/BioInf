import random
import structures
import collections


class BioSequence:
    def __init__(self, seq = 'ATCG', type = 'DNA', label = 'None'):
        self.sequence = seq.upper()
        self.type = type
        self.label = label
        if not self.__validate_seq():
            raise Exception('Given sequence is not valid')


    def __validate_seq(self):
        return set(set(self.sequence)).issubset(set(structures.DNA_Nucelotides))

    def generate_random_sequence(self, len = 50, type = 'DNA'):
        self.sequence = ''
        self.type = type
        if type == 'DNA':
            for i in range(len):
                self.sequence += random.choice(structures.DNA_Nucelotides)
        elif type == 'RNA':
            for i in range(len):
                self.sequence += random.choice(structures.RNA_Nucleotides)

    def __str__(self):
        return 'BioSequence Info\nLabel: ' + self.label + '\nType: ' + self.type + '\nSequence: ' + self.sequence

    def count_nucleotide_frequency(self):
        return dict(collections.Counter(self.sequence))
        # return pd.Series(list(seq)).value_counts().to_dict()

    def trainscribe(self):
        if self.type == 'DNA':
            self.sequence = self.sequence.replace('T', 'U')
            self.type = 'RNA'
        elif self.type == 'RNA':
            self.sequence = self.sequence.replace('U', 'T')
            self.type = 'DNA'
        else:
            pass

    def reverse_complement(self):
        if self.type == 'RNA':
            # raise Exception('RNA has no reverse complement')
            return 'RNA has no reverse complement'
        map = str.maketrans('ATCG', 'TAGC')
        return self.sequence.translate(map)[::-1]

    def gc_content(self):
        return ((self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence)) * 100

    def get_subseq_gc_contents(self, subseq_len):
        result = []
        for i in range(0, len(self.sequence) - subseq_len + 1, subseq_len):
            subseq = self.sequence[i:i + subseq_len]
            result.append(((subseq.count('G') + subseq.count('C')) / len(subseq)) * 100)
        return result

    def translate_to_aminoacid(self, start_index = 0):
        return [structures.DNA_Codons[self.sequence[index:index + 3]] for index in range(start_index, len(self.sequence) - 2, 3)]

    def get_specified_codon_usage(self, aminoacid):
        occurance_list = []
        for i in range(0, len(self.sequence) - 2, 3):
            tmp = self.sequence[i:i + 3]
            if structures.DNA_Codons[self.sequence[i:i + 3]] == aminoacid:
                occurance_list.append(tmp)

        result = dict(collections.Counter(occurance_list))
        total_weight = sum(result.values())
        for seq in result:
            result[seq] = round(result[seq] / total_weight, 2)
        return result

    def get_reading_frames_from_seq(self):
        frames = []
        frames.append(self.translate_to_aminoacid(0))
        frames.append(self.translate_to_aminoacid(1))
        frames.append(self.translate_to_aminoacid(2))
        tmp = BioSequence(seq = self.reverse_complement(), type = self.type)
        frames.append(tmp.translate_to_aminoacid(0))
        frames.append(tmp.translate_to_aminoacid(1))
        frames.append(tmp.translate_to_aminoacid(2))
        del tmp
        return frames

    def get_proteins_from_ORF(self, start_index=0, end_index=0, ordered=False):

        reading_frames = None
        if end_index > start_index:
            reading_frames = self.get_reading_frames_from_seq()
        else:
            reading_frames = self.get_reading_frames_from_seq()

        result = []
        for frame in reading_frames:
            found_proteins = BioSequence.get_protein_from_frame(frame)
            for protein in found_proteins:
                result.append(protein)

        if ordered:
            return sorted(result, key=len, reverse=True)  # reverse = true -> od najwiekszych do najmniejszych

        return result


    @staticmethod
    def get_protein_from_frame(frame):
        found_protein = []
        all_found_proteins = []
        for aminoacid in frame:
            if aminoacid == '_':
                if found_protein:
                    for x in found_protein:
                        all_found_proteins.append(x)
                    found_protein = []
            else:
                if aminoacid == 'M':
                    found_protein.append("")  # żeby for zadzałał
                for i in range(len(found_protein)):  # jeśli nie znalazł M to nic nie dodaje
                    found_protein[i] += aminoacid
        del found_protein
        return all_found_proteins

    @staticmethod
    def return_random_sequence(len = 50):
        sequence = ''
        for i in range(len):
            sequence += random.choice(structures.DNA_Nucelotides)
        return sequence

    @staticmethod
    def return_nucletide_counts(seq):
        return dict(collections.Counter(seq))

    @staticmethod
    def return_transctibtion(sequence, type = 'DNA'):
        seq = ''
        if type == 'DNA':
            seq = sequence.replace('T', 'U')
        elif type == 'RNA':
            seq = sequence.replace('U', 'T')
        else:
            return seq