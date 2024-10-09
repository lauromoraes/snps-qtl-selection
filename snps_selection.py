# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Test2

import os
from importlib.metadata import metadata


# GT: PL:GQ: DP:ADP: BSDP:ACN
# GT: Genotype
# PL: Phred-scaled genotype likelihoods rounded to the closest integer
# GQ: Genotype quality
# DP: Read depth
# ADP: Counts for observed alleles, including the reference allele
# BSDP: Number of base calls (depth) for the 4 nucleotides in called SNVs sorted as A,C,G,T
# ACN: Predicted copy number of each allele taking into account the prediction of number of copies of the region surrounding the variant


def create_outputs_folder(outputs_folder_name='snps-selection-outputs'):
    """
    Create a new folder to store the outputs of the script in the current working directory (cwd)
    and return the path to the folder created.

    :param outputs_folder_name:
    :return:
    """

    # Set the outputs folder path
    outputs_folder_path = os.path.join(os.getcwd(), outputs_folder_name)

    # Create a new folder to store the outputs
    os.system(f'mkdir -p {outputs_folder_path}')

    # Return the outputs folder path
    return outputs_folder_path


def get_major_nuc(info: str):
    '''
    Get the major nucleotide from the info field of a VCF file.
    :param info:
    :return:
    '''
    nucs = ['A', 'C', 'G', 'T']  # Nucleotides
    counts = [int(x) for x in info.split(',')]  # Get the counts of each nucleotide
    i = counts.index(max(counts))  # Get the index of the major nucleotide
    return nucs[i], counts  # Return the major nucleotide and the counts of each nucleotide


def correct_genotype(nuc_ref: str, nuc_alt: str, sample_info: str):
    '''
    Correct the Genotype (GT) field of a sample info of a row of a VCF file.
    :param nuc_ref:
    :param nuc_alt:
    :param sample_info:
    :return:
    '''
    fields = sample_info.split(':')  # Split the sample info by ':'
    # print(fields)
    valids = ['0', '1', '.']  # Valid values for the Genotype (GT) field
    ganotypes = fields[0].split('/')  # Split the Genotype (GT) field by '/'
    flag = True  # Flag to check if the genotype was corrected

    # If the are more than one genotype value, check if the genotype values are valid and correct the genotype value,
    # according to the reference and alternative nucleotides
    if len(ganotypes) > 1:  # Check if there are more than one genotype value
        for genotype in ganotypes:  # Iterate over the genotype values
            if genotype not in valids:  # Check if the genotype value is not valid
                print('Not valid {}'.format(genotype))
                flag = False  # Set the flag to False
                break  # Break the loop

        # If the genotype values are valid, determine the nucleotide with the highest count and set the genotype
        # value to 0 or 1, according to the reference and alternative nucleotides
        if flag:  # Check if the flag is True
            nuc, counts = get_major_nuc(fields[4])  # Get the major nucleotide and the counts of each nucleotide
            # if fields[0] == './.':
            #     print('CATCH', nuc, nuc_ref, nuc_alt)
            if nuc == nuc_ref:  # Check if the major nucleotide is equal to the reference nucleotide
                fields[0] = '0'  # Set the genotype value to 0
            elif nuc == nuc_alt:  # Check if the major nucleotide is equal to the alternative nucleotide
                fields[0] = '1'  # Set the genotype value to 1
            else:  # If the major nucleotide is different from the reference and alternative nucleotides
                print('BUG', fields[0], nuc, nuc_ref, nuc_alt, counts)
                return False

    corrected_samples_info = ':'.join(fields)  # Join the corrected genotype values
    # print('corrected_info', corrected_info)
    return corrected_samples_info  # Return the samples info with the corrected genotype values


def correct_variants(file_path: str):
    '''
    Correct the variants of a VCF file. The corrected variants are saved in a new file. The variants that could not
    be corrected are saved in a bug file.
    :param file_path:
    :return:
    '''
    new_file_path = 'ployd_corrected_' + file_path  # Define the new file path
    bug_file_path = 'bug_' + file_path  # Define the file path to store the bug variants
    new_lines = list()  # List to store the new lines
    bug_lines = list()  # List to store the bug lines
    with open(file_path, 'r') as f:  # Open the file
        for idx, line in enumerate(f.readlines()):  # Iterate over the lines of the file
            # print(idx)
            if not line.startswith('#'):  # Check if the line is not a header line nor a comment line (starts with '#')
                # print()
                flag = True  # Flag to check if the variant could be corrected
                fields = line.split('\t')  # Split the line by '\t'
                nuc_ref = fields[3]  # Get the reference nucleotide
                nuc_alt = fields[4]  # Get the alternative nucleotide
                if len(nuc_ref) > 1 or len(nuc_alt) > 1:  # Check if the reference or alternative nucleotides have length greater than 1
                    new_lines.append(line)  # Append the line to the new lines
                    continue  # Continue to the next line
                for i in range(-4, -2):  # Iterate over the parental sample fields
                    # print(fields[i])
                    corrected = correct_genotype(nuc_ref, nuc_alt, fields[i])  # Correct the ployd field of the sample info
                    if corrected:  # Check if the ployd field was corrected
                        fields[i] = corrected  # Set the corrected ployd field
                    else:  # If the ployd field was not corrected
                        flag = False  # Set the flag to False
                        # break
                if flag:  # Check if the flag is True
                    new_line = '\t'.join([str(x) for x in fields]) + '\n'  # Join the fields and add a '\n' character
                    new_lines.append(new_line)  # Append the new line to the new lines
                else:  # If the flag is False
                    new_lines.append(line)  # Append the line to the new lines
                    bug_lines.append(line)  # Append the line to the bug lines
            else:  # If the line is a comment line
                new_lines.append(line)  # Append the line to the new lines
                bug_lines.append(line)  # Append the line to the bug lines
    with open(new_file_path, 'w') as f:  # Write the new lines in the new file
        f.write(''.join(new_lines))  # Join the new lines and write in the new file
    with open(bug_file_path, 'w') as f:  # Write the bug lines in the bug file
        f.write(''.join(bug_lines))  # Join the bug lines and write in the bug file
    return new_file_path  # Return the path to the new file


def filter_special_cases(file_path):
    '''
    Filter out the special cases of a VCF file. The special cases are not saved in a new file.
    A special case is a variant that has a reference or alternative nucleotide with length greater than 1.
    
    :param file_path:
    :return:
    '''

    new_file_path = 'without_specials_cases_' + file_path  # Define the new file path
    filtered_out_file_path = 'special_cases_' + file_path  # Define the file path to store the special cases
    new_filtered_lines = list()  # List to store the lines filtered out
    new_lines = list()  # List to store the new lines

    with open(file_path, 'r') as f:  # Open the file
        for line in f.readlines():  # Iterate over the lines of the file
            line = line.strip()  # Remove the '\n' character
            if len(line) == 0:  # Check if the line is empty
                continue  # Continue to the next line
            if not line.startswith('#'):  # Check if the line is not a header line nor a comment line (starts with '#')
                fields = line.split('\t')  # Split the line by '\t'
                try:
                    nuc_ref = fields[3]  # Get the reference nucleotide
                except:
                    print(fields)
                nuc_alt = fields[4]  # Get the alternative nucleotide
                if len(nuc_ref) == 1 and len(
                        nuc_alt) == 1:  # Check if the reference and alternative nucleotides have length 1
                    new_lines.append(line)  # Append the line to the new lines
                else:  # If the reference or alternative nucleotides have length greater than 1
                    new_filtered_lines.append(line)  # Append the line to the filtered out lines

            elif line.startswith('#CHROM'):  # Check if the line is the header of the file (starts with '#CHROM')
                fields = line[:-1].split('\t')  # Remove the '\n' character and split the line by '\t'
                sample_names = fields[-4:]  # Get the sample names from the header
                fields.extend([f'%{x}' for x in sample_names])  # Add the new fields to the header
                new_line = '\t'.join(fields) + '\n'  # Join the fields and add a '\n' character
                new_lines.append(new_line)  # Append the new line to the new lines

            else:  # If the line is a comment line
                new_lines.append(line)  # Append the line to the new lines

    with open(new_file_path, 'w') as f:  # Write the new lines in the new file
        f.write('\n'.join(new_lines))  # Join the new lines and write in the new file
    with open(filtered_out_file_path, 'w') as f:  # Write the filtered out lines in the filtered out file
        f.write('\n'.join(new_filtered_lines))  # Join the filtered out lines and write in the filtered out file

    return new_file_path  # Return the path to the new file


def is_valid(fields, eval_type=0, threshold=None, avg=None, std=None):  # TODO separar em funcoes cada eval_typeblu
    infos = fields[-4:]  # Get the parental_sup, parental_inf, pool_sup and pool_rnd fields of the info
    genotype = [x.split(':')[0] for x in infos]  # Get
    print(genotype)

    if genotype[0] == genotype[1]:
        print('DISCARD:  parental_sup_counts ({}) and parental_inf ({}) equal.'.format(genotype[0], genotype[1]))
        return False  # Return False

    # Get counts of each nucleotide from the parental_sup field of the info
    parental_sup_counts = [int(x) for x in fields[-4].split(':')[4].split(',')]  # Get the counts of each nucleotide
    parental_sup_total_reads = float(sum(parental_sup_counts))  # Get the total number of reads

    # Get counts of each nucleotide from the parental_inf field of the info
    parental_alt_counts = [int(x) for x in fields[-3].split(':')[4].split(',')]  # Get the counts of each nucleotide
    parental_alt_total_reads = float(sum(parental_alt_counts))  # Get the total number of reads

    # Get counts of each nucleotide from the pool_sup field of the info
    pool_sup_counts = [int(x) for x in fields[-2].split(':')[4].split(',')]  # Get the counts of each nucleotide
    pool_sup_total_reads = float(sum(pool_sup_counts))  # Get the total number of reads

    # Get counts of each nucleotide from the pool_rnd field of the info
    pool_rnd_counts = [int(x) for x in fields[-1].split(':')[4].split(',')]  # Get the counts of each nucleotide
    pool_rnf_total_reads = float(sum(pool_rnd_counts))  # Get the total number of reads

    # Check if the total number of reads is less than 3 for any sample field
    if parental_sup_total_reads <= 3. or parental_alt_total_reads <= 3. or pool_sup_total_reads <= 3. or pool_rnf_total_reads <= 3.:
        print('DESCARTAR: numero de reads insuficientes.', pool_sup_counts)
        return False

    def has_diff(counts):
        aux = sorted(counts, reverse=True)
        return True if abs(aux[0] - aux[1]) <= 1 else False

    if has_diff(parental_sup_counts) or has_diff(parental_alt_counts):
        print('DESCARTAR: diferenca no numero de reads insuficientes.', parental_sup_counts, parental_alt_counts)
        return False

    #     # ref_nuc = fields[3]
    # alelo_sup = fields[3] if genotype[0] == '0' else fields[4]
    # print('alelo_sup', genotype[0], alelo_sup)
    # nucs = ['A', 'C', 'G', 'T']
    #
    # i_ref = get_allele
    #
    #
    # # Necessita haver pelo menos uma read semelhante ao parental superior
    # if eval_type == 0:
    #     if pool_sup_counts[nucs.index(alelo_sup)] > 0:
    #         print('INCLUIR: existem {} reads semelhantes ao parental_sup_counts ({}).'.format(pool_sup_counts[i_ref],
    #                                                                                           alelo_sup),
    #               pool_sup_counts)
    #     else:
    #         print('DESCARTAR: existem {} reads semelhantes ao parental_sup_counts ({}).'.format(pool_sup_counts[i_ref],
    #                                                                                             alelo_sup),
    #               pool_sup_counts)
    #         return False
    #
    # # Percentual de reads iguais ao parental superior deve ser igual ou maior que o limite
    # elif eval_type == 1:
    #
    #     parental_sup_counts = [int(x) for x in fields[-4].split(':')[4].split(',')]
    #     parental_sup_total_reads = float(sum(pool_sup_counts))
    #     # percent_205 = float(parental_sup_counts[i_ref]) / parental_sup_total_reads
    #
    #     pool_sup_total_reads = float(sum(pool_sup_counts))
    #     pool_sup_percent = float(pool_sup_counts[i_ref]) / pool_sup_total_reads
    #     if pool_sup_percent >= threshold:
    #         print(
    #             'INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e superior.'.format(
    #                 pool_sup_percent, alelo_sup), pool_sup_counts)
    #     else:
    #         print(
    #             'DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e inferior.'.format(
    #                 pool_sup_percent, alelo_sup), pool_sup_counts)
    #         return False
    #
    # # Percentual de reads iguais ao parental superior deve ser o maior
    # elif eval_type == 2:
    #     percents = [float(x) / pool_sup_total_reads for x in pool_sup_counts]
    #     greather = max(percents)
    #     percent = percents[i_ref]
    #     if percent == greather:
    #         print('INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. E o maior.'.format(
    #             percent,
    #             alelo_sup),
    #             pool_sup_counts)
    #     else:
    #         print(
    #             'DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. Nao e o maior.'.format(
    #                 percent, alelo_sup), pool_sup_counts)
    #         return False
    #
    # # Percentual de reads iguais ao parental superior deve ser um dos maiores
    # elif eval_type == 3:
    #     percents = [float(x) / pool_sup_total_reads for x in pool_sup_counts]
    #     greather = max(percents)
    #     percent = percents[i_ref]
    #     diff = abs(percent - greather)
    #     if diff <= threshold and diff > 0.:
    #         print(
    #             'INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. Esta no limite.'.format(
    #                 percent, alelo_sup), pool_sup_counts)
    #     else:
    #         print(
    #             'DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. Esta fora do limite.'.format(
    #                 percent, alelo_sup), pool_sup_counts)
    #         return False
    #
    # elif eval_type == 4:
    #
    #     parental_sup_counts = [int(x) for x in fields[-4].split(':')[4].split(',')]
    #     parental_sup_total_reads = float(sum(pool_sup_counts))
    #
    #     pool_sup_total_reads = float(sum(pool_sup_counts))
    #
    #     pool_rnd_counts = [int(x) for x in fields[-1].split(':')[4].split(',')]
    #     pool_rnf_total_reads = float(sum(pool_rnd_counts))
    #
    #     if parental_sup_total_reads <= 0 or pool_sup_total_reads <= 0 or pool_rnf_total_reads <= 0:
    #         print('Insuficient number of reads.\n\t205: {} | 207: {} | 208: {}'.format(parental_sup_total_reads,
    #                                                                                    pool_sup_total_reads,
    #                                                                                    pool_rnf_total_reads))
    #         return False
    #
    #     # percent_205 = float(parental_sup_counts[i_ref]) / parental_sup_total_reads
    #     pool_sup_percent = float(pool_sup_counts[i_ref]) / pool_sup_total_reads
    #     percent_208 = float(pool_rnd_counts[i_ref]) / pool_rnf_total_reads
    #
    #     condition_207 = pool_sup_percent >= threshold
    #     condition_208 = percent_208 >= (avg - std) and percent_208 <= (avg + std)
    #
    #     if condition_207 and condition_208:
    #         print(
    #             'INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e superior.'.format(
    #                 pool_sup_percent, alelo_sup), pool_sup_counts)
    #     else:
    #         print(
    #             'DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e inferior.'.format(
    #                 pool_sup_percent, alelo_sup), pool_sup_counts)
    #         return False
    #
    # return True



def calc_percents(field):
    counts = [int(x) for x in field.split(':')[4].split(',')]
    num_reads = float(sum(counts))
    percents = [round(x / num_reads, 2) for x in counts] if num_reads > 0 else [0 for _ in counts]
    to_str = ','.join([str(x) for x in percents])
    return to_str


def put_percents(fields):
    new_fields = [calc_percents(x) for x in fields[-4:]]
    fields.extend(new_fields)
    new_line = '\t'.join(fields) + '\n'
    return new_line


def get_lines(file_path: str):
    '''
    Get the metadata lines and the variant lines of a VCF file. The metadata lines are the lines that start with '#'.
    :param file_path:
    :return:
    '''
    metadata_lines = list()  # List to store the metadata lines
    variant_lines = None  # List to store the variant lines

    # Get the metadata lines of the file and the variant lines
    with open(file_path, 'r') as f:
        variant_lines = [line for line in f.readlines()]  # Get the lines of the file
        # print(file_path, len(variant_lines))
        for i, line in enumerate(variant_lines):  # Iterate over the lines of the file
            if line.startswith('#'):  # Check if the line starts with '#'; if True, the line is a metadata line
                metadata_lines.append(line)  # Append the metadata line to the metadata lines
            else:  # If the line does not start with '#', the line is a variant line
                break  # Break the loop
        variant_lines = variant_lines[i:]  # Get the lines of the file without the metadata lines

    return metadata_lines, variant_lines  # Return the metadata lines and the variant lines


def get_info_fields(variant_line: str, samples_idx: list[int] = [-4, -3, -2, -1]) -> tuple[str, str, list[list[str]]]:
    '''
    Get the info fields of a variant line of a VCF file.
    :param variant_line:
    :return:
    '''
    fields = variant_line.split('\t')  # Split the line by '\t'

    allele_ref = fields[3]  # Get the reference allele
    allele_alt = fields[4]  # Get the alternative allele

    parental_sup_info = [x for x in fields[samples_idx[0]].split(':')] # Get the parental_sup field of the info
    parental_inf_info = [x for x in fields[samples_idx[1]].split(':')] # Get the parental_inf field of the info
    pool_sup_info = [x for x in fields[samples_idx[2]].split(':')] # Get the pool_sup field of the info
    pool_rnd_info = [x for x in fields[samples_idx[3]].split(':')] # Get the pool_rnd field of the info

    sample_info = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]
    return allele_ref, allele_alt, sample_info


def get_counts(sample_info: list[str]) -> list[int]:
    '''
    Get the counts of each nucleotide from the sample info fields of a row of a VCF file.
    :param sample_info:
    :return:
    '''
    counts = [int(x) for x in sample_info[4].split(',')]  # Get the counts of each nucleotide
    return counts  # Return the counts of each nucleotide


def valid_line(sample_info: list[str]) -> bool:
    '''
    Check if a line of a VCF file is valid.
    :param sample_info: 
    :return:
    '''
    # Get the sample info fields of the line
    parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info = sample_info
    print('parental_sup_info', parental_sup_info)
    print('parental_inf_info', parental_inf_info)
    print('pool_sup_info', pool_sup_info)
    print('pool_rnd_info', pool_rnd_info)
    

    # Genotype verification of the parental_sup and parental_inf fields of the info
    genotype_sup = [x.split(':')[0] for x in parental_sup_info]  # Get the genotype values of the parental_sup field of the info
    genotype_inf = [x.split(':')[0] for x in parental_inf_info]  # Get the genotype values of the parental_inf field of the info
    if genotype_sup == genotype_inf:  # Check if the genotype values are equal
        print(f'DISCARD:  parental_sup_counts ({genotype_sup}) and parental_inf ({genotype_inf}) equal.')
        return False  # Return False if the genotype values are equal

    # Minimum number of reads verification
    for sample in sample_info:  # Iterate over the sample info fields
        counts = get_counts(sample)  # Get the counts of each nucleotide
        if sum(counts) <= 3:  # Check if the total number of reads is less than 3
            print('DISCARD: insuficient number of reads.')
            return False  # Return False if the total number of reads is less than 3

    # Check if the counts of the parental_sup and parental_inf fields of the info have a difference of at most 1
    parental_sup_counts = get_counts(
        parental_sup_info)  # Get the counts of each nucleotide from the parental_sup field of the info
    parental_inf_counts = get_counts(
        parental_inf_info)  # Get the counts of each nucleotide from the parental_inf field of the info

    parental_sup_ordered = sorted(parental_sup_counts,
                                  reverse=True)  # Sort the counts of the parental_sup field of the info in descending order
    parental_inf_ordered = sorted(parental_inf_counts,
                                  reverse=True)  # Sort the counts of the parental_inf field of the info in descending order

    if abs(parental_sup_ordered[0] - parental_inf_ordered[0]) <= 1:
        print(
            f'DISCARD: difference in the number of reads insufficient. parental_sup_counts({parental_sup_counts}) and parental_inf_counts({parental_inf_counts})')
        return False

    return True  # Return True if the line is valid


def least_one(pool_sup_info, allele_nuc):
    '''
    Check if there is at least one read similar to the parental_sup allele in the pool_sup.
    :param pool_sup_info:
    :param allele_nuc_idx:
    :return:
    '''
    allele_nuc_idx = get_nuc_index(allele_nuc)
    pool_sup_counts = get_counts(pool_sup_info)
    if not pool_sup_counts[allele_nuc_idx] > 0:
        return True
    return False


def percent_threshold(parental_sup_info, pool_sup_info, allele_nuc_idx, threshold):
    parental_sup_counts = get_counts(parental_sup_info)
    allele_count_on_parental_sup = parental_sup_counts[allele_nuc_idx]

    pool_sup_counts = get_counts(pool_sup_info)
    pool_sup_percent = float(pool_sup_counts[allele_nuc_idx]) / float(allele_count_on_parental_sup)

    if pool_sup_percent < threshold:
        print(f'DISCARD: the percentage of reads similar of the allele {allele_nuc_idx} in the pool_sup ('
              f'{pool_sup_percent}) is less than {threshold}.')

        return True

    pool_sup_total_reads = float(sum(pool_sup_counts))
    pool_sup_percent = float(pool_sup_counts[i_ref]) / pool_sup_total_reads
    if pool_sup_percent >= threshold:
        print(
            'INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e superior.'.format(
                pool_sup_percent, alelo_sup), pool_sup_counts)
    else:
        print(
            'DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e inferior.'.format(
                pool_sup_percent, alelo_sup), pool_sup_counts)
        return False

def filter_cases(eval_type: int, ref: str, alt: str, sample_info: list[list[str]], **kwargs) -> bool:

    parental_sup_genotype = sample_info[0][0] # Get the genotype of the parental_sup field of the info
    allele_nuc = get_allele(parental_sup_genotype, ref, alt) # Get the allele according to the genotype

    # Check if there is at least one read similar to the parental_sup allele in the pool_sup
    if eval_type == 0:
        pool_sup_info = sample_info[2] # Get the pool_sup field of the info
        print(allele_nuc, pool_sup_info)
        return least_one(pool_sup_info, allele_nuc)
    #
    elif eval_type == 1:
        parental_sup_info = sample_info[0] # Get the parental_sup field of the info
        pool_sup_info = sample_info[2] # Get the pool_sup field of the info
        return percent_threshold(parental_sup_info, pool_sup_info, allele_nuc, kwargs['threshold'])
    else:
        return False


def set_output_file_name(file_path: str, eval_type: int, **kwargs) -> str:
    '''
    Set the name of the output file according to the evaluation type.
    :param file_path:
    :param eval_type:
    :param kwargs:
    :return:
    '''
    new_file_path = None

    if eval_type == 0:
        new_file_path = 'filtered_' + 'least_one_' + file_path
    elif eval_type == 1:
        new_file_path = 'filtered_' + 'sup_limit_{}_'.format((str(kwargs['threshold']))) + file_path
    elif eval_type == 2:
        new_file_path = 'filtered_' + 'greater_' + file_path
    elif eval_type == 3:
        new_file_path = 'filtered_' + 'relaxed_greater_{}_'.format((str(kwargs['threshold']))) + file_path
    elif eval_type == 4:
        new_file_path = 'filtered_' + 'sup_limit_{}_avg_{}_std_{}_'.format((str(kwargs['threshold'])),
                                                                           str(kwargs['avg']),
                                                                           str(kwargs['std'])) + file_path
    return new_file_path


def get_allele(genotype: str, ref: str, alt: str) -> str:
    '''
    Get the allele according to the genotype.
    :param genotype:
    :param ref:
    :param alt:
    :return:
    '''

    alelo_sup = ref if genotype == '0' else alt
    return alelo_sup

def get_nuc_index(nuc: str) -> int:
    '''
    Get the nucleotide index.
    :param nuc:
    :return:
    '''
    return ['A', 'C', 'G', 'T'].index(nuc)

def filter_by_compare(file_path, eval_type=0, **kwargs):
    new_file_path = set_output_file_name(file_path, eval_type, **kwargs)  # Set the name of the output file
    if not new_file_path:  # Check if is a valid output file name
        return

    metadata_lines, variant_lines = get_lines(file_path)  # Get the metadata lines and the variant lines of the file

    output_lines = metadata_lines  # List to store the output lines; start with the metadata lines

    for line in variant_lines:  # Iterate over the variant lines

        line = line.strip() # Remove the '\n' character
        if len(line) == 0: # Check if the line is empty
            continue # Continue to the next line

        allele_ref, allele_alt, sample_info = get_info_fields(line)  # Get the info fields of the line

        # Check if the reference and alternative alleles have length 1 and if the line is valid
        if len(allele_ref) == 1 and len(allele_alt) == 1 and valid_line(sample_info):

            filtered_out = filter_cases(eval_type=eval_type, ref=allele_ref, alt=allele_alt,
                                        sample_info=sample_info, **kwargs)  # Apply the filter to the line
            if not filtered_out:  # Check if the line was not filtered out
                output_lines.append(line)  # Append the line to the output lines



# ==================================================================

            if is_valid(fields, eval_type=eval_type, threshold=threshold, avg=avg, std=std):
                new_line = put_percents(fields)
                print('Valid: ', new_line)
                output_lines.append(new_line)
            else:
                print('Invalid: ', line)
                pass

    # Define new file name to outputs
    if eval_type == 0:
        new_file_path = 'filtered_' + 'least_one_' + file_path
    elif eval_type == 1:
        new_file_path = 'filtered_' + 'sup_limit_{}_'.format((str(threshold))) + file_path
    elif eval_type == 2:
        new_file_path = 'filtered_' + 'greater_' + file_path
    elif eval_type == 3:
        new_file_path = 'filtered_' + 'relaxed_greater_{}_'.format((str(threshold))) + file_path
    elif eval_type == 4:
        new_file_path = 'filtered_' + 'sup_limit_{}_avg_{}_std_{}_'.format((str(threshold)), str(avg),
                                                                           str(std)) + file_path

    new_lines = list()
    with open(file_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                fields = line.rstrip().split('\t')
                nuc_ref = fields[3]
                nuc_alt = fields[4]
                if len(nuc_ref) == 1 and len(nuc_alt) == 1:
                    if is_valid(fields, eval_type=eval_type, threshold=threshold, avg=avg, std=std):
                        new_line = put_percents(fields)
                        print('Valid: ', new_line)
                        new_lines.append(new_line)
                    else:
                        print('Invalid: ', line)
                        pass
            else:
                # Append headers
                new_lines.append(line)
    with open(new_file_path, 'w') as f:
        f.write(''.join(new_lines))


# # file_path = 'variantsfile.vcf'
# file_path = 'variantsfile_Annotated.vcf'
#
# new_file_path = correct_variants(file_path)
# new_file_path = filter_special_cases(new_file_path)
#
# # Somente uma ocorrencio do parental_sup
# filter_by_compare(new_file_path, eval_type=0)
#
# # Ocorrencia do pool_sup maior que limiar
# threshold = .75
# filter_by_compare(new_file_path, eval_type=1, threshold=threshold)
#
# # Ocorrencia do parental_sup deve ser a maior
# filter_by_compare(new_file_path, eval_type=2)
#
# # Ocorrencia do parental_sup deve ser uma das maiores
# threshold = .1
# filter_by_compare(new_file_path, eval_type=3, threshold=threshold)
#
# # Ocorrencia do pool_sup maior que limiar e do pool_randomico entre a os desvios
# threshold = .75
# filter_by_compare(new_file_path, eval_type=4, threshold=threshold, avg=.5, std=.02)

# Function to read and parse the YAML configuration file
def read_yaml_config(yaml_config_file_path: str) -> dict:
    '''
    Read and parse the YAML configuration file.
    :param yaml_config_file_path:
    :return:
    '''
    import yaml  # Import the YAML library
    with open(yaml_config_file_path, 'r') as f:  # Open the file
        config = yaml.safe_load(f)  # Load the file
    return config  # Return the configuration

# Define Main fucntion of the script
if __name__ == '__main__':
    # Set and create the outputs folder
    # outputs_folder = create_outputs_folder()

    # file_path = 'variantsfile.vcf'
    file_path = 'variantsfile_Annotated.vcf' # Set the path to the VCF file

    new_file_path = correct_variants(file_path)
    new_file_path = filter_special_cases(new_file_path)

    configs = read_yaml_config('config.yaml')
    for config in configs:
        print(config)

    # Somente uma ocorrencio do parental_sup
    # filter_by_compare(new_file_path, eval_type=0)
