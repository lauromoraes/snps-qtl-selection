# -*- coding: utf-8 -*-
#!/usr/bin/env python

import os

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
    nucs = ['A', 'C', 'G', 'T'] # Nucleotides
    counts = [int(x) for x in info.split(',')] # Get the counts of each nucleotide
    i = counts.index(max(counts)) # Get the index of the major nucleotide
    return nucs[i], counts # Return the major nucleotide and the counts of each nucleotide


def correct_ployd(nuc_ref: str, nuc_alt: str, sample_info: str):
    '''
    Correct the ployd field of a sample info of a row of a VCF file.
    :param nuc_ref:
    :param nuc_alt:
    :param sample_info:
    :return:
    '''
    fields = sample_info.split(':') # Split the sample info by ':'
    # print(fields)
    valids = ['0', '1', '.'] # Valid values for the ployd field
    ployds = fields[0].split('/') # Get the ployd values
    flag = True # Flag to check if the ployd values are valid
    if len(ployds) > 1: # Check if there are more than one ployd value
        for ployd in ployds: # Iterate over the ployd values
            if ployd not in valids: # Check if the ployd value is valid
                print('Not valid {}'.format(ployd))
                flag = False # Set the flag to False
                break # Break the loop
        if flag: # Check if the flag is True
            nuc, counts = get_major_nuc(fields[4]) # Get the major nucleotide and the counts of each nucleotide
            # if fields[0] == './.':
            #     print('CATCH', nuc, nuc_ref, nuc_alt)
            if nuc == nuc_ref: # Check if the major nucleotide is equal to the reference nucleotide
                fields[0] = '0' # Set the ployd value to 0
            elif nuc == nuc_alt: # Check if the major nucleotide is equal to the alternative nucleotide
                fields[0] = '1' # Set the ployd value to 1
            else: # If the major nucleotide is different from the reference and alternative nucleotides
                print('BUG', fields[0], nuc, nuc_ref, nuc_alt, counts)
                return False

    corrected_samples_info = ':'.join(fields) # Join the corrected ployd values
    # print('corrected_info', corrected_info)
    return corrected_samples_info # Return the samples info with the corrected ployd values


def correct_variants(file_path: str):
    '''
    Correct the variants of a VCF file. The corrected variants are saved in a new file. The variants that could not
    be corrected are saved in a bug file.
    :param file_path:
    :return:
    '''
    new_file_path = 'ployd_corrected_' + file_path
    bug_file_path = 'bug_' + file_path
    new_lines = list()
    bug_lines = list()
    with open(file_path, 'r') as f:
        cnt = 0
        for line in f.readlines():
            # print(line)
            if not line.startswith('#'):
                # print()
                flag = True
                fields = line.split('\t')
                nuc_ref = fields[3]
                nuc_alt = fields[4]
                if len(nuc_ref) > 1 or len(nuc_alt) > 1:
                    new_lines.append(line)
                    continue
                for i in range(-4, -2):
                    # print(fields[i])
                    corrected = correct_ployd(nuc_ref, nuc_alt, fields[i])
                    if corrected:
                        fields[i] = corrected
                    else:
                        flag = False
                        # break
                if flag:
                    new_line = '\t'.join([str(x) for x in fields])
                    new_lines.append(new_line)
                else:
                    new_lines.append(line)
                    bug_lines.append(line)
                # print(line)
                # print(new_line)
                # print()
                # cnt += 1
                # if cnt == 3:
                #     break
            else:
                new_lines.append(line)
                bug_lines.append(line)
    with open(new_file_path, 'w') as f:
        f.write(''.join(new_lines))
    with open(bug_file_path, 'w') as f:
        f.write(''.join(bug_lines))
    return new_file_path


def filter_special_cases(file_path):
    '''
    Filter out the special cases of a VCF file. The special cases are not saved in a new file.
    A special case is a variant that has a reference or alternative nucleotide with length greater than 1.
    
    :param file_path:
    :return:
    '''

    new_file_path = 'without_specials_cases_' + file_path # Define the new file path
    filtered_out_file_path = 'special_cases_' + file_path # Define the file path to store the special cases
    new_filtered_lines = list() # List to store the lines filtered out
    new_lines = list() # List to store the new lines

    with open(file_path, 'r') as f: # Open the file
        for line in f.readlines(): # Iterate over the lines of the file

            if not line.startswith('#'): # Check if the line is not a header line nor a comment line (starts with '#')
                fields = line.split('\t') # Split the line by '\t'
                nuc_ref = fields[3] # Get the reference nucleotide
                nuc_alt = fields[4] # Get the alternative nucleotide
                if len(nuc_ref) == 1 and len(nuc_alt) == 1: # Check if the reference and alternative nucleotides have length 1
                    new_lines.append(line) # Append the line to the new lines
                else: # If the reference or alternative nucleotides have length greater than 1
                    new_filtered_lines.append(line) # Append the line to the filtered out lines

            elif line.startswith('#CHROM'): # Check if the line is the header of the file (starts with '#CHROM')
                fields = line[:-1].split('\t') # Remove the '\n' character and split the line by '\t'
                sample_names = fields[-4:] # Get the sample names from the header
                fields.extend([f'%{x}' for x in sample_names]) # Add the new fields to the header
                new_line = '\t'.join(fields) + '\n' # Join the fields and add a '\n' character
                new_lines.append(new_line) # Append the new line to the new lines

            else: # If the line is a comment line
                new_lines.append(line) # Append the line to the new lines


    with open(new_file_path, 'w') as f: # Write the new lines in the new file
        f.write(''.join(new_lines)) # Join the new lines and write in the new file
    with open(filtered_out_file_path, 'w') as f: # Write the filtered out lines in the filtered out file
        f.write(''.join(new_filtered_lines)) # Join the filtered out lines and write in the filtered out file

    return new_file_path # Return the path to the new file


def is_valid(fields, eval_type=0, threshold=None, avg=None, std=None):
    infos = fields[-4:]
    nucs = [x.split(':')[0] for x in infos]
    print(nucs)

    if nucs[0] == nucs[1]:
        print('DESCARTAR: parental_sup_counts ({}) e parental_inf ({}) iguais.'.format(nucs[0], nucs[1]))
        return False

    # ref_nuc = fields[3]
    alelo_sup = fields[3] if nucs[0] == '0' else fields[4]
    print('alelo_sup', nucs[0], alelo_sup)
    nucs = ['A', 'C', 'G', 'T']

    # Get counts of each nucleotide from the parental_sup field of the info
    parental_sup_counts = [int(x) for x in fields[-4].split(':')[4].split(',')] # Get the counts of each nucleotide
    parental_sup_total_reads = float(sum(parental_sup_counts)) # Get the total number of reads

    # Get counts of each nucleotide from the parental_inf field of the info
    parental_alt_counts = [int(x) for x in fields[-3].split(':')[4].split(',')] # Get the counts of each nucleotide
    parental_alt_total_reads = float(sum(parental_alt_counts)) # Get the total number of reads

    # Get counts of each nucleotide from the pool_sup field of the info
    pool_sup_counts = [int(x) for x in fields[-2].split(':')[4].split(',')] # Get the counts of each nucleotide
    pool_sup_total_reads = float(sum(pool_sup_counts)) # Get the total number of reads

    # Get counts of each nucleotide from the pool_rnd field of the info
    pool_rnd_counts = [int(x) for x in fields[-1].split(':')[4].split(',')] # Get the counts of each nucleotide
    pool_rnf_total_reads = float(sum(pool_rnd_counts)) # Get the total number of reads

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

    i_ref = nucs.index(alelo_sup)

    # Necessita haver pelo menos uma read semelhante ao parental superior
    if eval_type == 0:
        if pool_sup_counts[nucs.index(alelo_sup)] > 0:
            print('INCLUIR: existem {} reads semelhantes ao parental_sup_counts ({}).'.format(pool_sup_counts[i_ref], alelo_sup),
                  pool_sup_counts)
        else:
            print('DESCARTAR: existem {} reads semelhantes ao parental_sup_counts ({}).'.format(pool_sup_counts[i_ref], alelo_sup),
                  pool_sup_counts)
            return False

    # Percentual de reads iguais ao parental superior deve ser igual ou maior que o limite
    elif eval_type == 1:

        parental_sup_counts = [int(x) for x in fields[-4].split(':')[4].split(',')]
        parental_sup_total_reads = float(sum(pool_sup_counts))
        # percent_205 = float(parental_sup_counts[i_ref]) / parental_sup_total_reads

        pool_sup_total_reads = float(sum(pool_sup_counts))
        pool_sup_percent = float(pool_sup_counts[i_ref]) / pool_sup_total_reads
        if pool_sup_percent >= threshold:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e superior.'.format(
                pool_sup_percent, alelo_sup), pool_sup_counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e inferior.'.format(
                pool_sup_percent, alelo_sup), pool_sup_counts)
            return False

    # Percentual de reads iguais ao parental superior deve ser o maior
    elif eval_type == 2:
        percents = [float(x) / pool_sup_total_reads for x in pool_sup_counts]
        greather = max(percents)
        percent = percents[i_ref]
        if percent == greather:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. E o maior.'.format(percent,
                                                                                                               alelo_sup),
                  pool_sup_counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. Nao e o maior.'.format(
                percent, alelo_sup), pool_sup_counts)
            return False

    # Percentual de reads iguais ao parental superior deve ser um dos maiores
    elif eval_type == 3:
        percents = [float(x) / pool_sup_total_reads for x in pool_sup_counts]
        greather = max(percents)
        percent = percents[i_ref]
        diff = abs(percent - greather)
        if diff <= threshold and diff > 0.:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. Esta no limite.'.format(
                percent, alelo_sup), pool_sup_counts)
        else:
            print(
                'DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {}. Esta fora do limite.'.format(
                    percent, alelo_sup), pool_sup_counts)
            return False

    elif eval_type == 4:

        parental_sup_counts = [int(x) for x in fields[-4].split(':')[4].split(',')]
        parental_sup_total_reads = float(sum(pool_sup_counts))

        pool_sup_total_reads = float(sum(pool_sup_counts))

        pool_rnd_counts = [int(x) for x in fields[-1].split(':')[4].split(',')]
        pool_rnf_total_reads = float(sum(pool_rnd_counts))

        if parental_sup_total_reads <= 0 or pool_sup_total_reads <= 0 or pool_rnf_total_reads <= 0:
            print('Insuficient number of reads.\n\t205: {} | 207: {} | 208: {}'.format(parental_sup_total_reads, pool_sup_total_reads,
                                                                                       pool_rnf_total_reads))
            return False

        # percent_205 = float(parental_sup_counts[i_ref]) / parental_sup_total_reads
        pool_sup_percent = float(pool_sup_counts[i_ref]) / pool_sup_total_reads
        percent_208 = float(pool_rnd_counts[i_ref]) / pool_rnf_total_reads

        condition_207 = pool_sup_percent >= threshold
        condition_208 = percent_208 >= (avg - std) and percent_208 <= (avg + std)

        if condition_207 and condition_208:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e superior.'.format(
                pool_sup_percent, alelo_sup), pool_sup_counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e inferior.'.format(
                pool_sup_percent, alelo_sup), pool_sup_counts)
            return False

    return True


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


def filter_by_compare(file_path, eval_type=0, threshold=None, avg=None, std=None):
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


# file_path = 'variantsfile.vcf'
file_path = 'variantsfile_Annotated.vcf'

new_file_path = correct_variants(file_path)
new_file_path = filter_special_cases(new_file_path)

# Somente uma ocorrencio do parental_sup
filter_by_compare(new_file_path, eval_type=0)

# Ocorrencia do pool_sup maior que limiar
threshold = .75
filter_by_compare(new_file_path, eval_type=1, threshold=threshold)

# Ocorrencia do parental_sup deve ser a maior
filter_by_compare(new_file_path, eval_type=2)

# Ocorrencia do parental_sup deve ser uma das maiores
threshold = .1
filter_by_compare(new_file_path, eval_type=3, threshold=threshold)

# Ocorrencia do pool_sup maior que limiar e do pool_randomico entre a os desvios
threshold = .75
filter_by_compare(new_file_path, eval_type=4, threshold=threshold, avg=.5, std=.02)


# # Define Main fucntion of the script
# if __name__ == '__main__':
#     # Set and create the outputs folder
#     # outputs_folder = create_outputs_folder()
#
#     # file_path = 'variantsfile.vcf'
#     file_path = 'variantsfile_Annotated.vcf' # Set the path to the VCF file
#
#     new_file_path = correct_variants(file_path)
#     new_file_path = filter_special_cases(new_file_path)