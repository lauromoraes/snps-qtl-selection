import os
import sys
from typing import List

import yaml
import logging


class VariantLine(object):

    def __init__(self, line: str, samples_idx: List[int]) -> None:
        self.line = line.strip()
        if len(self.line) == 0:
            raise ValueError('The line must not be empty.')
        self.samples_idx = samples_idx
        if len(self.samples_idx) != 4:
            raise ValueError('The samples index must have 4 elements.')
        self.extract_info()

    def extract_info(self):
        fields = self.line.split('\t')
        samples_idx = self.samples_idx
        self.allele_ref = fields[3]
        self.allele_alt = fields[4]
        self.parental_sup_info = [x for x in fields[samples_idx[0]].split(':')]
        self.parental_inf_info = [x for x in fields[samples_idx[1]].split(':')]
        self.pool_sup_info = [x for x in fields[samples_idx[2]].split(':')]
        self.pool_rnd_info = [x for x in fields[samples_idx[3]].split(':')]
        self.sample_info = [self.parental_sup_info, self.parental_inf_info, self.pool_sup_info, self.pool_rnd_info]
        return True


class OutputManager:
    @staticmethod
    def create_outputs_folder(outputs_folder_name: str = 'vcf_filtered') -> str:
        outputs_folder_path = os.path.join(os.getcwd(), outputs_folder_name)
        os.system(f'mkdir -p {outputs_folder_path}')
        return outputs_folder_path


class VCFProcessor:
    def __init__(self, config_file_path: str, logger: logging.Logger, verbose: bool = False) -> None:
        '''
        Constructor of the VCFProcessor class
        :param config_file_path: the path of the configuration file
        :param logger: the logger object
        :param verbose:
        '''
        self.logger = logger
        self.config_file_path = config_file_path
        self.verbose = verbose

        self.config_obj = self.read_config(config_file_path)
        
        self.bug_lines = list()
        self.multiple_alleles_lines = list()
        self.not_single_genotype_lines = list()
        self.equal_genotype_lines = list()
        self.low_reads_lines = list()
        self.insufficient_diff_lines = list()

        self.variant_lines = None

    def read_config(self, config_file_path: str) -> dict:
        '''
        Method to read the configuration file and set the attributes of the class
        :param config_file_path:
        :return:
        '''
        with open(config_file_path, 'r') as f:
            try:
                config = yaml.safe_load(f)

                self.vcf_file_path = config['vcf_file_path']
                self.output_folder_path = OutputManager.create_outputs_folder(config['output_folder_path'])

                self.parental_sup_column = config['parental_sup_column']
                self.parental_inf_column = config['parental_inf_column']
                self.pool_sup_column = config['pool_sup_column']
                self.pool_rnd_column = config['pool_rnd_column']

                self.filter_names = [x['name'] for x in config['filters']]

                self.logger.info(f'Config file read successfully: {config_file_path}')

            except yaml.YAMLError as exc:
                print(exc)
                self.logger.error(f'Error reading config file: {config_file_path}')
                sys.exit(1)
        return config

    def load_vcf_file(self, verbose: bool = False) -> None:
        '''
        Method to read the VCF file and separate the metadata lines from the variant lines
        :param verbose:
        :return:
        '''
        with open(self.vcf_file_path, 'r') as f:
            try:
                lines = f.readlines()
                self.metadata_lines = [line for line in lines if line.startswith('#')]
                self.variant_lines = [line for line in lines if not line.startswith('#')]
                if verbose:
                    self.logger.info(f'VCF file read successfully: {self.vcf_file_path}')
                    self.logger.info(f'VCF file metadata lines: {len(self.metadata_lines)}')
                    self.logger.info(f'VCF file variant lines: {len(self.variant_lines)}')
            except Exception as e:
                self.logger.error(f'Error reading VCF file: {self.vcf_file_path}')
                sys.exit(1)

    def select_most_freq_genotype(self, sample_info: str, nuc_ref: str, nuc_alt: str, verbose: bool = False) -> str:
        '''
        Method to select the most frequent genotype in the sample info field and correct it to 0 or 1
        :param sample_info:
        :param nuc_ref:
        :param nuc_alt:
        :param verbose:
        :return:
        '''
        fields = sample_info.split(':')
        valids = ['0', '1', '.']
        genotypes = fields[0].split('/')
        flag = True

        for genotype in genotypes:
            if genotype not in valids:
                if verbose:
                    self.logger.info(f'Invalid genotype: {genotype}')
                flag = False
                break

        if flag:
            nuc, counts = self.get_major_nuc(fields[4])
            if nuc == nuc_ref:
                fields[0] = '0'
            elif nuc == nuc_alt:
                fields[0] = '1'
            else:
                if verbose:
                    self.logger.info(f'BUG: {fields[0]}, {nuc}, {nuc_ref}, {nuc_alt}, {counts}')
                return False

        corrected_samples_info = ':'.join(fields)
        return corrected_samples_info

    def correct_genotypes(self, verbose: bool = False) -> list[str]:
        '''
        Method to correct the genotypes of the variant lines
        :param verbose:
        :return:
        '''
        selected_rows = list()
        for i, row in enumerate(self.variant_lines):
            flag = True
            fields = row.split('\t')
            nuc_ref = fields[self.ref_idx]
            nuc_alt = fields[self.alt_idx]
            if len(nuc_ref) == 1 and len(nuc_alt) == 1:
                for j in range(-4, -2):
                    sample = fields[j]
                    corrected = self.select_most_freq_genotype(sample, nuc_ref, nuc_alt, verbose=verbose)
                    if corrected:
                        fields[j] = corrected
                    else:
                        flag = False
            else:
                if len(nuc_ref.split(',')) > 1 or len(nuc_alt.split(',')) > 1:
                    # self.logger.info(f'Variant with multiple alleles: {nuc_ref} - {nuc_alt}')
                    self.multiple_alleles_lines.append(row)
                else:
                    self.not_single_genotype_lines.append(row)
                continue

            new_line = '\t'.join([str(x) for x in fields])

            parental_inf_info = fields[self.parental_inf_idx].split(':')
            parental_sup_info = fields[self.parental_sup_idx].split(':')
            pool_sup_info = fields[self.pool_sup_idx].split(':')
            pool_rnd_info = fields[self.pool_rnd_idx].split(':')


            # Verify if parental_sup genotype and parental_inf genotype are equal
            genotype_sup = parental_sup_info[0]
            genotype_inf = parental_inf_info[0]

            if genotype_sup == genotype_inf:
                self.equal_genotype_lines.append(new_line)
                continue


            # Verify if allele counts are less than 3
            samples_fields = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]

            min_flag = False
            for sample in samples_fields:
                counts = self.get_counts(sample)
                if sum(counts) <= 3:
                    min_flag = True

            if min_flag:
                self.low_reads_lines.append(new_line)
                continue

            # Verify if the difference in the number of reads is less than 1
            parental_sup_counts = self.get_counts(parental_sup_info)
            parental_inf_counts = self.get_counts(parental_inf_info)
            parental_sup_ordered = sorted(parental_sup_counts, reverse=True)
            parental_inf_ordered = sorted(parental_inf_counts, reverse=True)

            # # TODO: review this case
            # if abs(parental_sup_ordered[0] - parental_inf_ordered[0]) <= 1:
            #     self.insufficient_diff_lines.append(new_line)
            #     continue

            # TODO: review this new case
            min_diff = 0.2
            diff_perc_sup = abs(parental_sup_ordered[0] - parental_sup_ordered[1]) / sum(parental_sup_counts)
            diff_perc_inf = abs(parental_inf_ordered[0] - parental_inf_ordered[1]) / sum(parental_inf_counts)

            if diff_perc_sup < min_diff or diff_perc_inf < min_diff:
                self.insufficient_diff_lines.append(new_line)
                continue

            if flag:
                selected_rows.append(new_line)
            else:
                self.bug_lines.append(new_line)

        self.variant_lines = selected_rows

        if verbose:
            self.logger.info(f'Corrected genotypes successfully')
            self.logger.info(f'Number of corrected lines: {len(self.variant_lines)}')
            self.logger.info(f'Number of not single genotype lines: {len(self.not_single_genotype_lines)}')
            self.logger.info(f'Number of multiple alleles lines: {len(self.multiple_alleles_lines)}')
            self.logger.info(f'Number of equal genotype lines: {len(self.equal_genotype_lines)}')
            self.logger.info(f'Number of low reads lines: {len(self.low_reads_lines)}')
            self.logger.info(f'Number of insufficient diff lines: {len(self.insufficient_diff_lines)}')
            self.logger.info(f'Number of bug lines: {len(self.bug_lines)}')


        return selected_rows

    def valid_line(self, sample_info: list[str]) -> bool:
        '''
        Method to verify if the line of a VCF file is valid
        :param sample_info: info of the samples
        :return: True if the line is valid, False otherwise
        '''
        parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info = sample_info
        genotype_sup = [x.split(':')[0] for x in parental_sup_info]
        genotype_inf = [x.split(':')[0] for x in parental_inf_info]

        # Genotype verification of the parental_sup and parental_inf fields of the info
        if genotype_sup == genotype_inf:
            print(f'DISCARD:  parental_sup_counts ({genotype_sup}) and parental_inf ({genotype_inf}) equal.')
            return False

        # Minimum number of reads verification
        for sample in sample_info:
            counts = self.get_counts(sample)
            if sum(counts) <= 3:
                print('DISCARD: insuficient number of reads.')
                return False

        parental_sup_counts = self.get_counts(parental_sup_info)
        parental_inf_counts = self.get_counts(parental_inf_info)
        parental_sup_ordered = sorted(parental_sup_counts, reverse=True)
        parental_inf_ordered = sorted(parental_inf_counts, reverse=True)

        # Check if the counts of the parental_sup and parental_inf fields of the info have a difference of at most 1
        if abs(parental_sup_ordered[0] - parental_inf_ordered[0]) <= 1:
            print(f'DISCARD: difference in the number of reads insufficient. parental_sup_counts({parental_sup_counts}) and parental_inf_counts({parental_inf_counts})')
            return False

        return True

    def write_preprocessed_files(self, verbose: bool = False) -> None:
        '''
        Method to write the preprocessed files in the output folder path
        :param verbose: if True, print the log messages
        :return: None
        '''
        file_types = ['corrected', 'bug', 'not_single_genotype', 'multiple_alleles', 'equal_genotype', 'low_reads', 'insufficient_diff']
        lines_list = [self.variant_lines, self.bug_lines, self.not_single_genotype_lines, self.multiple_alleles_lines, self.equal_genotype_lines, self.low_reads_lines, self.insufficient_diff_lines]
        for i, lines in enumerate(lines_list):
            basename = os.path.basename(self.vcf_file_path)
            preprocessed_folder = os.path.join(self.output_folder_path, 'preprocessed')
            os.system(f'mkdir -p {preprocessed_folder}')
            new_file_path = os.path.join(preprocessed_folder, f'{file_types[i]}_{basename}')
            with open(new_file_path, 'w') as f:
                f.write(''.join(self.metadata_lines))
                f.write(''.join(lines))
            if verbose:
                self.logger.info(f'{file_types[i].capitalize()} VCF file written successfully: {new_file_path}')


    def get_header_info(self, verbose: bool = False) -> str:
        header_line = self.metadata_lines[-1].strip()
        self.header_fields = header_line.split('\t')
        self.header_idx = {field: idx for idx, field in enumerate(self.header_fields)}

        self.parental_inf_idx = self.header_idx[self.parental_inf_column]
        self.parental_sup_idx = self.header_idx[self.parental_sup_column]
        self.pool_sup_idx = self.header_idx[self.pool_sup_column]
        self.pool_rnd_idx = self.header_idx[self.pool_rnd_column]

        self.alt_idx = self.header_idx['ALT']
        self.ref_idx = self.header_idx['REF']
        self.format = self.header_idx['FORMAT']

        if verbose:
            self.logger.info(f'VCFProcessor.get_header_info - Header fields: {self.header_fields}')
            self.logger.info(f'VCFProcessor.get_header_info - Parental sup column {self.parental_sup_column} and index: {self.parental_sup_idx}')
            self.logger.info(f'VCFProcessor.get_header_info - Parental inf column {self.parental_inf_column} and index: {self.parental_inf_idx}')
            self.logger.info(f'VCFProcessor.get_header_info - Pool sup column {self.pool_sup_column} and index: {self.pool_sup_idx}')
            self.logger.info(f'VCFProcessor.get_header_info - Pool rnd column {self.pool_rnd_column} and index: {self.pool_rnd_idx}')
            self.logger.info(f'VCFProcessor.get_header_info - Alt column and index: {self.alt_idx}')
            self.logger.info(f'VCFProcessor.get_header_info - Ref column and index: {self.ref_idx}')
            self.logger.info(f'VCFProcessor.get_header_info - Format column and index: {self.format}')
        
        return header_line

    # def get_info_fields(self, variant_line: str) -> tuple[str, str, list[list[str]]]:
    #     fields = variant_line.split('\t')
    #     allele_ref = fields[3]
    #     allele_alt = fields[4]
    #     parental_sup_info = [x for x in fields[self.parental_sup_column].split(':')]
    #     parental_inf_info = [x for x in fields[self.parental_inf_column].split(':')]
    #     pool_sup_info = [x for x in fields[self.pool_sup_column].split(':')]
    #     pool_rnd_info = [x for x in fields[self.pool_rnd_column].split(':')]
    #     sample_info = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]
    #     return allele_ref, allele_alt, sample_info

    def get_major_nuc(self, info: str):
        '''
        Method to get the nucleotide with the highest count in the info field
        
        :param info: the info string with the counts of the nucleotides
        :return: the nucleotide with the highest count and the counts of the nucleotides
        '''
        nucs = ['A', 'C', 'G', 'T']
        counts = [int(x) for x in info.split(',')]
        i = counts.index(max(counts))
        return nucs[i], counts

    def get_lines(self):
        '''
        Method to get the metadata and variant lines of a VCF file

        :return: the metadata lines and the variant lines
        '''
        metadata_lines = list()
        variant_lines = None

        with open(self.file_path, 'r') as f:
            variant_lines = [line for line in f.readlines()]
            for i, line in enumerate(variant_lines):
                if line.startswith('#'):
                    metadata_lines.append(line)
                else:
                    break
            variant_lines = variant_lines[i:]

        return metadata_lines, variant_lines

    def get_info_fields(self, variant_line: str, samples_idx: list[int] = [-4, -3, -2, -1]) -> tuple[str, str, list[list[str]]]:
        '''
        Method to extract and get the info fields of a variant line

        :param variant_line: the variant line of a VCF file to extract the info fields
        :param samples_idx: the indexes of the samples types in the variant line
        :return: the reference allele, the alternative allele and the sample info
        '''
        # TODO: create and return a new object from the class VariantLine (to be created)
        fields = variant_line.split('\t')
        allele_ref = fields[3]
        allele_alt = fields[4]
        parental_sup_info = [x for x in fields[samples_idx[0]].split(':')]
        parental_inf_info = [x for x in fields[samples_idx[1]].split(':')]
        pool_sup_info = [x for x in fields[samples_idx[2]].split(':')]
        pool_rnd_info = [x for x in fields[samples_idx[3]].split(':')]
        sample_info = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]
        return allele_ref, allele_alt, sample_info

    def get_counts(self, sample_info: list[str]) -> list[int]:
        '''
        Method to get the counts of a single sample info

        :param sample_info: the sample info to get the counts
        :return: the counts of the sample info
        '''
        counts = [int(x) for x in sample_info[4].split(',')]
        return counts


    def at_least(self, n_reads: int = 1) -> List[str]:
        output_lines = list()
        filtered_lines = list()
        lines = self.variant_lines

        if self.verbose:
            self.logger.info(f'Filtering variants with at least {n_reads} reads.')

        for l in lines:
            # Extract the info fields
            allele_ref, allele_alt, sample_info = self.get_info_fields(l)

            # Parental superior info
            parental_sup_info = sample_info[0]
            parental_sup_allele, parental_sup_counts = self.get_major_nuc(parental_sup_info[4])

            # Pool superior info
            pool_sup_info = sample_info[2]
            pool_sup_allele, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])

            # Verify if the parental_sup allele is the same as the pool_sup allele
            if parental_sup_allele != pool_sup_allele:
                filtered_lines.append(l)
                continue

            # Verify if the allele count on the pool_sup is greater than the minimum number of reads
            if max(pool_sup_counts) >= n_reads:
                output_lines.append(l)
            else:
                filtered_lines.append(l)

        return output_lines, filtered_lines

    # def get_allele(self, sample_info: list[str], allele_ref: str, allele_alt: str) -> str:
    #     '''
    #     Method to get the allele of a sample info based on the reference and alternative alleles
    #     :param sample_info: the sample info to get the allele
    #     :param allele_ref: nucleotide of the reference allele
    #     :param allele_alt: nucleotide of the alternative allele
    #     :return: the allele of the sample info
    #     '''
    #     allele_cod = sample_info[0]
    #     if allele_cod == '0':
    #         return allele_ref
    #     elif allele_cod == '1':
    #         return allele_alt
    #     return None

    def percent_threshold(self, threshold: float = 0.0) -> List[str]:

        if self.verbose:
            self.logger.info(f'Filtering variants with percentage of reads greater than {threshold}.')

        output_lines = list()
        filtered_lines = list()
        lines = self.variant_lines


        for l in lines:
            # Extract the info fields
            allele_ref, allele_alt, sample_info = self.get_info_fields(l)

            # Parental superior info
            parental_sup_info = sample_info[0]

            # Pool superior info
            pool_sup_info = sample_info[2]
            pool_sup_counts = self.get_counts(pool_sup_info)
            pool_sup_total_reads = sum(pool_sup_counts)

            # Allele and counts of the parental sup and pool sup
            parental_sup_allele_nuc, _ = self.get_major_nuc(parental_sup_info[4])
            pool_sup_alelle_nuc, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])

            # Verify if the pool_sup allele is the same as the parental_sup allele
            if parental_sup_allele_nuc != pool_sup_alelle_nuc:
                filtered_lines.append(l)
                continue

            # Calculate the ratio of the pool_sup allele
            pool_sup_allele_ratio = float(max(pool_sup_counts)) / float(pool_sup_total_reads)

            # Verify if the ratio is greater than the threshold
            if pool_sup_allele_ratio >= threshold:
                output_lines.append(l)
            else:
                filtered_lines.append(l)

        return output_lines, filtered_lines

    def parental_sup_greater(self) -> List[str]:

        if self.verbose:
            self.logger.info(f'Filter out variants that do not have the reference allele count greater than other alleles')

        output_lines = list()
        filtered_lines = list()
        lines = self.variant_lines


        for l in lines:
            allele_ref, allele_alt, sample_info = self.get_info_fields(l)

            # Parental superior info
            parental_sup_info = sample_info[0]

            # Pool superior info
            pool_sup_info = sample_info[2]

            # Allele and counts of the parental sup and pool sup
            parental_sup_allele_nuc, _ = self.get_major_nuc(parental_sup_info[4])
            pool_sup_alelle_nuc, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])

            # Verify if the pool_sup allele is the same as the parental_sup allele
            if parental_sup_allele_nuc != pool_sup_alelle_nuc:
                filtered_lines.append(l)
                continue

            if max(pool_sup_counts) == pool_sup_counts[self.get_nuc_index(pool_sup_alelle_nuc)]:
                output_lines.append(l)
            else:
                filtered_lines.append(l)

            # pool_sup_counts = self.get_counts(sample_info[2])
            # pool_sup_total_reads = sum(pool_sup_counts)
            #
            # # TODO: Tem que ser perecentual? Nao pode ser ser a contagem? A principio, nao havera diferenca
            # pool_sup_percents = [float(x) / float(pool_sup_total_reads) for x in pool_sup_counts]
            # if pool_sup_percents[self.get_nuc_index(allele_ref)] == max(pool_sup_percents):
            #     output_lines.append(l)
            # else:
            #     filtered_lines.append(l)

        return output_lines, filtered_lines


    def diff_from_greater(self, diff_max: float = 0.1, diff_min: float = None) -> List[str]:

        if self.verbose:
            self.logger.info(f'Accept variants that have the superior parental allele on superior pool frequency '
                             f'closer to the greatest frequency')

        output_lines = list()
        filtered_lines = list()
        lines = self.variant_lines


        for l in lines:
            allele_ref, allele_alt, sample_info = self.get_info_fields(l)

            # Parental superior info
            parental_sup_info = sample_info[0]

            # Pool superior info
            pool_sup_info = sample_info[2]

            # Allele and counts of the parental sup and pool sup
            parental_sup_allele_nuc, _ = self.get_major_nuc(parental_sup_info[4])
            pool_sup_alelle_nuc, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])
            pool_sup_total_reads = sum(pool_sup_counts)

            # Verify if the pool_sup allele is the same as the parental_sup allele
            if parental_sup_allele_nuc != pool_sup_alelle_nuc:
                filtered_lines.append(l)
                continue

            # Calculate the ratio of the pool_sup allele
            pool_sup_percents = [float(x) / float(pool_sup_total_reads) for x in pool_sup_counts]
            # Find the nucleotide with the greatest percentage
            greater_perc = max(pool_sup_percents)
            # Find the percentage of the parental superior allele
            pool_sup_allele_perc = pool_sup_percents[self.get_nuc_index(pool_sup_alelle_nuc)]
            # Calculate the difference between the greatest percentage and the parental superior allele percentage
            diff = abs(greater_perc - pool_sup_allele_perc)
            
            # TODO: Incluir o zero?

            # Verify if the difference is smaller than the threshold
            if diff_max >= diff > 0.0:
                output_lines.append(l)
            else:
                filtered_lines.append(l)

        return output_lines, filtered_lines

    def rnd_mean(self, threshold: float = 0.0, avg: float = 0.0, std: float = 0.0) -> List[str]:

        if self.verbose:
            self.logger.info(f'Verify if random pool allele equal to superior parental has a mean percentage is between the average and standard deviation')

        output_lines = list()
        filtered_lines = list()
        lines = self.variant_lines

        for l in lines:
            allele_ref, allele_alt, sample_info = self.get_info_fields(l)

            # Parental superior info
            parental_sup_info = sample_info[0]
            parental_sup_allele_nuc, _ = self.get_major_nuc(parental_sup_info[4])

            # TODO: Tem que verificar se allelo parental sup e pool sup sao iguais??

            pool_sup_counts = self.get_counts(sample_info[2])
            pool_rnd_counts = self.get_counts(sample_info[3])

            pool_sup_total_reads = sum(pool_sup_counts)
            pool_rnd_total_reads = sum(pool_rnd_counts)

            if pool_rnd_total_reads == 0 or pool_sup_total_reads == 0:
                self.logger.error(f'Error: pool_rnd_total_reads: {pool_rnd_total_reads} or pool_sup_total_reads: {pool_sup_total_reads} is zero.')
                filtered_lines.append(l)
                continue

            pool_sup_percents = [float(x) / float(pool_sup_total_reads) for x in pool_sup_counts]
            pool_rnd_percents = [float(x) / float(pool_rnd_total_reads) for x in pool_rnd_counts]

            ref_sup_perc = pool_sup_percents[self.get_nuc_index(parental_sup_allele_nuc)]
            ref_rnd_perc = pool_rnd_percents[self.get_nuc_index(parental_sup_allele_nuc)]

            cond1 = ref_sup_perc >= threshold
            cond2 = (avg + std) >= ref_rnd_perc >= (avg - std)

            if cond1 and cond2:
                output_lines.append(l)
            else:
                filtered_lines.append(l)

        return output_lines, filtered_lines


    def get_nuc_index(self, nuc: str) -> int:
        return ['A', 'C', 'G', 'T'].index(nuc)

    def apply_filters(self, verbose: bool = True) -> None:
        filters = self.config_obj['filters']


        for f in filters:
            self.logger.info('=' * 50)
            
            filter_name = f['name']
            output_lines = None

            if verbose:
                self.logger.info('< FILTER >')
                self.logger.info(f'{filter_name}')
                self.logger.info('< PARAMETERS >')
                self.logger.info(f'{f}')

            # Identify the filter to apply based on the name, get the parameters and apply the filter
            if filter_name == 'at_least':
                n_reads = int(f['n_reads'])
                output_lines, filtered_lines = self.at_least(n_reads=n_reads)

            elif filter_name == 'percent_threshold':
                threshold = float(f['threshold']) / 100.
                output_lines, filtered_lines = self.percent_threshold(threshold=threshold)

            elif filter_name == 'ref_greater':
                output_lines, filtered_lines = self.parental_sup_greater()

            elif filter_name == 'diff_from_greater':
                diff_max = float(f['diff_max']) / 100.
                output_lines, filtered_lines = self.diff_from_greater(diff_max=diff_max)

            elif filter_name == 'rnd_mean':
                threshold = float(f['threshold']) / 100.
                avg = float(f['avg_count']) / 100.
                std = float(f['std_dev']) / 100.
                output_lines, filtered_lines = self.rnd_mean(threshold=threshold, avg=avg, std=std)


            else:
                self.logger.error(f'Filter {filter_name} not found.')
                continue

            if verbose:
                self.logger.info('< SUMMARY >')
                self.logger.info(f'Number of lines before filter: {len(self.variant_lines)}')
                self.logger.info(f'Number of lines after filter: {len(output_lines)}')


            # Write the filtered VCF file
            if output_lines and (len(output_lines) > 0):
                output_file_path = os.path.join(self.output_folder_path, f'filtered_{filter_name}_{os.path.basename(self.vcf_file_path)}')
                with open(output_file_path, 'w') as f:
                    f.write(''.join(self.metadata_lines))
                    f.write(''.join(output_lines))

                filtered_file_path = os.path.join(self.output_folder_path, f'discarded_{filter_name}_{os.path.basename(self.vcf_file_path)}')
                with open(filtered_file_path, 'w') as f:
                    f.write(''.join(self.metadata_lines))
                    f.write(''.join(filtered_lines))
                if verbose:
                    self.logger.info('< OUTPUTS >')
                    self.logger.info(f'Filtered variants written successfully to {output_file_path}')
                    self.logger.info(f'Discarded variants written to {filtered_file_path}')
            else:
                self.logger.warning(f'No lines to write in the filtered VCF file: {filter_name}_{os.path.basename(self.vcf_file_path)}')

            self.logger.info('-' * 50)
