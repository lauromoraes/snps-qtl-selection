import os
import sys
from typing import List

import yaml
import logging

class OutputManager:
    @staticmethod
    def create_outputs_folder(outputs_folder_name: str = 'vcf_filtered') -> str:
        outputs_folder_path = os.path.join(os.getcwd(), outputs_folder_name)
        os.system(f'mkdir -p {outputs_folder_path}')
        return outputs_folder_path


class VCFProcessor:
    def __init__(self, config_file_path: str, logger: logging.Logger) -> None:
        self.logger = logger
        self.config_file_path = config_file_path
        self.config_obj = self.read_config(config_file_path)
        
        self.bug_lines = list()
        self.multiple_alleles_lines = list()
        self.not_single_genotype_lines = list()
        self.equal_genotype_lines = list()
        self.low_reads_lines = list()
        self.insufficient_diff_lines = list()

        

    def read_config(self, config_file_path: str) -> dict:
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

    def select_most_freq_genotype(self, sample_info: str, nuc_ref: str, nuc_alt: str,  verbose: bool = False) -> str:
        fields = sample_info.split(':')
        valids = ['0', '1', '.']
        genotypes = fields[0].split('/')
        flag = True

        if len(genotypes) > 1:
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
                    # self.logger.info(f'Variant with multiple alleles: {fields}')
                    self.multiple_alleles_lines.append(row)
                else:
                    self.not_single_genotype_lines.append(row)
                continue

            parental_inf_info = fields[self.parental_inf_idx].split(':')
            parental_sup_info = fields[self.parental_sup_idx].split(':')
            pool_sup_info = fields[self.pool_sup_idx].split(':')
            pool_rnd_info = fields[self.pool_rnd_idx].split(':')


            # Verify if parental_sup genotype and parental_inf genotype are equal
            genotype_sup = parental_sup_info[0]
            genotype_inf = parental_inf_info[0]

            if genotype_sup == genotype_inf:
                self.equal_genotype_lines.append(row)
                continue


            # Verify if allele counts are less than 3
            samples_fields = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]

            for sample in samples_fields:
                counts = self.get_counts(sample)
                if sum(counts) <= 3:
                    self.low_reads_lines.append(row)
                    continue

            # Verify if the difference in the number of reads is less than 1
            parental_sup_counts = self.get_counts(parental_sup_info)
            parental_inf_counts = self.get_counts(parental_inf_info)
            parental_sup_ordered = sorted(parental_sup_counts, reverse=True)
            parental_inf_ordered = sorted(parental_inf_counts, reverse=True)

            if abs(parental_sup_ordered[0] - parental_inf_ordered[0]) <= 1:
                self.insufficient_diff_lines.append(row)
                continue


            if flag:
                new_line = '\t'.join([str(x) for x in fields])
                selected_rows.append(new_line)
            else:
                self.bug_lines.append(row)

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
        parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info = sample_info
        genotype_sup = [x.split(':')[0] for x in parental_sup_info]
        genotype_inf = [x.split(':')[0] for x in parental_inf_info]
        if genotype_sup == genotype_inf:
            print(f'DISCARD:  parental_sup_counts ({genotype_sup}) and parental_inf ({genotype_inf}) equal.')
            return False

        for sample in sample_info:
            counts = self.get_counts(sample)
            if sum(counts) <= 3:
                print('DISCARD: insuficient number of reads.')
                return False

        parental_sup_counts = self.get_counts(parental_sup_info)
        parental_inf_counts = self.get_counts(parental_inf_info)
        parental_sup_ordered = sorted(parental_sup_counts, reverse=True)
        parental_inf_ordered = sorted(parental_inf_counts, reverse=True)

        if abs(parental_sup_ordered[0] - parental_inf_ordered[0]) <= 1:
            print(f'DISCARD: difference in the number of reads insufficient. parental_sup_counts({parental_sup_counts}) and parental_inf_counts({parental_inf_counts})')
            return False

        return True

    def write_preprocessed_files(self, verbose: bool = False) -> None:
        file_types = ['corrected', 'bug', 'not_single_genotype', 'multiple_alleles', 'equal_genotype', 'low_reads', 'insufficient_diff']
        lines_list = [self.variant_lines, self.bug_lines, self.not_single_genotype_lines, self.multiple_alleles_lines, self.equal_genotype_lines, self.low_reads_lines, self.insufficient_diff_lines]
        for i, lines in enumerate(lines_list):
            basename = os.path.basename(self.vcf_file_path)
            new_file_path = os.path.join(self.output_folder_path, f'{file_types[i]}_{basename}')
            with open(new_file_path, 'w') as f:
                f.write(''.join(self.metadata_lines))
                f.write(''.join(lines))
            if verbose:
                self.logger.info(f'{file_types[i].capitalize()} VCF file written successfully: {new_file_path}')

    def filter(self, verbose: bool = False) -> None:
        filters = self.config_obj['filters']
        for f in filters:
            if f['name'] == 'least_one':
                self.filter_least_one(f['threshold'], verbose=verbose)
            elif f['name'] == 'sup_limit':
                self.filter_sup_limit(f['threshold'], verbose=verbose)
            elif f['name'] == 'greater':
                self.filter_greater(verbose=verbose)
            elif f['name'] == 'relaxed_greater':
                self.filter_relaxed_greater(f['threshold'], verbose=verbose)
            elif f['name'] == 'sup_limit_avg_std':
                self.filter_sup_limit_avg_std(f['threshold'], f['avg'], f['std'], verbose=verbose)
            else:
                self.logger.error(f'Filter {f["name"]} not found.')


    # Method to set the samples indexes in the VCF file using parentals and pools
    def set_samples_idx(self):
        samples_cols = [-4, -3, -2, -1]
        self.parental_sup_idx 


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

    def get_info_fields(self, variant_line: str) -> tuple[str, str, list[list[str]]]:
        fields = variant_line.split('\t')
        allele_ref = fields[3]
        allele_alt = fields[4]
        parental_sup_info = [x for x in fields[self.parental_sup_column].split(':')]
        parental_inf_info = [x for x in fields[self.parental_inf_column].split(':')]
        pool_sup_info = [x for x in fields[self.pool_sup_column].split(':')]
        pool_rnd_info = [x for x in fields[self.pool_rnd_column].split(':')]
        sample_info = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]
        return allele_ref, allele_alt, sample_info

    def get_major_nuc(self, info: str):
        nucs = ['A', 'C', 'G', 'T']
        counts = [int(x) for x in info.split(',')]
        i = counts.index(max(counts))
        return nucs[i], counts

    # def correct_genotype(self, nuc_ref: str, nuc_alt: str, sample_info: str):
    #     fields = sample_info.split(':')
    #     valids = ['0', '1', '.']
    #     ganotypes = fields[0].split('/')
    #     flag = True
    #
    #     if len(ganotypes) > 1:
    #         for genotype in ganotypes:
    #             if genotype not in valids:
    #                 print('Not valid {}'.format(genotype))
    #                 flag = False
    #                 break
    #
    #         if flag:
    #             nuc, counts = self.get_major_nuc(fields[4])
    #             if nuc == nuc_ref:
    #                 fields[0] = '0'
    #             elif nuc == nuc_alt:
    #                 fields[0] = '1'
    #             else:
    #                 print('BUG', fields[0], nuc, nuc_ref, nuc_alt, counts)
    #                 return False
    #
    #     corrected_samples_info = ':'.join(fields)
    #     return corrected_samples_info
    #
    # def correct_variants(self):
    #     new_file_path = 'ployd_corrected_' + self.file_path
    #     bug_file_path = 'bug_' + self.file_path
    #     new_lines = list()
    #     bug_lines = list()
    #     with open(self.file_path, 'r') as f:
    #         for idx, line in enumerate(f.readlines()):
    #             if not line.startswith('#'):
    #                 flag = True
    #                 fields = line.split('\t')
    #                 nuc_ref = fields[3]
    #                 nuc_alt = fields[4]
    #                 if len(nuc_ref) > 1 or len(nuc_alt) > 1:
    #                     new_lines.append(line)
    #                     continue
    #                 for i in range(-4, -2):
    #                     corrected = self.correct_genotype(nuc_ref, nuc_alt, fields[i])
    #                     if corrected:
    #                         fields[i] = corrected
    #                     else:
    #                         flag = False
    #                 if flag:
    #                     new_line = '\t'.join([str(x) for x in fields]) + '\n'
    #                     new_lines.append(new_line)
    #                 else:
    #                     new_lines.append(line)
    #                     bug_lines.append(line)
    #             else:
    #                 new_lines.append(line)
    #                 bug_lines.append(line)
    #     with open(new_file_path, 'w') as f:
    #         f.write(''.join(new_lines))
    #     with open(bug_file_path, 'w') as f:
    #         f.write(''.join(bug_lines))
    #     return new_file_path
    #
    # def filter_special_cases(self):
    #     new_file_path = 'without_specials_cases_' + self.file_path
    #     filtered_out_file_path = 'special_cases_' + self.file_path
    #     new_filtered_lines = list()
    #     new_lines = list()
    #
    #     with open(self.file_path, 'r') as f:
    #         for line in f.readlines():
    #             line = line.strip()
    #             if len(line) == 0:
    #                 continue
    #             if not line.startswith('#'):
    #                 fields = line.split('\t')
    #                 try:
    #                     nuc_ref = fields[3]
    #                 except:
    #                     print(fields)
    #                 nuc_alt = fields[4]
    #                 if len(nuc_ref) == 1 and len(nuc_alt) == 1:
    #                     new_lines.append(line)
    #                 else:
    #                     new_filtered_lines.append(line)
    #             elif line.startswith('#CHROM'):
    #                 fields = line[:-1].split('\t')
    #                 sample_names = fields[-4:]
    #                 fields.extend([f'%{x}' for x in sample_names])
    #                 new_line = '\t'.join(fields) + '\n'
    #                 new_lines.append(new_line)
    #             else:
    #                 new_lines.append(line)
    #
    #     with open(new_file_path, 'w') as f:
    #         f.write('\n'.join(new_lines))
    #     with open(filtered_out_file_path, 'w') as f:
    #         f.write('\n'.join(new_filtered_lines))
    #
    #     return new_file_path

    def get_lines(self):
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
        counts = [int(x) for x in sample_info[4].split(',')]
        return counts


    def least_one(self, pool_sup_info, allele_nuc):
        allele_nuc_idx = self.get_nuc_index(allele_nuc)
        pool_sup_counts = self.get_counts(pool_sup_info)
        if not pool_sup_counts[allele_nuc_idx] > 0:
            return True
        return False

    def percent_threshold(self, parental_sup_info, pool_sup_info, allele_nuc_idx, threshold):
        parental_sup_counts = self.get_counts(parental_sup_info)
        allele_count_on_parental_sup = parental_sup_counts[allele_nuc_idx]
        pool_sup_counts = self.get_counts(pool_sup_info)
        pool_sup_percent = float(pool_sup_counts[allele_nuc_idx]) / float(allele_count_on_parental_sup)

        if pool_sup_percent < threshold:
            print(f'DISCARD: the percentage of reads similar of the allele {allele_nuc_idx} in the pool_sup ({pool_sup_percent}) is less than {threshold}.')
            return True

        pool_sup_total_reads = float(sum(pool_sup_counts))
        pool_sup_percent = float(pool_sup_counts[i_ref]) / pool_sup_total_reads
        if pool_sup_percent >= threshold:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e superior.'.format(pool_sup_percent, alelo_sup), pool_sup_counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantes ao parental_sup_counts ({}) e de {} que e inferior.'.format(pool_sup_percent, alelo_sup), pool_sup_counts)
            return False

    def filter_cases(self, eval_type: int, ref: str, alt: str, sample_info: list[list[str]], **kwargs) -> bool:
        parental_sup_genotype = sample_info[0][0]
        allele_nuc = self.get_allele(parental_sup_genotype, ref, alt)

        if eval_type == 0:
            pool_sup_info = sample_info[2]
            print(allele_nuc, pool_sup_info)
            return self.least_one(pool_sup_info, allele_nuc)
        elif eval_type == 1:
            parental_sup_info = sample_info[0]
            pool_sup_info = sample_info[2]
            return self.percent_threshold(parental_sup_info, pool_sup_info, allele_nuc, kwargs['threshold'])
        else:
            return False

    def set_output_file_name(self, eval_type: int, **kwargs) -> str:
        new_file_path = None
        if eval_type == 0:
            new_file_path = 'filtered_' + 'least_one_' + self.file_path
        elif eval_type == 1:
            new_file_path = 'filtered_' + 'sup_limit_{}_'.format((str(kwargs['threshold']))) + self.file_path
        elif eval_type == 2:
            new_file_path = 'filtered_' + 'greater_' + self.file_path
        elif eval_type == 3:
            new_file_path = 'filtered_' + 'relaxed_greater_{}_'.format((str(kwargs['threshold']))) + self.file_path
        elif eval_type == 4:
            new_file_path = 'filtered_' + 'sup_limit_{}_avg_{}_std_{}_'.format((str(kwargs['threshold'])), str(kwargs['avg']), str(kwargs['std'])) + self.file_path
        return new_file_path

    def get_allele(self, genotype: str, ref: str, alt: str) -> str:
        alelo_sup = ref if genotype == '0' else alt
        return alelo_sup

    def get_nuc_index(self, nuc: str) -> int:
        return ['A', 'C', 'G', 'T'].index(nuc)

    def filter_by_compare(self, eval_type=0, **kwargs):
        new_file_path = self.set_output_file_name(eval_type, **kwargs)
        if not new_file_path:
            return

        metadata_lines, variant_lines = self.get_lines()
        output_lines = metadata_lines

        for line in variant_lines:
            line = line.strip()
            if len(line) == 0:
                continue

            allele_ref, allele_alt, sample_info = self.get_info_fields(line)
            if len(allele_ref) == 1 and len(allele_alt) == 1 and self.valid_line(sample_info):
                filtered_out = self.filter_cases(eval_type=eval_type, ref=allele_ref, alt=allele_alt, sample_info=sample_info, **kwargs)
                if not filtered_out:
                    output_lines.append(line)

        with open(new_file_path, 'w') as f:
            f.write(''.join(output_lines))
