# -*- coding: utf-8 -*-
# !/usr/bin/env python

import os
import sys
import logging
import yaml

from VCFProcessor import *


def print_welcome():
    print('VCF Processor v1.0')
    print('Author: L. Â. G. de Moraes')
    print('Date: 2024-10-01')
    print('Description: This program reads a YAML configuration file with the parameters to filter a VCF file.')
    print('Usage: python main.py <config_file>')
    print('Example: python main.py config.yaml')
    print('-' * 30, end='')


def setup_logger():
    # Create global logger object and set level to INFO to show only INFO messages and above in the console output
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    # Create console handler and set level to INFO
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # Add formatter to console_handler
    console_handler.setFormatter(formatter)
    # Add console_handler to logger
    logger.addHandler(console_handler)


def setup_arguments():
    logging.info('Setting up arguments...')
    logging.info('Number of arguments: %d', len(sys.argv))
    logging.info('Arguments: %s', sys.argv)

    # Check the number of arguments
    if len(sys.argv) != 2:
        logging.error('Invalid number of arguments.')
        sys.exit(1)

    # Get the YAML file
    yaml_file = sys.argv[1]

    # Check if the file exists and is a YAML file
    if not os.path.isfile(yaml_file):
        logging.error('Configuration file does not exist.', yaml_file)
        sys.exit(1)

    # Check if the file is a YAML file
    if not yaml_file.endswith('.yaml') and not yaml_file.endswith('.yml'):
        logging.error('Configuration file is not a YAML file.', yaml_file)
        sys.exit(1)

    return yaml_file


# Function to get the configuration file parameters
def get_config_parameters(config_file):
    logging.info('Getting configuration parameters...')
    # Read the YAML file
    with open(config_file, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
            logging.info('Configuration file loaded successfully.')
            logging.info('Configuration: %s', config)
        except yaml.YAMLError as exc:
            logging.error('Error loading configuration file.', exc)
            sys.exit(1)

    return config


if __name__ == '__main__':
    # Setup logger
    setup_logger()

    # log the start of the program
    logging.info('Program started.')

    # Check the number of arguments and get the VCF file and filter type
    config_file = setup_arguments()

    # lOG the arguments
    logging.info('Configuration file: %s', config_file)

    # Get the configuration parameters
    logging.info('Getting configuration parameters...')
    config = get_config_parameters(config_file)

    # Print welcome message
    print_welcome()

    logging.info('Creating VCFProcessor object...')
    # Create VCFProcessor object
    vcf_processor = VCFProcessor(config_file_path=config_file, logger=logging)
    vcf_processor.load_vcf_file(verbose=True)
    vcf_processor.get_header_info(verbose=True)
    vcf_processor.correct_genotypes(verbose=True)
    vcf_processor.write_preprocessed_files(verbose=True)

    # vcf_processor.filter(verbose=True)

    # Log the end of the program
    logging.info('Program finished successfully.')

    # Exit the program
    sys.exit(0)
