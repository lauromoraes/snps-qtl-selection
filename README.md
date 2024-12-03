# SNP Selection Workflow from Sequencing Data

## Overview
This project includes a basic workflow for processing a [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format) file (`.vcf`) as input and generates additional *vcf* files by filtering variants based on various criteria.


## Config files



Example:
```
vcf_file_path: "./variantsfile_Annotated.vcf"
output_folder_path: "./vcf_filtered"

parental_sup_column: "Samplen206_"
parental_inf_column: "Samplen205_"
pool_sup_column: "Samplen207_"
pool_rnd_column: "Samplen208_"

filters:
  - name: "at_least"
    description: "Filter out variants that do not have the reference allele count greater than or equal to the specified number"
    n_reads: "1"
  - name: "percent_threshold"
    description: "Filter out variants that the reference allele count is less than the specified percentage of the total allele count"
    threshold: "75"
  - name: "ref_greater"
    description: "Filter out variants that do not have the reference allele count greater than other alleles"
  - name: "diff_from_greater"
    description: "Filter out variants that do not have the reference allele count closer to the greater allele count"
    diff_max: "10"
  - name: "rnd_mean"
    description: "Filter out variants"
    threshold: "75"
    avg_count : "50"
    std_dev : "2"
```

# Usage

Example of how to run the program using the `config.yaml` configuration file

```
python main.py config.yml
```
