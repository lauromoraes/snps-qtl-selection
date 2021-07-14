def get_major_nuc(info):
    nucs = ['A', 'C', 'G', 'T']
    counts = [int(x) for x in info.split(',')]
    i = counts.index(max(counts))
    return nucs[i], counts


def correct_ployd(nuc_ref, nuc_alt, info):
    fields = info.split(':')
    # print(fields)
    valids = ['0', '1', '.']
    ployds = fields[0].split('/')
    flag = True
    if len(ployds)>1:
        for ployd in ployds:
            if ployd not in valids:
                print('Not valid {}'.format(ployd))
                flag = False
                break
        if flag:
            nuc, counts = get_major_nuc(fields[4])
            # if fields[0] == './.':
            #     print('CATCH', nuc, nuc_ref, nuc_alt)
            if nuc == nuc_ref:
                fields[0] = '0'
            elif nuc == nuc_alt:
                fields[0] = '1'
            else:
                print('BUG', fields[0], nuc, nuc_ref, nuc_alt, counts)
                return False
        corrected_info = ':'.join(fields)
        # print('corrected_info', corrected_info)
    return corrected_info


def correct_variants(file_path):
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
                if len(nuc_ref)>1 or len(nuc_alt)>1:
                    new_lines.append(line)
                    continue
                for i in range(-4,-2):
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
    new_file_path = 'without_specials_' + file_path
    new_lines = list()
    with open(file_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                fields = line.split('\t')
                nuc_ref = fields[3]
                nuc_alt = fields[4]
                if len(nuc_ref)==1 and len(nuc_alt)==1:
                    new_lines.append(line)
            elif line.startswith('#CHROM'):
                fields = line[:-1].split('\t')
                fields.extend([ '%_{}_'.format(x) for x in (205,206,207,208) ])
                print(fields)
                new_line = '\t'.join(fields)+'\n'
                new_lines.append(new_line)
                print(new_line)
            else:
                new_lines.append(line)
    with open(new_file_path, 'w') as f:
        f.write(''.join(new_lines))
    return new_file_path


def is_valid(fields, eval_type=0, threshold=None, avg=None, std=None):
    infos = fields[-4:]
    nucs = [ x.split(':')[0] for x in infos ]
    print(nucs)

    if nucs[0] == nucs[1]:
        print('DESCARTAR: parental_sup ({}) e parental_inf ({}) iguais.'.format(nucs[0], nucs[1]))
        return False

    # ref_nuc = fields[3]
    alelo_sup = fields[3] if nucs[0]=='0' else fields[4]
    print('alelo_sup', nucs[0], alelo_sup)
    nucs = ['A', 'C', 'G', 'T']
    counts = [int(x) for x in fields[-2].split(':')[4].split(',')]
    
    counts_205 = [int(x) for x in fields[-4].split(':')[4].split(',')]
    num_reads_205 = float(sum(counts_205))
    counts_206 = [int(x) for x in fields[-3].split(':')[4].split(',')]
    num_reads_206 = float(sum(counts_206))
    counts_207 = [int(x) for x in fields[-2].split(':')[4].split(',')]
    num_reads_207 = float(sum(counts_207))
    counts_208 = [int(x) for x in fields[-1].split(':')[4].split(',')]
    num_reads_208 = float(sum(counts_208))
    

    num_reads = float(sum(counts))
    if num_reads_205 <= 3. or num_reads_206 <= 3. or num_reads_207 <= 3. or num_reads_208 <= 3.:
        print('DESCARTAR: numero de reads insuficientes.', counts)
        return False

    def has_diff(counts):
    	aux = sorted(counts, reverse=True)
    	return True if abs(aux[0] - aux[1]) <= 1 else False

    if has_diff(counts_205) or has_diff(counts_206):
    	print('DESCARTAR: diferenca no numero de reads insuficientes.', counts_205, counts_206)
    	return False

    i_ref = nucs.index(alelo_sup)

    # Necessita haver pelo menos uma read semelhante ao parental superior
    if eval_type == 0:
        if counts[nucs.index(alelo_sup)] > 0:
            print('INCLUIR: existem {} reads semelhantes ao parental_sup ({}).'.format(counts[i_ref], alelo_sup), counts)
        else:
            print('DESCARTAR: existem {} reads semelhantes ao parental_sup ({}).'.format(counts[i_ref], alelo_sup), counts)
            return False

    # Percentual de reads iguais ao parental superior deve ser igual ou maior que o limite
    elif eval_type == 1:

        counts_205 = [int(x) for x in fields[-4].split(':')[4].split(',')]
        num_reads_205 = float(sum(counts))
        percent_205 = float(counts_205[i_ref]) / num_reads_205

        num_reads_207 = float(sum(counts))
        percent_207 = float(counts[i_ref]) / num_reads_207
        if percent_205 >= threshold and percent_207 >= threshold:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup ({}) e de {} que e superior.'.format(percent_207, alelo_sup), counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantess ao parental_sup ({}) e de {} que e inferior.'.format(percent_207, alelo_sup), counts)
            return False

    # Percentual de reads iguais ao parental superior deve ser o maior
    elif eval_type == 2:
        percents = [ float(x) / num_reads for x in counts ]
        greather = max(percents)
        percent = percents[i_ref]
        if percent == greather:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup ({}) e de {}. E o maior.'.format(percent, alelo_sup), counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantes ao parental_sup ({}) e de {}. Nao e o maior.'.format(percent, alelo_sup), counts)
            return False

    # Percentual de reads iguais ao parental superior deve ser um dos maiores
    elif eval_type == 3:
        percents = [ float(x) / num_reads for x in counts ]
        greather = max(percents)
        percent = percents[i_ref]
        diff = abs(percent - greather)
        if diff <= threshold and diff > 0.:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup ({}) e de {}. Esta no limite.'.format(percent, alelo_sup), counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantes ao parental_sup ({}) e de {}. Esta fora do limite.'.format(percent, alelo_sup), counts)
            return False

    elif eval_type == 4:

        counts_205 = [int(x) for x in fields[-4].split(':')[4].split(',')]
        num_reads_205 = float(sum(counts))

        num_reads_207 = float(sum(counts))

        counts_208 = [int(x) for x in fields[-1].split(':')[4].split(',')]
        num_reads_208 = float(sum(counts_208))

        if num_reads_205 <= 0 or num_reads_207 <= 0 or num_reads_208 <= 0:
            print('Insuficient number of reads.\n\t205: {} | 207: {} | 208: {}'.format(num_reads_205, num_reads_207, num_reads_208))
            return False

        percent_205 = float(counts_205[i_ref]) / num_reads_205
        percent_207 = float(counts[i_ref]) / num_reads_207
        percent_208 = float(counts_208[i_ref]) / num_reads_208

        condition_207 = percent_205 >= threshold and percent_207 >= threshold
        condition_208 = percent_208 >= (avg-std) and percent_208 <= (avg+std)

        if condition_207 and condition_208:
            print('INCLUIR: o percentual de reads semelhantes ao parental_sup ({}) e de {} que e superior.'.format(percent_207, alelo_sup), counts)
        else:
            print('DESCARTAR: o percentual de reads semelhantes ao parental_sup ({}) e de {} que e inferior.'.format(percent_207, alelo_sup), counts)
            return False

    return True

def calc_percents(field):
    counts = [int(x) for x in field.split(':')[4].split(',')]
    num_reads = float(sum(counts))
    percents = [round(x/num_reads,2) for x in counts] if num_reads > 0 else [0 for _ in counts]
    to_str = ','.join([str(x) for x in percents])
    return to_str


def put_percents(fields):
    new_fields = [calc_percents(x) for x in fields[-4:]]
    fields.extend(new_fields)
    new_line = '\t'.join(fields)+'\n'
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
        new_file_path = 'filtered_' + 'sup_limit_{}_avg_{}_std_{}_'.format((str(threshold)), str(avg), str(std)) + file_path

    new_lines = list()
    with open(file_path, 'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                fields = line.rstrip().split('\t')
                nuc_ref = fields[3]
                nuc_alt = fields[4]
                if len(nuc_ref)==1 and len(nuc_alt)==1:
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