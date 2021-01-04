import pandas as pd
import numpy as np
import gseapy as gp
import mygene
import sys
import time
import seaborn as sns
mg = mygene.MyGeneInfo()


# WORKING WITH _impact.txt FILES (extraction of gene names impacted at
# diferent levels, intersection of gene lists, generation of tables with
# gene names and their associated GO
def read_impact_data(path):
    data = pd.read_csv(
        path,
        sep='\t',
        usecols=[
            'CHROM',
            'POS',
            'REF',
            'ALT',
            'IMPACT',
            'GENE',
            'EFFECT',
            'LOF_PERC',
            'NMD_PERC'],
        names=[
            'CHROM',
            'POS',
            'REF',
            'ALT',
            'IMPACT',
            'GENE',
            'EFFECT',
            'LOF_GENE',
            'LOF_PERC',
            'NMD_GENE',
            'NMD_PERC'],
        header=0)
    return data


def subset_byimpact(data, impactlist):
    """returns a dataframe - a subset of data based on a list
    with snpEff impact values: 'MODIFIER', 'LOW', 'MODERATE', 'HIGH'"""
    reslist = []
    for i in impactlist:
        print(i)
        portion = data[data['IMPACT'] == i]
        reslist.append(portion)
    subset = pd.concat(reslist)
    return subset


def add_genenames(data):
    """adds names of genes based on the 'GENE' column"""
    genelist = []
    for i in data['GENE']:
        genelist.append(i.split('_')[0])
    data['GENE_NAME'] = genelist
    return data


def write_genenames(path, genelist):
    """writes names of genes from a list into a file;
    apply after generate_genelist"""
    with open(path, 'w') as file:
        for i in genelist:
            line = i + "\n"
            file.writelines(line)


def generate_genelist(selection):
    """generates a list of gene names without duplications;
    apply after add_genenames"""
    genelist = []
    for i in selection['GENE_NAME']:
        if i not in genelist:
            genelist.append(i)
    return genelist


def go_get_genenames(path, impactlist, write=True, generate_list=False):
    """generates a list of or writes to a .txt file names of the genes with the specified impact from snpEff annotation,
    'MODIFIER', 'LOW', 'MODERATE', 'HIGH',
    given in a list of the desired impacts"""
    data = read_impact_data(path)
    subset = subset_byimpact(data, impactlist)
    subset_gnames = add_genenames(subset)
    genelist = generate_genelist(subset_gnames)
    if write:
        ssname = ''
        for i in impactlist:
            if i != 'MODIFIER':
                ssname += i[0]
            else:
                ssname += 'B'
        newpath = path[:-4] + '_genenames_' + ssname + '.txt'
        write_genenames(newpath, genelist)
    return genelist


def go_get_genenames_all_impact_levels(path):
    """writes gene names into a .txt file
    with all possible levels of filtration based on impact"""
    list_of_impactlists = [['MODIFIER', 'LOW', 'MODERATE', 'HIGH'], [
        'LOW', 'MODERATE', 'HIGH'], ['MODERATE', 'HIGH'], ['HIGH']]
    for i in list_of_impactlists:
        go_get_genenames(path, i)


def intersect_2_genelists(list1, list2):
    """returns a list of genes present in both gene lists"""
    intersection = []
    for i in list1:
        if i in list2:
            intersection.append(i)
    return intersection


def print_genelist(genelist):
    for i in genelist:
        print(i)


class NoOutputOptionSpecified(Exception):
    pass


def get_go_annotations(genelist):
    global annotations_all
    geneanno = {}
    for i in genelist:
        anno = ''
        if i in annotations_all:
            geneanno[i] = annotations_all[i]
        else:
            while anno == '':
                try:
                    anno = mg.getgene(i, fields='go')
                    break
                except BaseException:
                    print('Connection refused be the server...')
                    print('ZZZZzzzz')
                    time.sleep(5)
                    print('I am awaken again!')
                    continue
            if anno is not None:
                if 'go' in anno:
                    geneanno[i] = anno['go']
                    annotations_all[i] = anno['go']
                else:
                    print('No annotation for the following gene ise found:', i)
                    geneanno[i] = {'NA'}
            else:
                geneanno[i] = {'NA'}
    return geneanno


def generate_annotation_summary(geneanno):
    """generates a summary table from
    the gene ontology annotation"""
    summary_list = []
    for i in geneanno:
        gene = i
        if geneanno[i] != {'NaN'}:
            for j in geneanno[i]:
                gocategory = j
                if isinstance(geneanno[i][j], dict):
                    goid = geneanno[i][j]['id']
                    goterm = geneanno[i][j]['term']
                    goevidence = geneanno[i][j]['evidence']
                    summary_list.append(
                        [gene, gocategory, goid, goterm, goevidence])
                else:
                    for x in range(len(geneanno[i][j])):
                        goid = geneanno[i][j][x]['id']
                        goterm = geneanno[i][j][x]['term']
                        goevidence = geneanno[i][j][x]['evidence']
                        summary_list.append(
                            [gene, gocategory, goid, goterm, goevidence])
    summary_df = pd.DataFrame(
        summary_list,
        columns=[
            'GENE_NAME',
            'GO_CATEGORY',
            'GO_ID',
            'GO_TERM',
            'GO_EVIDENCE'])
    return summary_df


def write_go_annotation_table(path, impactlist, summary_df):
    '''writes the go annotation table to a .csv file'''
    ssname = ''
    for i in impactlist:
        if i != 'MODIFIER':
            ssname += i[0]
        else:
            ssname += 'B'
    newpath = path[:-4] + '_genenames_' + ssname + '_gotable.csv'
    summary_df.to_csv(newpath, index=False)


def get_go_annotations_init():
    instruction = input(
        'input format: filepath impact(-a, -b, -l, -m, -h, -lm etc.) mode(-ggol, -ggot, -wgl, -wgot): ')
    instruction = instruction.split(' ')
    """first a path (a name of the file) must be specified, then - either -a to generate .txt files with
gene names for all possible levels of impact filtration or specify the desired filtration
level explicitly: -b-for modifier impact, -l-for low impact, -m-for moderate impact and -h-for high impact,
-lm, -ml, -lmh etc. options are available too,
this level should be written right after the path, if need to return list of genes or annotation table, use explicit definition;
then use:
-wgl to write the genes to a .txt file AND/OR
-wgot to write the gene names together with their go id and corresponding terms
-ggl to return the list of genes OR
-ggot to return the dataframe with go annotated genes"""
    ggl = False
    ggot = True
    wgl = True
    wgot = True
    global annotations_all
    annotations_all = {}
    if instruction[1] == '-a':
        go_get_genenames_all_impact_levels(instruction[0])
    else:
        impact_levels = []
        for i in instruction[1][1:]:
            if i == 'b':
                impact_levels.append('MODIFIER')
            elif i == 'l':
                impact_levels.append('LOW')
            elif i == 'm':
                impact_levels.append('MODERATE')
            elif i == 'h':
                impact_levels.append('HIGH')
            else:
                raise NoOutputOptionSpecified(
                    'Specify levels of impact: -a-for all possible levels of filtration, -b-for modifier impact, -l-for low impact, -m-for moderate impact and -h-for high impact')
        if '-wgl' not in instruction:
            wgl = False
        elif '-wgot' not in instruction:
            wgot = False
        elif '-ggl' in instruction:
            ggl = True
            ggot = False
        glist = go_get_genenames(
            instruction[0],
            impactlist=impact_levels,
            write=wgl,
            generate_list=ggl)
        if ggot or wgot:
            ann = get_go_annotations(glist)
            summary = generate_annotation_summary(ann)
            if wgot:
                write_go_annotation_table(
                    instruction[0],
                    impactlist=impact_levels,
                    summary_df=summary)
        if ggl:
            return glist
        elif ggot:
            return summary


def get_go_annotations_RAW_init():
    instruction = vars(__builtins__).get('raw_input', input)
    instruction = sys.argv[1:]
    ggl = False
    ggot = True
    wgl = True
    wgot = True
    global annotations_all
    annotations_all = {}
    if instruction[1] == '-a':
        go_get_genenames_all_impact_levels(instruction[0])
    else:
        impact_levels = []
        for i in instruction[1][1:]:
            if i == 'b':
                impact_levels.append('MODIFIER')
            elif i == 'l':
                impact_levels.append('LOW')
            elif i == 'm':
                impact_levels.append('MODERATE')
            elif i == 'h':
                impact_levels.append('HIGH')
            else:
                raise NoOutputOptionSpecified(
                    'Specify levels of impact: -a-for all possible levels of filtration, -b-for modifier impact, -l-for low impact, -m-for moderate impact and -h-for high impact')
        if '-wgl' not in instruction:
            wgl = False
        elif '-wgot' not in instruction:
            wgot = False
        elif '-ggl' in instruction:
            ggl = True
            ggot = False
        glist = go_get_genenames(
            instruction[0],
            impactlist=impact_levels,
            write=wgl,
            generate_list=ggl)
        if ggot or wgot:
            ann = get_go_annotations(glist)
            summary = generate_annotation_summary(ann)
            if wgot:
                write_go_annotation_table(
                    instruction[0],
                    impactlist=impact_levels,
                    summary_df=summary)
        if ggl:
            return glist
        elif ggot:
            return summary


# WORKING WITH GENE SETS
def capture_gene_IDs_from_genesets_assign_colors(paths, cat_names, colors=[]):
    '''generates a dictionary of genes with different function according to GO from gene sets .txt files
    (can be downloaded from https://www.yeastgenome.org/go/, for example) and assigns colors to them, if specified,
    give list of paths to the .txt files, then list of cat_names, and finally a list of colors;
    see default color names, for example, in etc/colors.conf of Circos distribution'''
    geneset_dict = {}
    for i in range(len(paths)):
        gs_df = pd.read_csv(paths[i], sep='\t', skiprows=8)
        ss = pd.DataFrame(gs_df['Gene Systematic Name'])
        ss['color'] = [colors[i] for j in range(len(ss))]
        name_of_csv = str(cat_names[i]) + '.geneset'
        ss.to_csv(name_of_csv, sep="\t", index=False)
        geneset_dict[cat_names[i]] = []
        geneset_dict[cat_names[i]].append(list(ss['Gene Systematic Name']))
        geneset_dict[cat_names[i]].append(colors[i])
    for process in geneset_dict:
        geneset_dict[process][0] = list(np.unique(geneset_dict[process][0]))
    return geneset_dict


def capture_gene_coordinates_from_gff(path_to_gff, prefix='NA'):
    '''input-path to the gff file (without FASTA part), returns a dataframe with gene coordinates;
    can write it into a file if prefix is given'''
    gff = pd.read_csv(
        path_to_gff,
        sep='\t',
        skiprows=18,
        names=[
            'seqID',
            'Source',
            'SeqType',
            'Start',
            'Stop',
            'Score',
            'Strand',
            'Phase',
            'Attributes'])
    gff = gff[gff['SeqType'] == 'gene']
    gff['geneID'] = list(map(lambda x: x.split(
        ';')[0].split('=')[1], gff['Attributes']))
    gff = gff[['seqID', 'geneID', 'Start', 'Stop']]
    if prefix != 'NA':
        name_of_csv = str(prefix) + '.gene_coords'
        gff.to_csv(name_of_csv, sep="\t", index=False)
    return gff


def get_gene_subset_coords_from_gff(geneset_dict, path_to_gff):
    '''generates a txt file with the genes and their colors from the dictionary if they appear in
    the specified gff'''
    gff = capture_gene_coordinates_from_gff(path_to_gff)
    general_dict = {}
    for process in geneset_dict:
        process_dict = {
            'seqID': [],
            'start': [],
            'stop': [],
            'genes': [],
            'color': []}
        for gene in geneset_dict[process][0]:
            if gene in list(gff['geneID']):
                process_dict['genes'].append(gene)
                line = gff[gff['geneID'] == gene]
                process_dict['start'].append(int(line['Start']))
                process_dict['stop'].append(int(line['Stop']))
                process_dict['seqID'].append(str(list(line['seqID'])[0]))
                process_dict['color'].append(geneset_dict[process][1])
        process_df = pd.DataFrame(process_dict)
        general_dict[process] = process_dict
        name_of_csv = str(process) + '_in_genome.gene_coords'
        process_df.to_csv(name_of_csv, sep="\t", index=False)
    return general_dict


def get_manysamples_gene_subset_coords_from_gff(geneset_dict, path_to_gff):
    '''generates txt files with the genes and their colors from the dictionary with samples and processes
    if they appear in the specified gff'''
    gff = capture_gene_coordinates_from_gff(path_to_gff)
    general_dict = {}
    for sample in geneset_dict:
        general_dict[sample] = {}
        for process in geneset_dict[sample]:
            process_dict = {
                'seqID': [],
                'start': [],
                'stop': [],
                'genes': [],
                'color': []}
            for gene in geneset_dict[sample][process][0]:
                if gene in list(gff['geneID']):
                    process_dict['genes'].append(gene)
                    line = gff[gff['geneID'] == gene]
                    process_dict['start'].append(int(line['Start']))
                    process_dict['stop'].append(int(line['Stop']))
                    process_dict['seqID'].append(str(list(line['seqID'])[0]))
                    process_dict['color'].append(
                        geneset_dict[sample][process][1])
            process_df = pd.DataFrame(process_dict)
            general_dict[sample][process] = process_dict
            name_of_csv = str(sample) + '_' + str(process) + '.gene_coords'
            process_df.to_csv(name_of_csv, sep="\t", index=False)
    return general_dict


def generate_impact_geneset_dict(paths, prefixes):
    '''use to generate dictionaries with genes;
    paths is a list of paths to .txt files with vlists of genes;
    prefixes - is a list of prefixes for later use'''
    gene_subsets = {}
    for i in range(len(paths)):
        gene_list = pd.read_csv(paths[i], names=['geneID'])
        gene_list = list(gene_list['geneID'])
        gene_subsets[prefixes[i]] = []
        gene_subsets[prefixes[i]].append(gene_list)
    return gene_subsets


def intersect_two_genesets(gs1, gs2):
    '''generates a dictionary with genes present in both gene set dictionaries for later use
    and writes them into a .txt file,
    the gs1 is compared with gs2;
    colors are given from gs2 (which is supposed to represent processes or other ontology categories)'''
    intersection_gs = {}
    for category1 in gs1:
        intersection_gs[category1] = {}
        for gene1 in gs1[category1][0]:
            for category2 in gs2:
                for gene2 in gs2[category2][0]:
                    if gene1 == gene2:
                        if category2 not in intersection_gs[category1]:
                            intersection_gs[category1][category2] = []
                        intersection_gs[category1][category2].append(gene1)
        for process in intersection_gs[category1]:
            intersection_gs[category1][process] = [
                list(intersection_gs[category1][process])]
            intersection_gs[category1][process].append(gs2[process][1])
    return intersection_gs


def get_unique_igs(sampleID, igs):
    """sample ID - string? representing sample ID, for which unique genes must be found"""
    genes = list(np.unique(igs[sampleID][0]))
    for i in igs:
        if i != sampleID:
            genes_other = igs[i][0]
            for z in genes:
                if z in genes_other:
                    genes.remove(z)
    return genes


def get_common_igs(sampleIDs, igs):
    """sampleIDs - is a list of sample IDs, for which common genes must be found"""
    firstID = sampleIDs[0]
    common_genes = list(np.unique(igs[firstID][0]))
    for i in sampleIDs[1:]:
        gset = igs[i][0]
        common_genes = (set(common_genes) & set(gset))
    return list(common_genes)


def unique_igs(igs):
    """generates dictionaries with unique"""
    unique_igs = {}
    for i in igs.keys():
        unique_igs[i + 'u'] = []
        unique_igs[i + 'u'].append(get_unique_igs(i, igs))
    return unique_igs


def draw_gene_counts(categories_igs):
    genes_numbers_dict = {}
    genes_numbers_dict['sample_id'] = list(categories_igs.keys())
    for i in categories_igs.keys():
        for j in categories_igs[i]:
            if j not in genes_numbers_dict:
                genes_numbers_dict[j] = [
                    0 for z in range(len(categories_igs.keys()))]
    for n in range(len(categories_igs.keys())):
        key = list(categories_igs.keys())[n]
        for z in genes_numbers_dict.keys():
            if z in categories_igs[key]:
                genes_numbers_dict[z][n] = len(categories_igs[key][z][0])
    gene_counts_df = pd.DataFrame.from_dict(genes_numbers_dict)
    gene_counts_df = pd.melt(
        gene_counts_df,
        id_vars="sample_id",
        var_name="category",
        value_name="gene_counts")
    g = sns.catplot(
        x='category',
        y='gene_counts',
        hue='sample_id',
        kind='bar',
        data=gene_counts_df)
    g.set_xticklabels(rotation=30, ha="right")
    g.savefig("gene_counts.jpeg")


def capture_gene_categories():
    """General usage: "-mode filepaths=(list files with lists of genes separated by ",") annotations=(list files with lists of genes of different categories, one file for one categoty) category_names=(for each annotation file specify category) colors=(for each category specify colors) gff=(the genome's .gff file) sampleID=(sample IDs)"
First the mode must be specified: -get_igs (generate impacted gene sets) - generates .txt files with genes from the specified lists and their colors;
get_unique_igs - generates .txt files with unique impacted genes from each of the 2 input files (specify 2 input files), -draw_igs - generates a histogram showing numbers of genes from each category in each sample, -draw_unique_igs - the same, but with unique variants."""
    instruction = input('input format: mode(-gigs, -get_unique_igs, -draw_igs, -draw_unique_igs) filepaths=fpA,fpB,fpC,etc annotations=a1,a2,a3,etc category_names=cat1,cat2,cat3,etc colors=c1,c2,c3,etc gff=gff_file sID=sampleIDA,sampleIDB,sampleIDC,etc: ')
    instruction = instruction.split(' ')
    mode = instruction[0]
    filepaths = instruction[1][10:].split(',')
    annotations = instruction[2][12:].split(',')
    category_names = instruction[3][15:].split(',')
    colors = instruction[4][7:].split(',')
    gff = instruction[5][4:]
    sID = instruction[6][4:].split(',')
    colorlist = []
    for i in colors:
        colorlist.append('fill_color=' + i)
    gs_dict = capture_gene_IDs_from_genesets_assign_colors(
        annotations, category_names, colors=colorlist)
    gen_dict = get_gene_subset_coords_from_gff(gs_dict, gff)
    igs = generate_impact_geneset_dict(filepaths, sID)
    if instruction[0] == '-gigs':
        categories_igs = intersect_two_genesets(igs, gs_dict)
        get_manysamples_gene_subset_coords_from_gff(categories_igs, gff)
    elif instruction[0] == '-get_unique_igs':
        u_igs_dict = unique_igs(igs)
        categuries_u_igs = intersect_two_genesets(u_igs_dict, gs_dict)
        get_manysamples_gene_subset_coords_from_gff(categuries_u_igs, gff)
    elif instruction[0] == '-draw_igs':
        categories_igs = intersect_two_genesets(igs, gs_dict)
        draw_gene_counts(categories_igs)
    elif instruction[0] == '-draw_unique_igs':
        u_igs_dict = unique_igs(igs)
        categuries_u_igs = intersect_two_genesets(u_igs_dict, gs_dict)
        draw_gene_counts(categuries_u_igs)


# GETTING COORDINATES FOR A LIST OF GENES
def generate_impact_geneset_dict_withcolors(paths, colors, prefixes):
    '''use to generate dictionaries with genes and colors from lists of genes in .txt format;
    paths is a list of paths to .txt files with vlists of genes;
    prefixes - is a list of prefixes for later use;
    colors - list of colors;
    to use later with get_process_genes_coords'''
    gene_subsets = {}
    for i in range(len(paths)):
        gene_list = pd.read_csv(paths[i], names=['geneID'])
        gene_list = list(gene_list['geneID'])
        color = colors[i]
        gene_subsets[prefixes[i]] = []
        gene_subsets[prefixes[i]].append(gene_list)
        gene_subsets[prefixes[i]].append(colors[i])
    return gene_subsets


def get_coords():
    '''given a list of genes, generates a .txt file with their coordinates from the given gff'''
    instruction = input('input format: filepath color prefix path_to_gff: ')
    instruction = instruction.split(' ')
    path = [instruction[0]]
    color = ['fill_color=' + instruction[1]]
    prefix = [instruction[2]]
    path_to_gff = instruction[3]
    genes = generate_impact_geneset_dict_withcolors(path, color, prefix)
    get_gene_subset_coords_from_gff(genes, path_to_gff)


# DROPPING DUPLICATES
def drop_duplicates(path_to_tab, colindex, ncols):
    '''removes rows with duplicate values in a certain column,
    colnumber - 0-based, '''
    colnames = [i for i in range(int(ncols))]
    tab = pd.read_csv(path_to_tab, sep='\t', names=colnames)
    restab = tab.drop_duplicates(tab.columns[int(colindex)])
    filename = str(path_to_tab) + '_nodup.txt'
    restab.to_csv(filename, sep='\t', index=False)


def drop_duplicates_init():
    inputv = input(
        'enter path, column index (0-based), number of columns (one whitespace-separated): ').split(' ')
    drop_duplicates(inputv[0], inputv[1], inputv[2])


# GETTING COMMON GENES
def get_common_genes_init():
    inputv = input(
        'enter filepath(s) (separated by a comma) color prefix(es) (state for each sample, separated by a comma) gff: ').split(' ')
    filepaths = inputv[0].split(',')
    color = ['fill_color=' + inputv[1]]
    prefixes = inputv[2].split(',')
    gff = inputv[3]
    igs = generate_impact_geneset_dict(filepaths, prefixes)
    list_of_common_genes = get_common_igs(prefixes, igs)
    fileprefix = ''
    filename = ''
    for i in prefixes:
        filename += i + '_'
        fileprefix += i + '_'
    filename += 'common.txt'
    write_genenames(filename, list_of_common_genes)
    genes = generate_impact_geneset_dict_withcolors(
        [filename], color, [fileprefix])
    get_gene_subset_coords_from_gff(genes, gff)


# WORKING WITH SV GFF FILES
def read_gff_for_circos(path):
    data = pd.read_csv(
        path,
        sep='\t',
        skiprows=1,
        names=[
            'seqID',
            'Source',
            'Type',
            'Start',
            'End',
            'Score',
            'Strand',
            'Phase',
            'Attributes'])  # skiprows=1 identifies that the first row must be skipped
    return data


def output_vartype(data, vartype, min_length, prefix, w=1):
    '''subsets the .gff data by the type of the variant and its minimum length'''
    data_vartype = data[data['Type'] == vartype]
    data_vartype['Length'] = list(
        map(lambda x: int(x.split(';')[1].split('=')[1]), data_vartype['Attributes']))
    data_vartype = data_vartype[data_vartype['Length'] > min_length]
    data_vartype = data_vartype[['seqID', 'Start', 'End']]
    prefix_col = [prefix for i in range(len(data_vartype.seqID))]
    data_vartype['Label'] = prefix_col
    if w == 1:
        name_of_csv = str(prefix) + '_' + str(min_length) + \
            '_' + str(vartype + '.out')
        data_vartype.to_csv(name_of_csv, sep="\t")
    return data_vartype


def generate_SV_files():
    '''Specify a path first, then type of SV ('REPEAT','CNV','DEL','INS'),
        then min length and prefix'''
    instruction = input(
        'input format: path type_of_SV(REPEAT, CNV, DEL, INS) min_length_of_SV(bp) prefix: ')
    instruction = instruction.split(' ')
    path = instruction[0]
    vartype = instruction[1]
    min_length = int(instruction[2])
    prefix = instruction[3]
    data = read_gff_for_circos(path)
    output_vartype(data, vartype, min_length, prefix, w=1)


# GENERATE TOTAL SUMMARY VERSION WITH LISTS IN THE OUTPUT DATAFRAME
def total_variant_summary_init():
    '''generates a table with summary for variants with the specified impact in the genes present in the
    target_genes .txt file;
    prepare the _impact.txt files with the command: cat your_vcf.vcf | perl vcfEffOnePerLine.pl | SnpSift extractFields - CHROM POS REF ALT "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].EFFECT" "LOF[*].GENE" "LOF[*].PERC" "NMD[*].GENE" "NMD[*].PERC"> output_impact.txt'''
    instruction = input(
        'input format: paths, impacts (-b, -l, -m, -h, -lm etc.), target_genes, FileIDs, sampleIDs, seq_platforms, aligners, callers, pipelineIDs: ')
    instruction = instruction.split(' ')
    paths = instruction[0].split(',')
    impactlvl = instruction[1]
    target_genes = instruction[2]
    FileID = instruction[3].split(',')
    sampleID = instruction[4].split(',')
    platform = instruction[5].split(',')
    aligner = instruction[6].split(',')
    caller = instruction[7].split(',')
    pipelineID = instruction[8].split(',')
    global annotations_all
    annotations_all = {}
    resultant_df = generate_variant_summary(
        path=paths[0],
        impactlvl=impactlvl,
        target_genes=target_genes,
        FileID=FileID[0],
        pipelineID=pipelineID[0],
        sampleID=sampleID[0],
        platform=platform[0],
        aligner=aligner[0],
        caller=caller[0])
    for i in range(len(paths) - 1):
        sample_df = generate_variant_summary(path=paths[i + 1],
                                             impactlvl=impactlvl,
                                             target_genes=target_genes,
                                             FileID=FileID[i + 1],
                                             pipelineID=pipelineID[i + 1],
                                             sampleID=sampleID[i + 1],
                                             platform=platform[i + 1],
                                             aligner=aligner[i + 1],
                                             caller=caller[i + 1])
        resultant_df = resultant_df.append(sample_df, ignore_index=True)
    resultant_df = resultant_df[['File_ID',
                                 'Sample_ID',
                                 'Pipeline_ID',
                                 'Platform',
                                 'Aligner',
                                 'Caller',
                                 'CHROM',
                                 'POS',
                                 'REF',
                                 'ALT',
                                 'GENE_NAME',
                                 'IMPACT',
                                 'EFFECT',
                                 'LOF_PERC',
                                 'GO_ID',
                                 'GO_annotations',
                                 'Variant_ID',
                                 'Filename']]
    additional_cols = {}
    for f in FileID:
        additional_cols['Present_in_' + f] = []
    Is_unique_inpipe_col = []
    Is_unique_acrosspipes_col = []
    for row in resultant_df.itertuples():
        var = row.Variant_ID
        pipeline = row.Pipeline_ID
        varsamples = list(np.unique(
            list(resultant_df[resultant_df['Variant_ID'] == var]['File_ID'])))
        var_pipeline_samples = list(np.unique(list(resultant_df[(
            resultant_df.Variant_ID == var) & (resultant_df.Pipeline_ID == pipeline)]['Sample_ID'])))
        var_real_samples = list(np.unique(
            list(resultant_df[(resultant_df.Variant_ID == var)]['Sample_ID'])))
        varsamples.sort()
        var_pipeline_samples.sort()
        var_real_samples.sort()
        for col in additional_cols:
            if col[11:] in varsamples:
                additional_cols[col].append(1)
            else:
                additional_cols[col].append(0)
        if len(var_pipeline_samples) <= 1:
            Is_unique_inpipe_col.append(1)
        else:
            Is_unique_inpipe_col.append(0)
        if len(var_real_samples) <= 1:
            Is_unique_acrosspipes_col.append(1)
        else:
            Is_unique_acrosspipes_col.append(0)
    for newcol in additional_cols:
        resultant_df[newcol] = additional_cols[newcol]
    resultant_df['Is_unique_forthesample_within_pipeline'] = Is_unique_inpipe_col
    resultant_df['Is_unique_forthesample_acrosspipes'] = Is_unique_acrosspipes_col
    resultant_df.reset_index(inplace=True, drop=True)
    resultant_df.to_csv('_'.join(sampleID) + '.csv', index=False)
    return resultant_df


def generate_variant_summary(
        path,
        impactlvl,
        target_genes,
        FileID,
        pipelineID,
        sampleID,
        platform,
        aligner,
        caller):
    impact_levels = []
    for i in impactlvl:
        if i == 'b':
            impact_levels.append('MODIFIER')
        elif i == 'l':
            impact_levels.append('LOW')
        elif i == 'm':
            impact_levels.append('MODERATE')
        elif i == 'h':
            impact_levels.append('HIGH')
    print('impact levels: ')
    data = read_impact_data(path)
    data_subset = subset_byimpact(data, impact_levels)
    data_subset_with_genenames = add_genenames(data_subset)
    gene_list_df = pd.read_csv(target_genes, names=['geneID'])
    target_genelist = list(gene_list_df['geneID'])
    print('file: ', path)
    print('searching for the genes:', ', '.join(target_genelist))
    data_subset_with_target_genes = data_subset_with_genenames[data_subset_with_genenames['GENE_NAME'].isin(
        target_genelist)]
    data_subset_with_target_genes['Pipeline_ID'] = [
        pipelineID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['File_ID'] = [
        FileID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Sample_ID'] = [
        sampleID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Platform'] = [
        platform for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Aligner'] = [
        aligner for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Caller'] = [
        caller for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Filename'] = [
        path for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    genes = list(data_subset_with_target_genes['GENE_NAME'])
    annotations = get_go_annotations(genes)
    annotations_summary_df = generate_annotation_summary(annotations)
    annotations_summary_df_grouped_bygenes = annotations_summary_df.groupby(
        'GENE_NAME')
    annotations_col = []
    GO_ID_col = []
    for i in data_subset_with_target_genes['GENE_NAME']:
        if i in list(annotations_summary_df['GENE_NAME']):
            annotations_col.append(
                (list(
                    np.unique(
                        list(
                            annotations_summary_df_grouped_bygenes.get_group(i)['GO_TERM'])))))
            GO_ID_col.append(
                (list(
                    np.unique(
                        list(
                            annotations_summary_df_grouped_bygenes.get_group(i)['GO_ID'])))))
        else:
            annotations_col.append('NA')
            GO_ID_col.append('NA')
    data_subset_with_target_genes['GO_ID'] = GO_ID_col
    data_subset_with_target_genes['GO_ID']
    data_subset_with_target_genes['GO_annotations'] = annotations_col
    variant_ID = []
    for index, row in data_subset_with_target_genes.iterrows():
        variant_ID.append(str(row['CHROM']) +
                          str(row['POS']) + str(row['ALT']))
    data_subset_with_target_genes['Variant_ID'] = variant_ID
    return data_subset_with_target_genes


# GENERATE TOTAL SUMMARY VERSION WITHOUT LISTS IN THE OUTPUT DATAFRAME
def total_variant_summary_init_without_lists():
    '''generates a table with summary for variants with the specified impact in the genes present in the
    target_genes .txt file;
    prepare the _impact.txt files with the command: cat your_vcf.vcf | perl vcfEffOnePerLine.pl | SnpSift extractFields - CHROM POS REF ALT "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].EFFECT" "LOF[*].GENE" "LOF[*].PERC" "NMD[*].GENE" "NMD[*].PERC"> output_impact.txt'''
    instruction = input(
        'input format: paths, impacts (-b, -l, -m, -h, -lm etc.), target_genes, FileIDs, sampleIDs, seq_platforms, aligners, callers, pipelineIDs: ')
    instruction = instruction.split(' ')
    paths = instruction[0].split(',')
    impactlvl = instruction[1]
    target_genes = instruction[2]
    FileID = instruction[3].split(',')
    sampleID = instruction[4].split(',')
    platform = instruction[5].split(',')
    aligner = instruction[6].split(',')
    caller = instruction[7].split(',')
    pipelineID = instruction[8].split(',')
    global annotations_all
    annotations_all = {}
    resultant_df = generate_variant_summary_without_lists(
        path=paths[0],
        impactlvl=impactlvl,
        target_genes=target_genes,
        FileID=FileID[0],
        pipelineID=pipelineID[0],
        sampleID=sampleID[0],
        platform=platform[0],
        aligner=aligner[0],
        caller=caller[0])
    for i in range(len(paths) - 1):
        sample_df = generate_variant_summary_without_lists(path=paths[i + 1],
                                                           impactlvl=impactlvl,
                                                           target_genes=target_genes,
                                                           FileID=FileID[i + 1],
                                                           pipelineID=pipelineID[i + 1],
                                                           sampleID=sampleID[i + 1],
                                                           platform=platform[i + 1],
                                                           aligner=aligner[i + 1],
                                                           caller=caller[i + 1])
        resultant_df = resultant_df.append(sample_df, ignore_index=True)
    resultant_df = resultant_df[['File_ID',
                                 'Sample_ID',
                                 'Pipeline_ID',
                                 'Platform',
                                 'Aligner',
                                 'Caller',
                                 'CHROM',
                                 'POS',
                                 'REF',
                                 'ALT',
                                 'GENE_NAME',
                                 'IMPACT',
                                 'EFFECT',
                                 'LOF_PERC',
                                 'GO_ID',
                                 'GO_annotations',
                                 'Variant_ID',
                                 'Filename']]
    additional_cols = {}
    for f in FileID:
        additional_cols['Present_in_' + f] = []
    Is_unique_inpipe_col = []
    Is_unique_acrosspipes_col = []
    for row in resultant_df.itertuples():
        var = row.Variant_ID
        pipeline = row.Pipeline_ID
        varsamples = list(np.unique(
            list(resultant_df[resultant_df['Variant_ID'] == var]['File_ID'])))
        var_pipeline_samples = list(np.unique(list(resultant_df[(
            resultant_df.Variant_ID == var) & (resultant_df.Pipeline_ID == pipeline)]['Sample_ID'])))
        var_real_samples = list(np.unique(
            list(resultant_df[(resultant_df.Variant_ID == var)]['Sample_ID'])))
        varsamples.sort()
        var_pipeline_samples.sort()
        var_real_samples.sort()
        for col in additional_cols:
            if col[11:] in varsamples:
                additional_cols[col].append(1)
            else:
                additional_cols[col].append(0)
        if len(var_pipeline_samples) <= 1:
            Is_unique_inpipe_col.append(1)
        else:
            Is_unique_inpipe_col.append(0)
        if len(var_real_samples) <= 1:
            Is_unique_acrosspipes_col.append(1)
        else:
            Is_unique_acrosspipes_col.append(0)
    for newcol in additional_cols:
        resultant_df[newcol] = additional_cols[newcol]
    resultant_df['Is_unique_forthesample_within_pipeline'] = Is_unique_inpipe_col
    resultant_df['Is_unique_forthesample_acrosspipes'] = Is_unique_acrosspipes_col
    resultant_df.reset_index(inplace=True, drop=True)
    resultant_df.to_csv('_'.join(sampleID) + 'nolists.csv', index=False)
    return resultant_df


def generate_variant_summary_without_lists(
        path,
        impactlvl,
        target_genes,
        FileID,
        pipelineID,
        sampleID,
        platform,
        aligner,
        caller):
    impact_levels = []
    for i in impactlvl:
        if i == 'b':
            impact_levels.append('MODIFIER')
        elif i == 'l':
            impact_levels.append('LOW')
        elif i == 'm':
            impact_levels.append('MODERATE')
        elif i == 'h':
            impact_levels.append('HIGH')
    print('impact levels: ')
    data = read_impact_data(path)
    data_subset = subset_byimpact(data, impact_levels)
    data_subset_with_genenames = add_genenames(data_subset)
    gene_list_df = pd.read_csv(target_genes, names=['geneID'])
    target_genelist = list(gene_list_df['geneID'])
    print('file: ', path)
    print('searching for the genes:', ', '.join(target_genelist))
    data_subset_with_target_genes = data_subset_with_genenames[data_subset_with_genenames['GENE_NAME'].isin(
        target_genelist)]
    data_subset_with_target_genes['Pipeline_ID'] = [
        pipelineID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['File_ID'] = [
        FileID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Sample_ID'] = [
        sampleID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Platform'] = [
        platform for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Aligner'] = [
        aligner for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Caller'] = [
        caller for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Filename'] = [
        path for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    genes = list(data_subset_with_target_genes['GENE_NAME'])
    annotations = get_go_annotations(genes)
    annotations_summary_df = generate_annotation_summary(annotations)
    annotations_summary_df_grouped_bygenes = annotations_summary_df.groupby(
        'GENE_NAME')
    annotations_col = []
    GO_ID_col = []
    for i in data_subset_with_target_genes['GENE_NAME']:
        if i in list(annotations_summary_df['GENE_NAME']):
            annotations_col.append(
                ','.join(
                    list(
                        np.unique(
                            list(
                                annotations_summary_df_grouped_bygenes.get_group(i)['GO_TERM'])))))
            GO_ID_col.append(
                ','.join(
                    list(
                        np.unique(
                            list(
                                annotations_summary_df_grouped_bygenes.get_group(i)['GO_ID'])))))
        else:
            annotations_col.append('NA')
            GO_ID_col.append('NA')
    data_subset_with_target_genes['GO_ID'] = GO_ID_col
    data_subset_with_target_genes['GO_ID']
    data_subset_with_target_genes['GO_annotations'] = annotations_col
    variant_ID = []
    for index, row in data_subset_with_target_genes.iterrows():
        variant_ID.append(str(row['CHROM']) +
                          str(row['POS']) + str(row['ALT']))
    data_subset_with_target_genes['Variant_ID'] = variant_ID
    return data_subset_with_target_genes


# VERSION WITHOUT GENELIST BUT WITH ANNOTATIONS
def total_variant_summary_notarget_init():
    '''generates a table with summary for variants with the specified impact in the genes present in the
    target_genes .txt file;
    prepare the _impact.txt files with the command: cat your_vcf.vcf | perl vcfEffOnePerLine.pl | SnpSift extractFields - CHROM POS REF ALT "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].EFFECT" "LOF[*].GENE" "LOF[*].PERC" "NMD[*].GENE" "NMD[*].PERC"> output_impact.txt'''
    instruction = input(
        'input format: paths, impacts (-a, -b, -l, -m, -h, -lm etc.), FileIDs, sampleIDs, seq_platforms, aligners, callers, pipelineIDs: ')
    instruction = instruction.split(' ')
    paths = instruction[0].split(',')
    impactlvl = instruction[1]
    FileID = instruction[2].split(',')
    sampleID = instruction[3].split(',')
    platform = instruction[4].split(',')
    aligner = instruction[5].split(',')
    caller = instruction[6].split(',')
    pipelineID = instruction[7].split(',')
    global annotations_all
    annotations_all = {}
    resultant_df = generate_variant_notarget_summary(
        path=paths[0],
        impactlvl=impactlvl,
        FileID=FileID[0],
        pipelineID=pipelineID[0],
        sampleID=sampleID[0],
        platform=platform[0],
        aligner=aligner[0],
        caller=caller[0])
    for i in range(len(paths) - 1):
        sample_df = generate_variant_notarget_summary(path=paths[i + 1],
                                                      impactlvl=impactlvl,
                                                      FileID=FileID[i + 1],
                                                      pipelineID=pipelineID[i + 1],
                                                      sampleID=sampleID[i + 1],
                                                      platform=platform[i + 1],
                                                      aligner=aligner[i + 1],
                                                      caller=caller[i + 1])
        resultant_df = resultant_df.append(sample_df, ignore_index=True)
        resultant_df.to_csv('intermediate.csv', index=False)
    resultant_df = resultant_df[['File_ID',
                                 'Sample_ID',
                                 'Pipeline_ID',
                                 'Platform',
                                 'Aligner',
                                 'Caller',
                                 'CHROM',
                                 'POS',
                                 'REF',
                                 'ALT',
                                 'GENE_NAME',
                                 'IMPACT',
                                 'EFFECT',
                                 'LOF_PERC',
                                 'GO_ID',
                                 'GO_annotations',
                                 'Variant_ID',
                                 'Filename']]
    additional_cols = {}
    for f in FileID:
        additional_cols['Present_in_' + f] = []
    Is_unique_inpipe_col = []
    Is_unique_acrosspipes_col = []
    for row in resultant_df.itertuples():
        var = row.Variant_ID
        pipeline = row.Pipeline_ID
        varsamples = list(np.unique(
            list(resultant_df[resultant_df['Variant_ID'] == var]['File_ID'])))
        var_pipeline_samples = list(np.unique(list(resultant_df[(
            resultant_df.Variant_ID == var) & (resultant_df.Pipeline_ID == pipeline)]['Sample_ID'])))
        var_real_samples = list(np.unique(
            list(resultant_df[(resultant_df.Variant_ID == var)]['Sample_ID'])))
        varsamples.sort()
        var_pipeline_samples.sort()
        var_real_samples.sort()
        for col in additional_cols:
            if col[11:] in varsamples:
                additional_cols[col].append(1)
            else:
                additional_cols[col].append(0)
        if len(var_pipeline_samples) <= 1:
            Is_unique_inpipe_col.append(1)
        else:
            Is_unique_inpipe_col.append(0)
        if len(var_real_samples) <= 1:
            Is_unique_acrosspipes_col.append(1)
        else:
            Is_unique_acrosspipes_col.append(0)
    for newcol in additional_cols:
        resultant_df[newcol] = additional_cols[newcol]
    resultant_df['Is_unique_forthesample_within_pipeline'] = Is_unique_inpipe_col
    resultant_df['Is_unique_forthesample_acrosspipes'] = Is_unique_acrosspipes_col
    resultant_df.reset_index(inplace=True, drop=True)
    resultant_df.to_csv('_'.join(sampleID) + '.csv', index=False)
    return resultant_df


def generate_variant_notarget_summary(
        path,
        impactlvl,
        FileID,
        pipelineID,
        sampleID,
        platform,
        aligner,
        caller):
    impact_levels = []
    for i in impactlvl:
        if i == 'b':
            impact_levels.append('MODIFIER')
        elif i == 'l':
            impact_levels.append('LOW')
        elif i == 'm':
            impact_levels.append('MODERATE')
        elif i == 'h':
            impact_levels.append('HIGH')
    print('impact levels: ', ', '.join(impact_levels))
    data = read_impact_data(path)
    data_subset = subset_byimpact(data, impact_levels)
    data_subset_with_target_genes = add_genenames(
        data_subset)  # in this case all genes are target
    print('file: ', path)
    data_subset_with_target_genes['Pipeline_ID'] = [
        pipelineID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['File_ID'] = [
        FileID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Sample_ID'] = [
        sampleID for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Platform'] = [
        platform for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Aligner'] = [
        aligner for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Caller'] = [
        caller for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    data_subset_with_target_genes['Filename'] = [
        path for i in range(len(data_subset_with_target_genes['GENE_NAME']))]
    genes = list(data_subset_with_target_genes['GENE_NAME'])
    annotations = get_go_annotations(genes)
    annotations_summary_df = generate_annotation_summary(annotations)
    annotations_summary_df_grouped_bygenes = annotations_summary_df.groupby(
        'GENE_NAME')
    annotations_col = []
    GO_ID_col = []
    for i in data_subset_with_target_genes['GENE_NAME']:
        if i in list(annotations_summary_df['GENE_NAME']):
            annotations_col.append(
                (list(
                    np.unique(
                        list(
                            annotations_summary_df_grouped_bygenes.get_group(i)['GO_TERM'])))))
            GO_ID_col.append(
                (list(
                    np.unique(
                        list(
                            annotations_summary_df_grouped_bygenes.get_group(i)['GO_ID'])))))
        else:
            annotations_col.append('NA')
            GO_ID_col.append('NA')
    data_subset_with_target_genes['GO_ID'] = GO_ID_col
    data_subset_with_target_genes['GO_ID']
    data_subset_with_target_genes['GO_annotations'] = annotations_col
    variant_ID = []
    for index, row in data_subset_with_target_genes.iterrows():
        variant_ID.append(str(row['CHROM']) +
                          str(row['POS']) + str(row['ALT']))
    data_subset_with_target_genes['Variant_ID'] = variant_ID
    return data_subset_with_target_genes


# CONTINUING GENERAL REPORT
def continue_summary():
    instruction = input(
        'input format: paths, impacts (-a, -b, -l, -m, -h, -lm etc.), FileIDs, sampleIDs, seq_platforms, aligners, callers, pipelineIDs, current_result.csv: ')
    instruction = instruction.split(' ')
    paths = instruction[0].split(',')
    impactlvl = instruction[1]
    FileID = instruction[2].split(',')
    sampleID = instruction[3].split(',')
    platform = instruction[4].split(',')
    aligner = instruction[5].split(',')
    caller = instruction[6].split(',')
    pipelineID = instruction[7].split(',')
    current_res = instruction[8]
    resultant_df = pd.read_csv(current_res, sep=',')
    global annotations_all
    annotations_all = {}
    for i in range(len(paths)):
        sample_df = generate_variant_notarget_summary(
            path=paths[i],
            impactlvl=impactlvl,
            FileID=FileID[i],
            pipelineID=pipelineID[i],
            sampleID=sampleID[i],
            platform=platform[i],
            aligner=aligner[i],
            caller=caller[i])
        resultant_df = resultant_df.append(sample_df, ignore_index=True)
        resultant_df.to_csv('intermediate_coontinue.csv', index=False)
    resultant_df = resultant_df[['File_ID',
                                 'Sample_ID',
                                 'Pipeline_ID',
                                 'Platform',
                                 'Aligner',
                                 'Caller',
                                 'CHROM',
                                 'POS',
                                 'REF',
                                 'ALT',
                                 'GENE_NAME',
                                 'IMPACT',
                                 'EFFECT',
                                 'LOF_PERC',
                                 'GO_ID',
                                 'GO_annotations',
                                 'Variant_ID',
                                 'Filename']]
    additional_cols = {}
    for f in FileID:
        additional_cols['Present_in_' + f] = []
    Is_unique_inpipe_col = []
    Is_unique_acrosspipes_col = []
    for row in resultant_df.itertuples():
        var = row.Variant_ID
        pipeline = row.Pipeline_ID
        varsamples = list(np.unique(
            list(resultant_df[resultant_df['Variant_ID'] == var]['File_ID'])))
        var_pipeline_samples = list(np.unique(list(resultant_df[(
            resultant_df.Variant_ID == var) & (resultant_df.Pipeline_ID == pipeline)]['Sample_ID'])))
        var_real_samples = list(np.unique(
            list(resultant_df[(resultant_df.Variant_ID == var)]['Sample_ID'])))
        varsamples.sort()
        var_pipeline_samples.sort()
        var_real_samples.sort()
        for col in additional_cols:
            if col[11:] in varsamples:
                additional_cols[col].append(1)
            else:
                additional_cols[col].append(0)
        if len(var_pipeline_samples) <= 1:
            Is_unique_inpipe_col.append(1)
        else:
            Is_unique_inpipe_col.append(0)
        if len(var_real_samples) <= 1:
            Is_unique_acrosspipes_col.append(1)
        else:
            Is_unique_acrosspipes_col.append(0)
    for newcol in additional_cols:
        resultant_df[newcol] = additional_cols[newcol]
    resultant_df['Is_unique_forthesample_within_pipeline'] = Is_unique_inpipe_col
    resultant_df['Is_unique_forthesample_acrosspipes'] = Is_unique_acrosspipes_col
    resultant_df.reset_index(inplace=True, drop=True)
    resultant_df.to_csv('_'.join(sampleID) + '.csv', index=False)
    return resultant_df


# GENERAL INIT
def general_launch():
    instruction = input(
        'input mode: get_go_annotation or geneset_job or get_coords or generate_SV_files or drop_dupl or get_common_genes or total_variant_summary_fortargets or total_variant_summary_notargets or continue_summary: ')
    instruction = instruction.split(' ')
    if instruction[0] == 'get_go_annotation':
        get_go_annotations_init()
    elif instruction[0] == 'geneset_job':
        capture_gene_categories()
    elif instruction[0] == 'generate_SV_files':
        generate_SV_files()
    elif instruction[0] == 'get_coords':
        get_coords()
    elif instruction[0] == 'drop_dupl':
        drop_duplicates_init()
    elif instruction[0] == 'get_common_genes':
        get_common_genes_init()
    elif instruction[0] == 'total_variant_summary_fortargets':
        total_variant_summary_init()
    elif instruction[0] == 'total_variant_summary_notargets':
        total_variant_summary_notarget_init()
    elif insruction[0] == 'continue_summary':
        continue_summary()


general_launch()
