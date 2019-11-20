#!/usr/bin/env python
import os
import argparse
import pandas as pd

# set_up the parser for input
parser = argparse.ArgumentParser('filters annovar output with japan params')
parser.add_argument('-TM2_limit', type=int, help='test value for argparse')
parser.add_argument('input', type=str, help='input from annovar')
parser.add_argument('output', type=str, help='output file to filtered/..')

# read arguments
args = parser.parse_args()
i, o = args.input, args.output
output_base = os.path.splitext(o)[0]

# params
TM2_limit = args.TM2_limit

# load the file as a pandas Dataframe
df = pd.read_csv(i, sep='\t', low_memory=False, na_values='.', dtype={'Chr':str, 'Start':int, 'End':int}).sort_values(['Chr', 'Start'])


def daniel_filter(data, output):
    '''
    creates filtered output using the daniel filtering
    input: pd.dataframe of unfiltered annovar output
    output: 
    - filtered/sample_tumor-normal_daniel.csv
    '''
    ############## FILTERS ##########################################

    filter_stats = {}
    
    filter_stats['initial'] = len(data.index)


    #filter exonic, splicing and exonic;splicing (but not ncRNA_splicing, _exonic,...)
    exon = data['func'] == 'exonic'
    splicing = data['func'] == 'splicing'
    ex_spl = data['func'] == 'exonic;splicing'
    exon_splicing = exon | splicing | ex_spl
    exonic_data = data[exon_splicing]
    filter_stats.update({'exonic':len(exonic_data.index)})


    #take out reads with tumor_reads2 < 4
    TM2 = exonic_data['TR2'] >= TM2_limit
    TM2plus = exonic_data['TR2+'] > 0
    TM2minus = exonic_data['TR2-'] > 0
    TM2_data = exonic_data[TM2 & TM2plus & TM2minus]
    filter_stats.update({'TM2':len(TM2_data.index)})


    # take out calls that have an entry for GenomicSuperDups-DB
    no_superD = TM2_data['superDups'].isnull()
    no_dups = TM2_data[no_superD]
    filter_stats.update({'superDups_filtered':len(no_dups.index)})


    #take out calls that have an entry for dbsnp and not for cosmic
    cosmic70 = no_dups['cosmic70'].notnull()
    no_dbsnp = no_dups['snp138'].isnull()
    db_filter = cosmic70 | no_dbsnp
    db_filtered = no_dups[db_filter]
    filter_stats.update({'dbsnp_filtered':len(db_filtered.index)})

    db_filtered.to_csv(f"{output}.csv", index=False)


    #filter Somatic calls and remove synonymous variants
    somatic = db_filtered['somatic_status'] == 'Somatic'
    not_synonymous = db_filtered['exonic_func'] != 'synonymous SNV'
    sn = somatic & not_synonymous
    somatic_data = db_filtered[sn]
    filter_stats.update({'somatic':len(somatic_data.index)})

    #filter Germline calls and remove synonymous variants
    germline = db_filtered['somatic_status'] == 'Germline'
    not_synonymous = db_filtered['exonic_func'] != 'synonymous SNV'
    gn = germline & not_synonymous
    germline_data = db_filtered[gn]
    filter_stats.update({'germline':len(somatic_data.index)})

    #filter LOH calls and remove synonymous variants
    loh = db_filtered['somatic_status'] == 'LOH'
    not_synonymous = db_filtered['exonic_func'] != 'synonymous SNV'
    ln = germline & not_synonymous
    loh_data = db_filtered[ln]
    filter_stats.update({'LOH':len(somatic_data.index)})

# write selected output to files
    somatic_data.to_csv(f"{output}_somatic.csv", index=False)
    germline_data.to_csv(f"{output}_germline.csv", index=False)
    loh_data.to_csv(f"{output}_loh.csv", index=False)

daniel_filter(df, output_base)