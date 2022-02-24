#!/usr/bin/env python

__author__ = 'Leland Taylor'
__date__ = '2020-07-09'
__version__ = '0.0.1'

# import os
import argparse
import pandas as pd


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Format a file to
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '--cb_results_tsvs',
        action='store',
        dest='cb_results_tsvs',
        required=True,
        help='Comma-separated list of result tsv files with paths to\
            the final cellbender 10x results.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='output_file',
        default='file_paths_10x',
        help='Output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Read in the list of cellranger results files per sample...
    # These are actually file that contain the full paths of the CellBender
    # results
    res_files = options.cb_results_tsvs.split(',')

    # Join everything together into a dataframe
    res_pd_list = []
    for i in res_files:
        res_pd_list.append(pd.read_csv(i, sep='\t'))
    df = pd.concat(res_pd_list, axis=0, sort=False).reset_index(drop=True)

    # Write a different output file for each FPR and parameter sets
    col_groups = ['fpr', 'cb_params']
    df_grouped = df.groupby(col_groups)
    for name, group in df_grouped:
        # If an experiment_id occurs >1 in the data then a problem has happened
        if group['experiment_id'].value_counts().max() > 1:
            raise Exception('ERROR: duplicate experiment_id.')
        group = group.drop(col_groups, axis=1)

        # Now write the list of files
        group.to_csv(
            '{}-{}-FPR_{}.tsv'.format(
                options.output_file,
                name[1].replace('-', '__'),
                str(name[0]).replace('.', 'pt')
            ),
            sep='\t',
            index=False
        )


if __name__ == '__main__':
    main()
