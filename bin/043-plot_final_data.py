#!/usr/bin/env python
__version__ = '0.0.3'

import argparse
import os
import scanpy as sc
import plotnine as plt9
# import pandas as pd


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Plots number of cells per sample in AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5ad', '--h5ad_file',
        action='store',
        dest='h5ad',
        required=True,
        help='H5AD file.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='output_file',
        default='',
        help='Basename of output png file. Will have .png appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    anndata_file = options.h5ad
    out_file_base = options.output_file

    # Get basename of the output file
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(anndata_file.rstrip('h5ad').rstrip('\\.'))
        )

    # Load the AnnData file and set matrix to be log data
    adata = sc.read_h5ad(filename=anndata_file)

    # One could also do the below... cell_passes_qc is set during the
    # merge in the CDI pipeline based on param filters which should equal
    # the same thing above
    filt = (adata.obs['cell_passes_qc'] == True)
    adata = adata[filt, :]

    # Get the total number of cells per sample
    df_ncells = adata.obs['experiment_id'].value_counts()
    df_ncells = df_ncells.rename_axis(
        'experiment_id'
    ).reset_index(name='n_cells')
    df_ncells['dummy_value'] = ''

    # Plot the number of cells for each sample
    gplt = plt9.ggplot(df_ncells, plt9.aes(
        x='dummy_value',
        y='n_cells'
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_boxplot(outlier_alpha=0)
    gplt = gplt + plt9.geom_jitter(
        # plt9.aes(color='cell_filter_type'),
        alpha=0.15,
        width=0.25
    )
    # gplt = gplt + plt9.geom_text(vjust=1.6, color='white', size=3.5)
    gplt = gplt + plt9.scale_y_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.scale_fill_brewer(
        palette='Dark2',
        type='qual'
    )
    gplt = gplt + plt9.labs(
        title='',
        y='Number of cells',
        x='',
        fill=''
    )
    gplt.save(
        '{}-n_cells_boxplot.png'.format(out_file_base),
        dpi=300/2,
        width=4,
        height=5
    )

    # Another boxplot faceted by disease status
    if 'disease_status' in adata.obs.columns:
        # Get the total number of cells per sample faceted by disease_status
        df_ncells = adata.obs.groupby(['experiment_id'])[
            'disease_status'
        ].value_counts().reset_index(name='n_cells')

        # Plot the number of cells for each sample
        gplt = plt9.ggplot(df_ncells, plt9.aes(
            x='disease_status',
            y='n_cells'
        ))
        gplt = gplt + plt9.theme_bw()
        gplt = gplt + plt9.geom_boxplot(outlier_alpha=0)
        gplt = gplt + plt9.geom_jitter(
            # plt9.aes(color='cell_filter_type'),
            alpha=0.15,
            width=0.25
        )
        # gplt = gplt + plt9.geom_text(vjust=1.6, color='white', size=3.5)
        gplt = gplt + plt9.scale_y_continuous(
            trans='log10',
            labels=comma_labels,
            minor_breaks=0
        )
        gplt = gplt + plt9.scale_fill_brewer(
            palette='Dark2',
            type='qual'
        )
        gplt = gplt + plt9.labs(
            title='',
            y='Number of cells',
            x='',
            fill=''
        )
        gplt.save(
            '{}-n_cells_boxplot-disease_status.png'.format(out_file_base),
            dpi=300/2,
            width=4,
            height=5
        )

    # Now generate the cell type counts matrix
    # Make long dataframe of experiment_id cluster number_of_cells
    df_cell_counts = adata.obs.groupby(
        ['experiment_id', 'predicted_celltype']
    ).size().reset_index(name='n_cells')
    # df_cell_counts = df_cell_counts.pivot(
    #     index='experiment_id',
    #     columns='predicted_celltype',
    #     values='n_cells'
    # )
    df_cell_counts['less_than_five_cells'] = df_cell_counts['n_cells'] < 5

    gplt = plt9.ggplot(df_cell_counts, plt9.aes(
        x='experiment_id',
        y='predicted_celltype'
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=90))
    gplt = gplt + plt9.scale_color_brewer(
        palette='Dark2',
        type='qual'
    )
    gplta = gplt + plt9.geom_point(
        plt9.aes(
            # alpha='n_cells',
            color='less_than_five_cells',
            size='n_cells'
        ),
        alpha=0.75
    )
    gplta = gplta + plt9.labs(
        title='',
        x='Sample',
        y='Cell types',
        color='Less than 5 cells',
        size='Number of cells',
        alpha='Number of cells'
    )
    gplta.save(
        '{}-n_cells_dotplot-celltype.png'.format(out_file_base),
        width=30,
        height=10,
        limitsize=False
    )


if __name__ == '__main__':
    main()
