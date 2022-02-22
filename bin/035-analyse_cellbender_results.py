#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-10-27'
__version__ = '0.0.1'


import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import argparse
import tables
import pandas as pd
import numpy as np
import scipy.sparse
from typing import Dict
import scanpy as sc
import anndata
import plotnine as plt9
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec

# To resolve strange TclError for interactive job
import matplotlib
matplotlib.use('Agg')  # Agg for png and pdf for pdf

# Nice large palette.
COLORS_LARGE_PALLETE = [
    '#0F4A9C', '#3F84AA', '#C9EBFB', '#8DB5CE', '#C594BF', '#DFCDE4',
    '#B51D8D', '#6f347a', '#683612', '#B3793B', '#357A6F', '#989898',
    '#CE778D', '#7F6874', '#E09D37', '#FACB12', '#2B6823', '#A0CC47',
    '#77783C', '#EF4E22', '#AF1F26'
]


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(
    file: str,
    analyzed_barcodes_only: bool = True
) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count
            matrix. True to load a limited set of barcodes: only those
            analyzed by the algorithm. This allows relevant latent
            variables to be loaded properly into adata.obs and adata.obsm,
            rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """
    d = dict_from_h5(file)
    X = scipy.sparse.csc_matrix(
        (d.pop('data'), d.pop('indices'), d.pop('indptr')),
        shape=d.pop('shape')
    ).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print(
                'Warning: analyzed_barcodes_only=True, but the key ',
                '"barcodes_analyzed_inds" or "barcode_indices_for_latents" ',
                'is missing from the h5 file. ',
                'Will output all barcodes, and proceed as if ',
                'analyzed_barcodes_only=False'
            )

    print(d.keys())

    # Construct the count matrix.
    if 'gene_names' in d.keys():
        gene_symbols = d.pop('gene_names').astype(str)
    else:
        gene_symbols = d.pop('name').astype(str)
    adata = anndata.AnnData(
        X=X,
        obs={'barcode': d.pop('barcodes').astype(str)},
        var={
            'gene_ids': d.pop('id').astype(str),
            'gene_symbols': gene_symbols
        }
    )
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_ids', inplace=True)

    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print(
                'Unable to load data into AnnData: ', key, value, type(value)
            )

    if analyzed_barcodes_only:
        cols = adata.obs.columns[
            adata.obs.columns.str.startswith('barcodes_analyzed')
            | adata.obs.columns.str.startswith('barcode_indices')
        ]
        for col in cols:
            try:
                del adata.obs[col]
            except Exception:
                pass

    return adata


def load_anndata_from_input_and_output(
    input_dir_10x: str,
    cellbender_h5: str,
    analyzed_barcodes_only: bool = True
) -> 'anndata.AnnData':
    """Load remove-background output count matrix into an anndata object.

    Args:
        input_dir: Raw 10x dir.
        output_file: Output h5 file created by remove-background (can be
            filtered or not).
        analyzed_barcodes_only: Argument passed to anndata_from_h5().
            False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count
            matrix. True to load a limited set of barcodes: only those
            analyzed by the algorithm. This allows relevant latent variables
            to be loaded properly into adata.obs and adata.obsm, rather than
            adata.uns.

    Return:
        adata_out: AnnData object with counts before and after
            remove-background, as well as inferred latent variables from
            remove-background.
    """
    # Load input data.
    adata_10x = sc.read_10x_mtx(
        path=input_dir_10x,
        # var_names='gene_symbols',
        var_names='gene_ids',
        make_unique=False
    )

    # Load remove-background output data.
    # We need to do this because of bug here:
    # https://github.com/broadinstitute/CellBender/issues/57
    adata_out = anndata_from_h5(
        cellbender_h5,
        analyzed_barcodes_only=analyzed_barcodes_only
    )

    # Subset the raw dataset to the relevant barcodes.
    adata_10x = adata_10x[adata_out.obs.index]

    # Put count matrices into 'layers' in anndata for clarity.
    adata_out.layers['counts_raw'] = adata_10x.X.copy()
    adata_out.layers['counts_cellbender'] = adata_out.X.copy()

    # Pre-compute a bit of metadata.
    # adata_out.var['n_cellranger'] = np.array(
    #     adata_out.layers['cellranger'].sum(axis=0)
    # ).squeeze()
    # adata_out.var['n_cellbender'] = np.array(
    #     adata_out.layers['cellbender'].sum(axis=0)
    # ).squeeze()

    return adata_out


def _make_data_plot_difference(adata, X, gene_list):
    filt = np.isin(
        adata.var.gene_symbols.values,
        gene_list
    )
    df_plt = pd.DataFrame(data=X[:, filt].todense())
    df_plt.index = adata.obs.index
    df_plt.columns = adata.var.gene_symbols.values[filt]
    df_plt = df_plt[gene_list]
    df_plt = pd.melt(df_plt, ignore_index=False)
    df_plt['gene_symbols'] = pd.Categorical(
        df_plt['variable'],
        categories=gene_list
    )
    return df_plt


def plot_difference(
    adata,
    plot_name='cellbender_results'
):
    # Get the differences in counts per cell
    X_raw_minus_cb = adata.layers[
        'counts_raw'
    ] - adata.layers['counts_cellbender']
    X_dif = abs(X_raw_minus_cb)

    # Get the top most different genes
    df_diff_genes = pd.DataFrame(data=adata.var.gene_symbols.values)
    df_diff_genes['ensembl_id'] = adata.var.index
    df_diff_genes['gene_symbols'] = adata.var.gene_symbols.values
    df_diff_genes['dif_across_cells'] = np.asarray(
        X_dif.sum(axis=0)
    ).reshape(-1)
    df_diff_genes = df_diff_genes.sort_values(
        'dif_across_cells',
        ascending=False
    )

    # Select the top 100 genes and plot the difference in counts across
    # cells where x axis = gene, y axis = difference, and point = cell.
    top_n_genes = 100
    df_plt = _make_data_plot_difference(
        adata,
        X_raw_minus_cb,
        df_diff_genes['gene_symbols'].head(n=top_n_genes)
    )
    # print(df_plt.head())
    gplt = plt9.ggplot(df_plt)
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_boxplot(
        plt9.aes(x='gene_symbols', y='value'),
        alpha=0.25
        #outlier_shape=''
    )
    gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=90))
    gplt = gplt + plt9.labs(
        x='',
        y='Raw counts - cellbender adjusted',
        title='Top {} most different genes'.format(top_n_genes)
    )
    gplt.save(
        '{}-count_difference-boxplot.png'.format(plot_name),
        #dpi=300,
        width=14,
        height=4
    )
    gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=0))
    gplt = gplt + plt9.coord_flip()
    gplt.save(
        '{}-count_difference-boxplot_vertical.png'.format(plot_name),
        #dpi=300,
        width=4,
        height=14
    )

    # Same plot but the abs difference on log scale
    df_plt = _make_data_plot_difference(
        adata,
        X_dif,
        df_diff_genes['gene_symbols'].head(n=top_n_genes)
    )
    df_plt['value'] += 1
    # print(df_plt.head())
    gplt = plt9.ggplot(df_plt)
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_boxplot(
        plt9.aes(x='gene_symbols', y='value'),
        alpha=0.25
        # outlier_shape=''
    )
    gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=90))
    gplt = gplt + plt9.labs(
        x='',
        y='Abs(raw counts - cellbender adjusted)',
        title='Top {} most different genes'.format(top_n_genes)
    )
    gplt = gplt + plt9.scale_y_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    gplt.save(
        '{}-abs_count_difference-boxplot.png'.format(plot_name),
        #dpi=300,
        width=14,
        height=4
    )
    gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=0))
    gplt = gplt + plt9.coord_flip()
    gplt.save(
        '{}-abs_count_difference-boxplot_vertical.png'.format(plot_name),
        #dpi=300,
        width=4,
        height=14
    )


def plot_ambient_by_difference(
    adata,
    plot_name='cellbender_results'
):

    # Compute the total amount of expression of each gene
    adata.var['total_gene_counts_raw'] = np.array(
        adata.layers['counts_raw'].sum(axis=0)
    ).squeeze()
    adata.var['total_gene_counts_cellbender'] = np.array(
        adata.layers['counts_cellbender'].sum(axis=0)
    ).squeeze()

    adata.var['difference_total_gene_counts_raw_cellbender'] = adata.var[
        'total_gene_counts_raw'
    ] - adata.var['total_gene_counts_cellbender']

    # Make the plot
    gplt = plt9.ggplot(adata.var)
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_point(
        plt9.aes(
            x='ambient_expression',
            y='difference_total_gene_counts_raw_cellbender'
        ),
        alpha=0.25
    )
    gplt = gplt + plt9.labs(
        x='Ambient RNA signature',
        y='Counts removed by cellbender',
        title='Ambient RNA signature removal per gene'
    )
    # gplt = gplt + plt9.scale_y_continuous(
    #     trans='log10',
    #     labels=comma_labels,
    #     minor_breaks=0
    # )
    gplt.save(
        '{}-ambient_signature-scatter.png'.format(plot_name),
        #dpi=300,
        width=5,
        height=5
    )

    # Add gene names to the plot
    gplt = plt9.ggplot(adata.var)
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_text(
        plt9.aes(
            x='ambient_expression',
            y='difference_total_gene_counts_raw_cellbender',
            label='gene_symbols'
        ),
        alpha=0.25
    )
    gplt = gplt + plt9.labs(
        x='Ambient RNA signature',
        y='Counts removed by cellbender',
        title='Ambient RNA signature removal per gene'
    )
    # gplt = gplt + plt9.scale_y_continuous(
    #     trans='log10',
    #     labels=comma_labels,
    #     minor_breaks=0
    # )
    gplt.save(
        '{}-ambient_signature-scatter_genenames.png'.format(plot_name),
        #dpi=300,
        width=5,
        height=5
    )


def basic_cluster(adata, count_layer):
    # set the right counts
    adata.X = adata.layers[count_layer]

    # Label mitochondrial encoded trancripts
    # This includes:
    # * Mitochondrially encoded protein coding genes
    # * Mitochondrially encoded transcribed regions
    # * Mitochondrially encoded RNAs
    #
    # The below gene list was downloaded from on 3 Aug 2020:
    # https://www.genenames.org/data/genegroup/#!/group/1972
    gene_group__mito_transcript = [
        'MT-7SDNA',
        'MT-ATP6',
        'MT-ATP8',
        'MT-ATT',
        'MT-CO1',
        'MT-CO2',
        'MT-CO3',
        'MT-CSB1',
        'MT-CSB2',
        'MT-CSB3',
        'MT-CYB',
        'MT-HPR',
        'MT-HSP1',
        'MT-HSP2',
        'MT-LIPCAR',
        'MT-LSP',
        'MT-ND1',
        'MT-ND2',
        'MT-ND3',
        'MT-ND4',
        'MT-ND4L',
        'MT-ND5',
        'MT-ND6',
        'MT-OHR',
        'MT-OLR',
        'MT-RNR1',
        'MT-RNR2',
        'MT-RNR3',
        'MT-TA',
        'MT-TAS',
        'MT-TC',
        'MT-TD',
        'MT-TE',
        'MT-TER',
        'MT-TF',
        'MT-TFH',
        'MT-TFL',
        'MT-TFX',
        'MT-TFY',
        'MT-TG',
        'MT-TH',
        'MT-TI',
        'MT-TK',
        'MT-TL1',
        'MT-TL2',
        'MT-TM',
        'MT-TN',
        'MT-TP',
        'MT-TQ',
        'MT-TR',
        'MT-TS1',
        'MT-TS2',
        'MT-TT',
        'MT-TV',
        'MT-TW',
        'MT-TY'
    ]
    adata.var['gene_group__mito_transcript'] = [
        x in gene_group__mito_transcript for x in adata.var['gene_symbols']
    ]
    # Label ribosomal protein genes
    # Ribosomal protein: A ribosomal protein is any of
    # the proteins that, in conjunction with rRNA, make up the ribosomal
    # subunits involved in the cellular process of translation. A large
    # part of the knowledge about these organic molecules has come from the
    # study of E. coli ribosomes. Most ribosomal proteins have been
    # isolated and specific antibodies have been produced. These, together
    # with electronic microscopy and the use of certain reactives, have
    # allowed for the determination of the topography of the proteins in
    # ribosome. E.coli, other bacteria and Archaea have a 30S small subunit
    # and a 50S large subunit, whereas humans and yeasts have a 40S small
    # subunit and a 60S large subunit. Equivalent subunits are frequently
    # numbered differently between bacteria, Archaea, yeasts and humans.
    #
    # The below gene list was downloaded from on 3 Aug 2020:
    # https://www.genenames.org/data/genegroup/#!/group/1054
    gene_group__ribo_protein = set([
        'DAP3',
        'FAU',
        'MRPL1',
        'MRPL1',
        'MRPL10',
        'MRPL10',
        'MRPL11',
        'MRPL11',
        'MRPL12',
        'MRPL12',
        'MRPL13',
        'MRPL13',
        'MRPL14',
        'MRPL14',
        'MRPL15',
        'MRPL15',
        'MRPL16',
        'MRPL16',
        'MRPL17',
        'MRPL17',
        'MRPL18',
        'MRPL18',
        'MRPL19',
        'MRPL19',
        'MRPL2',
        'MRPL2',
        'MRPL20',
        'MRPL20',
        'MRPL21',
        'MRPL21',
        'MRPL22',
        'MRPL22',
        'MRPL23',
        'MRPL23',
        'MRPL24',
        'MRPL24',
        'MRPL27',
        'MRPL27',
        'MRPL28',
        'MRPL28',
        'MRPL3',
        'MRPL3',
        'MRPL30',
        'MRPL30',
        'MRPL32',
        'MRPL32',
        'MRPL33',
        'MRPL33',
        'MRPL34',
        'MRPL34',
        'MRPL35',
        'MRPL35',
        'MRPL36',
        'MRPL36',
        'MRPL37',
        'MRPL37',
        'MRPL38',
        'MRPL38',
        'MRPL39',
        'MRPL39',
        'MRPL4',
        'MRPL4',
        'MRPL40',
        'MRPL40',
        'MRPL41',
        'MRPL41',
        'MRPL42',
        'MRPL42',
        'MRPL43',
        'MRPL43',
        'MRPL44',
        'MRPL44',
        'MRPL45',
        'MRPL45',
        'MRPL46',
        'MRPL46',
        'MRPL47',
        'MRPL47',
        'MRPL48',
        'MRPL48',
        'MRPL49',
        'MRPL49',
        'MRPL50',
        'MRPL50',
        'MRPL51',
        'MRPL51',
        'MRPL52',
        'MRPL52',
        'MRPL53',
        'MRPL53',
        'MRPL54',
        'MRPL54',
        'MRPL55',
        'MRPL55',
        'MRPL57',
        'MRPL57',
        'MRPL58',
        'MRPL9',
        'MRPS10',
        'MRPS10',
        'MRPS11',
        'MRPS11',
        'MRPS12',
        'MRPS12',
        'MRPS14',
        'MRPS14',
        'MRPS15',
        'MRPS15',
        'MRPS16',
        'MRPS16',
        'MRPS17',
        'MRPS17',
        'MRPS18A',
        'MRPS18A',
        'MRPS18B',
        'MRPS18B',
        'MRPS18C',
        'MRPS18C',
        'MRPS2',
        'MRPS2',
        'MRPS21',
        'MRPS21',
        'MRPS22',
        'MRPS22',
        'MRPS23',
        'MRPS23',
        'MRPS24',
        'MRPS24',
        'MRPS25',
        'MRPS25',
        'MRPS26',
        'MRPS26',
        'MRPS27',
        'MRPS27',
        'MRPS28',
        'MRPS28',
        'MRPS30',
        'MRPS30',
        'MRPS31',
        'MRPS31',
        'MRPS33',
        'MRPS33',
        'MRPS34',
        'MRPS34',
        'MRPS35',
        'MRPS35',
        'MRPS36',
        'MRPS36',
        'MRPS5',
        'MRPS6',
        'MRPS7',
        'MRPS9',
        'RPL10',
        'RPL10A',
        'RPL10L',
        'RPL11',
        'RPL12',
        'RPL13A',
        'RPL14',
        'RPL15',
        'RPL17',
        'RPL18A',
        'RPL19',
        'RPL21',
        'RPL22',
        'RPL23',
        'RPL23A',
        'RPL24',
        'RPL26',
        'RPL26L1',
        'RPL27',
        'RPL27A',
        'RPL28',
        'RPL29',
        'RPL3',
        'RPL30',
        'RPL31',
        'RPL32',
        'RPL34',
        'RPL35',
        'RPL35A',
        'RPL36',
        'RPL36A',
        'RPL36AL',
        'RPL37',
        'RPL37A',
        'RPL38',
        'RPL39',
        'RPL39L',
        'RPL3L',
        'RPL4',
        'RPL41',
        'RPL5',
        'RPL6',
        'RPL7',
        'RPL7A',
        'RPL7L1',
        'RPL8',
        'RPL9',
        'RPLP0',
        'RPLP1',
        'RPLP2',
        'RPS10',
        'RPS11',
        'RPS12',
        'RPS13',
        'RPS14',
        'RPS15',
        'RPS15A',
        'RPS16',
        'RPS17',
        'RPS18',
        'RPS19',
        'RPS2',
        'RPS20',
        'RPS21',
        'RPS23',
        'RPS24',
        'RPS25',
        'RPS26',
        'RPS27',
        'RPS27A',
        'RPS27L',
        'RPS28',
        'RPS29',
        'RPS3',
        'RPS3A',
        'RPS4X',
        'RPS4Y1',
        'RPS4Y2',
        'RPS5',
        'RPS6',
        'RPS7',
        'RPS8',
        'RPS9',
        'UBA52'
    ])
    adata.var['gene_group__ribo_protein'] = [
        x in gene_group__ribo_protein for x in adata.var['gene_symbols']
    ]

    # calculate basic qc metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[
            'gene_group__mito_transcript',
            'gene_group__ribo_protein'
        ],
        inplace=True
    )

    # filter low quality cells
    filter_queries = [
        'n_genes_by_counts < 100',
        # 'pct_counts_gene_group__ribo_protein > 50',
        'pct_counts_gene_group__mito_transcript > 25'
    ]
    adata.obs['cell_passes_qc'] = True
    for filter_query in filter_queries:
        cells_to_remove = adata.obs.query(filter_query).index
        adata.obs.loc[cells_to_remove, 'cell_passes_qc'] = False
        adata = adata[
            np.invert(adata.obs.index.isin(cells_to_remove)),
            :
        ]

    # Add a raw counts layer.
    adata.layers['counts_analyzed'] = adata.X.copy()

    # normalise
    sc.pp.filter_genes(adata, min_cells=5)
    # Total-count normalize (library-size correct) the data matrix X to
    # counts per million, so that counts become comparable among cells.
    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
        exclude_highly_expressed=False,
        key_added='normalization_factor',  # add to adata.obs
        inplace=True
    )
    # Logarithmize the data: X = log(X + 1) where log = natural logorithm.
    # Numpy has a nice function to undo this np.expm1(adata.X).
    sc.pp.log1p(adata)

    # calculate HVG
    # Calculate the highly variable genes on the log1p(norm) data.
    # Do so for each sample and then merge - this avoids the selection of
    # batch-specific, highly variable genes.
    sc.pp.highly_variable_genes(
        adata,
        # min_mean=0.0125,
        # max_mean=3,
        # min_disp=0.5,
        flavor='seurat',
        n_top_genes=2000,  # 2000 = SeuratFindVariableFeatures
        inplace=True
    )

    # Exclude mitocondrial and ribosomal genes from highly variable gene set.
    adata.var.loc[
        adata.var['gene_group__mito_transcript'],
        ['highly_variable']
    ] = False
    adata.var.loc[
        adata.var['gene_group__ribo_protein'],
        ['highly_variable']
    ] = False

    # Scale the data to unit variance.
    # This effectively weights each gene evenly. Otherwise
    # genes with higher expression values will naturally have higher
    # variation that will be captured by PCA
    sc.pp.scale(
        adata,
        zero_center=False,  # If true, sparse becomes dense
        max_value=None,
        copy=False
    )

    # calculate PCs
    n_pcs = 15
    sc.tl.pca(
        adata,
        n_comps=n_pcs,
        zero_center=True,  # Set to true for standard PCA
        svd_solver='arpack',  # arpack reproducible when zero_center = True
        use_highly_variable=True,
        copy=False,
        random_state=np.random.RandomState(0),
        chunked=False
    )

    # calculate neighbors
    sc.pp.neighbors(
        adata,
        use_rep='X_pca',
        n_pcs=n_pcs,
        copy=False,
        random_state=0
    )

    # calculate umap
    sc.tl.umap(
        adata,
        n_components=2,
        copy=False,
        random_state=0
    )

    # cluster the data
    cluster_method = 'leiden'
    if cluster_method == 'leiden':
        sc.tl.leiden(
            adata,
            copy=False,
            random_state=0
        )
    elif cluster_method == 'louvain':
        sc.tl.louvain(
            adata,
            flavor='vtraag',
            copy=False,
            random_state=0
        )

    # Also save the clusters to the same spot so we know where they will be.
    adata.uns['cluster'] = adata.uns[cluster_method]
    adata.uns['cluster']['params']['method'] = cluster_method
    adata.obs['cluster'] = adata.obs[cluster_method]

    return adata


# This function is based off of scanpy:
# https://github.com/theislab/scanpy/blob/master/scanpy/plotting/_tools/scatterplots.py
def panel_grid(hspace, wspace, ncols, num_panels):
    """Init plot."""
    n_panels_x = min(ncols, num_panels)
    n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02
    # each panel will have the size of rcParams['figure.figsize']
    fig = plt.figure(
        figsize=(
            n_panels_x * rcParams['figure.figsize'][0] * (1 + wspace),
            n_panels_y * rcParams['figure.figsize'][1],
        )
    )
    left = 0.2 / n_panels_x
    bottom = 0.13 / n_panels_y
    gs = gridspec.GridSpec(
        nrows=n_panels_y,
        ncols=n_panels_x,
        left=left,
        right=1 - (n_panels_x - 1) * left - 0.01 / n_panels_x,
        bottom=bottom,
        top=1 - (n_panels_y - 1) * bottom - 0.1 / n_panels_y,
        hspace=hspace,
        wspace=wspace
    )
    return fig, gs


def save_plot_facet(
    adatas,
    out_file_base,
    color_var=None,
    colors_quantitative=True,
    colors_large_palette=COLORS_LARGE_PALLETE,
    drop_legend=-1
):
    """Save a plot."""
    n_params = len(adatas)
    fig, grid = panel_grid(
        hspace=0.125*n_params,
        wspace=None,
        ncols=4,
        num_panels=len(adatas)
    )
    i__ax = 0
    for adata in adatas:

        # Set the facet title.
        plt_title = adata.uns['adata_name']
        plt_title = plt_title.rstrip()

        # Get the proper axis for this plot.
        ax = plt.subplot(grid[i__ax])

        # We could avoid this line by setting basis=i__umap, but then not
        # consistent axis labels (e.g., X_umap__n_neighbors_151).
        # adata.obsm['X_umap'] = adata.obsm['X_umap']

        # If not colors_quantitative, then boot up categorical plots.
        legend_loc = 'right margin'
        color_palette = 'viridis'
        if colors_quantitative is False:
            # Cast to category - required for booleans.
            adata.obs[color_var] = adata.obs[color_var].astype('category')
            n_categories = len(adata.obs[color_var].cat.categories)
            color_palette = None
            if n_categories <= len(plt.get_cmap('Dark2').colors):
                color_palette = 'Dark2'
            elif n_categories <= len(colors_large_palette):
                color_palette = colors_large_palette
            if drop_legend >= 0 and n_categories >= drop_legend:
                legend_loc = None

        sc.pl.umap(
            adata=adata,
            color=color_var,
            palette=color_palette,
            alpha=0.4,
            title=plt_title,
            legend_loc=legend_loc,
            ax=ax,
            show=False
        )

        i__ax += 1

    fig.savefig(
        '{}-color__{}.png'.format(out_file_base, color_var),
        #dpi=300,
        bbox_inches='tight'
    )


def main():
    """Run CLI."""
    # (3) If gene list provided, will plot markers across clusters.
    parser = argparse.ArgumentParser(
        description="""
            Reads original H5 file and CellBender results. Generates a few QC
            plots:
            (1) Subtracts the CellBender counts matrix from the original
                counts matrix to see what was removed by CellBender.
            (2) Compares clusters before and after.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-txd', '--tenxdata_path',
        action='store',
        dest='txd',
        required=True,
        help='Path to data in 10x format.'
    )

    parser.add_argument(
        '-h5cb', '--h5_cellbender',
        action='store',
        dest='h5_cb',
        required=True,
        help='H5 AnnData file.'
    )

    # parser.add_argument(
    #     '-gcb', '--genome_cellbender',
    #     action='store',
    #     dest='g',
    #     default="background_removed",
    #     help='Filter expression to genes within this genome.\
    #         For CellBender v1, this will be "background_removed".\
    #         For CellBender v2, this will be "matrix".'
    # )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='output_file',
        default='cellbender_results',
        help='Basename of output files.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=1,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    #sc.set_figure_params(dpi_save=300)

    out_base = options.output_file

    # See what the results file looks likely
    # import tables
    # with tables.open_file('output_filtered.h5') as f:
    #     print(f)  # display the structure of the h5 file

    adata = load_anndata_from_input_and_output(
        input_dir_10x=options.txd,
        cellbender_h5=options.h5_cb,
        analyzed_barcodes_only=True
    )

    # Plot the difference in counts across cells for the top 100 most
    # different genes
    plot_difference(adata, plot_name=out_base)
    # Make a scatter of the inferred ambient RNA by number counts removed
    plot_ambient_by_difference(adata, plot_name=out_base)

    # Do basic clustering before and after cellbender
    # Do basic clustering before and after cellbender
    # NOTE: because we filter for low quality cells the dims no longer match
    adata_clustered_raw = basic_cluster(
        adata.copy(),
        'counts_raw'
    )
    adata_clustered_raw.uns['adata_name'] = 'raw'
    print("Finished adata_clustered_raw")
    adata_clustered_cb = basic_cluster(
        adata.copy(),
        'counts_cellbender'
    )
    adata_clustered_cb.uns['adata_name'] = 'cellbender'
    print("Finished adata_clustered_cb")

    # Plot UMAPS of both side by side
    save_plot_facet(
        adatas=[adata_clustered_raw, adata_clustered_cb],
        out_file_base='{}-umap'.format(out_base),
        # color_var=color_var,
        # colors_quantitative=True,
        drop_legend=False
    )
    save_plot_facet(
        adatas=[adata_clustered_raw, adata_clustered_cb],
        out_file_base='{}-umap'.format(out_base),
        color_var='cluster',
        # colors_quantitative=True,
        drop_legend=False
    )

    # Dot plot of gene expression across cell types


def dev():
    f_10x = "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data/ti/healthy/OTARscRNA8356110/raw_feature_bc_matrix"
    f_cb = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/leland_keras/cellbender/output_filtered.h5"
    # adata_10x = sc.read_10x_mtx(
    #     path=f_10x,
    #     # var_names='gene_symbols',
    #     var_names='gene_ids',
    #     make_unique=False
    # )
    # adata_cb = sc.read_10x_h5(filename=f_cb, genome='matrix')

    adata = load_anndata_from_input_and_output(
        input_dir_10x=f_10x,
        cellbender_h5=f_cb,
        analyzed_barcodes_only=True
    )
    # Plot the difference in counts across cells for the top 100 most
    # different genes
    plot_difference(adata)
    print("Finished plot difference")
    # Make a scatter of the inferred ambient RNA by number counts removed
    plot_ambient_by_difference(adata)

    # Do basic clustering before and after cellbender
    # NOTE: because we filter for low quality cells the dims no longer match
    adata_clustered_raw = basic_cluster(
        adata.copy(),
        'counts_raw'
    )
    adata_clustered_raw.uns['adata_name'] = 'raw'
    print("Finished adata_clustered_raw")
    adata_clustered_cb = basic_cluster(
        adata.copy(),
        'counts_cellbender'
    )
    adata_clustered_cb.uns['adata_name'] = 'cellbender'
    print("Finished adata_clustered_cb")

    # Do plot umap clustering before and after cellbender
    # adata_clustered_raw.obsm['X_umap__raw'] = adata_clustered_raw.obsm[
    #     'X_umap'
    # ].copy()
    # adata_clustered_raw.obsm['X_umap__cellbender'] = adata_clustered_cb.obsm[
    #     'X_umap'
    # ]
    save_plot_facet(
        adatas=[adata_clustered_raw, adata_clustered_cb],
        out_file_base='test',
        # color_var=color_var,
        # colors_quantitative=True,
        drop_legend=False
    )
    save_plot_facet(
        adatas=[adata_clustered_raw, adata_clustered_cb],
        out_file_base='test2',
        color_var='cluster',
        # colors_quantitative=True,
        drop_legend=False
    )

    # Dotplot of marker genes
    # _ = sc.pl.dotplot(
    #     adata=adata,
    #     var_names=marker_genes_found,
    #     groupby='cluster',
    #     gene_symbols='gene_symbols',
    #     dendrogram=run_dendrogram,
    #     show=False,
    #     use_raw=False,
    #     log=False,
    #     color_map='Blues',
    #     save='-{}-{}-{}.png'.format(
    #         out_file_base,
    #         i_out,
    #         data_scale
    #     )
    # )
    # _ = sc.pl.dotplot(
    #     adata=adata,
    #     var_names=marker_genes_found,
    #     groupby='cluster',
    #     gene_symbols='gene_symbols',
    #     dendrogram=run_dendrogram,
    #     show=False,
    #     standard_scale='var',  # Scale color between 0 and 1
    #     use_raw=False,
    #     color_map='Blues',
    #     save='-{}-{}-{}_standardized.png'.format(
    #         out_file_base,
    #         i_out,
    #         data_scale
    #     )
    # )


if __name__ == '__main__':
    #dev()
    main()
