#!/usr/bin/env python

__author__ = 'Leland Taylor'
__date__ = '2020-07-03'
__version__ = '0.0.1'

import argparse
#import csv
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import random
import numpy as np
import pandas as pd

import plotnine as plt9
from kneed import KneeLocator
from distutils.version import LooseVersion
from scipy.interpolate import UnivariateSpline


import scanpy as sc


# To resolve strange TclError for interactive job
import matplotlib
matplotlib.use('Agg')  # Agg for png and pdf for pdf

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)

# Set valid methods for estimators
valid_methods = [
    'dropletutils::barcoderanks::inflection',
    'dropletutils::barcoderanks::knee',
    'kneedle::spline=None',
    'kneedle::spline=interp1d',
    'kneedle::spline=polynomial',
    'expected',
    'expected_umicutoff'
]


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def estimate_cutoffs_plot(
    output_file,
    df_plt,
    df_estimate_ncells,
    df_fit=None,
    scale_x_log10=False,
    save_plot=True,
    add_text=False
):
    """Plot UMI counts by sorted cell barcodes."""
    if min(df_plt['umi_counts']) <= 0:
        fix_log_scale = min(df_plt['umi_counts']) + 1
        df_plt['umi_counts'] = df_plt['umi_counts'] + fix_log_scale
    if add_text:
        df_estimate_ncells['add_text_y'] = np.random.randint(
            low=df_plt['umi_counts'].min()-25,
            high=df_plt['umi_counts'].max()-25,
            size=df_estimate_ncells.shape[0]
        )
    gplt = plt9.ggplot()
    gplt = gplt + plt9.theme_bw()
    if len(df_plt) <= 50000:
        gplt = gplt + plt9.geom_point(
            mapping=plt9.aes(x='barcode', y='umi_counts'),
            data=df_plt,
            alpha=0.05,
            size=0.1
        )
    else:
        gplt = gplt + plt9.geom_line(
            mapping=plt9.aes(x='barcode', y='umi_counts'),
            data=df_plt,
            alpha=0.25,
            size=0.75,
            color='grey'
        )
    gplt = gplt + plt9.geom_vline(
        mapping=plt9.aes(xintercept='n_cells', color='method'),
        data=df_estimate_ncells,
        alpha=0.75,
        linetype='dashdot'
    )
    if add_text:
        gplt = gplt + plt9.geom_text(
            mapping=plt9.aes(
                x='n_cells',
                y='add_text_y',
                label='n_cells',
                color='method'
            ),
            data=df_estimate_ncells,
            alpha=0.75
        )
    gplt = gplt + plt9.scale_color_brewer(
        palette='Dark2',
        type='qual'
    )
    if scale_x_log10:
        gplt = gplt + plt9.scale_x_continuous(
            trans='log10',
            labels=comma_labels,
            minor_breaks=0
        )
    else:
        gplt = gplt + plt9.scale_x_continuous(
            labels=comma_labels,
            minor_breaks=0
        )
    gplt = gplt + plt9.scale_y_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.labs(
        title='',
        y='UMI counts',
        x='Barcode index, sorted by UMI count',
        color='Cutoff'
    )
    # Add the fit of the droplet utils model
    if df_fit:
        gplt = gplt + plt9.geom_line(
            mapping=plt9.aes(x='x', y='y'),
            data=df_fit,
            alpha=1,
            color='yellow'
        )
    if save_plot:
        gplt.save(
            '{}.png'.format(output_file),
            dpi=300,
            width=5,
            height=4
        )
    return gplt


def dropletutils_cutoff(df_analysis):
    """Calculate the knee and inflection point as in DropletUtils."""
    # Calculate the knee and inflection point as in DropletUtils::barcodeRanks
    # https://github.com/MarioniLab/DropletUtils/blob/master/R/barcodeRanks.R
    # Numerical differentiation to identify bounds for spline fitting.
    # The upper/lower bounds are defined at the plateau and inflection.
    #
    # The lower bound on the total UMI count, at or below which all barcodes
    # are assumed to correspond to empty droplets
    # df_analysis = df.loc[df['umi_counts'] >= lower_bound, :]
    if min(df_analysis['umi_counts']) <= 0:
        fix_log_scale = min(df_analysis['umi_counts']) + 1
        df_analysis['umi_counts'] = df_analysis['umi_counts'] + fix_log_scale
    d1n = np.diff(
        np.log10(df_analysis['umi_counts'])
    )/np.diff(np.log10(df_analysis['barcode']))
    # NOTE: these are the indexes of the edge of the plot
    right_edge = np.argmin(d1n)
    left_edge = np.argmax(d1n[0:right_edge])

    # We restrict to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation
    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee
    # point
    x = np.log10(
        df_analysis.iloc[left_edge:right_edge, :]['barcode'].values
    )
    y = np.log10(
        df_analysis.iloc[left_edge:right_edge, :]['umi_counts'].values
    )
    dropletutils_fit = UnivariateSpline(x, y)
    fit1 = dropletutils_fit.derivative(n=1)
    fit2 = dropletutils_fit.derivative(n=2)
    curvature = fit1(x) / np.power(1 + np.power(fit2(x), 2), 1.5)
    dropletutils_knee = 10 ** y[np.argmin(curvature)]
    dropletutils_inflection = df_analysis.iloc[right_edge, :]['umi_counts']
    dropletutils_knee_ncells = 10 ** x[np.argmin(curvature)]
    dropletutils_inflection_ncells = df_analysis.iloc[right_edge, :]['barcode']

    result = [
        {
            'method': 'dropletutils::barcoderanks::knee',
            'umi_counts_cutoff': dropletutils_knee,
            'n_cells': int(dropletutils_knee_ncells)
        },
        {
            'method': 'dropletutils::barcoderanks::inflection',
            'umi_counts_cutoff': dropletutils_inflection,
            'n_cells': int(dropletutils_inflection_ncells)
        }
    ]

    return right_edge, left_edge, result


def kneedle_cutoff(df_analysis, verbose=True):
    """Get kneedle estimate for what is a cell."""
    results = []
    # Get kneedle estimate using raw data OR fit smoothing spline. Also
    # maximize sensitivity.
    # NOTE: Much better estimates with kneedle (that is similar to cellranger3)
    # if this is done the original value space - not log space.
    kneedle_dict = {}
    for fit in [None, 'interp1d', 'polynomial']:
        key = 'kneedle::spline={}'.format(fit)
        sensitivity = 5000
        while True:
            if fit:
                kneedle_dict[key] = KneeLocator(
                    df_analysis['barcode'].values,
                    df_analysis['umi_counts'].values,
                    curve='convex',
                    direction='decreasing',
                    S=sensitivity,
                    interp_method=fit
                )
            else:
                kneedle_dict[key] = KneeLocator(
                    df_analysis['barcode'].values,
                    df_analysis['umi_counts'].values,
                    curve='convex',
                    direction='decreasing',
                    S=sensitivity
                )
            if kneedle_dict[key].knee is None:
                sensitivity -= 100
            else:
                if verbose:
                    print('S:\t{}\nknee:\t{}\nelbow:\t{}'.format(
                        sensitivity,
                        round(kneedle_dict[key].knee, 3),
                        round(kneedle_dict[key].elbow, 3)
                    ))
                results.append({
                    'method': key,
                    'umi_counts_cutoff': df_analysis.loc[
                        df_analysis['barcode'].values == int(
                            kneedle_dict[key].knee
                        ),
                        'umi_counts'
                    ].values[0],
                    'n_cells': int(kneedle_dict[key].knee),
                    'sensitivity': sensitivity
                    # 'elbow': 10 ** kneedle_dict[key].elbow  # same as knee
                })
                break
    return results


def get_final_estimates(
    df_estimates,
    method='dropletutils::barcoderanks::inflection',
    subtract_n_factor=0,
    verbose=False
):
    # Check to see if list. If list then it must be three values:
    # (1) estimator 1
    # (2) estimator 2
    # (3) factor to split the difference between the estimates
    split_factor = 1.0

    if not isinstance(method, str):     # assume it is a list
        if len(method) == 1:
            pass
        elif len(method) == 3:
            split_factor = float(method.pop())
        else:
            raise Exception(
                'Error in method ({}). \
                Should be:\testimate1,estimate2,split_factor'.format(
                    ','.join(method)
                )
            )
    else:
        method = [method]

    # Check the methods
    apply_subtract_factor = True
    for met in method:
        if met not in valid_methods:
            raise Exception('Invalid method:\t{}.'.format(met))
        # Do not apply the subtract factor if we have an "expected" method
        if met in ['expected', 'expected_umicutoff']:
            apply_subtract_factor = False

    # Now calculate the final estimate
    estimates = []
    for met in method:
        filt = df_estimates['method'] == met
        estimates.append(df_estimates.loc[filt, 'n_cells'].values[0])
    if len(estimates) == 1:
        estimated_n = estimates[0]
    elif len(estimates) == 2:
        estimated_n = abs(estimates[0] - estimates[1]) * split_factor
        estimated_n += min(estimates)
    else:
        raise Exception('Error in computing estimate.')

    # If method not an exact one, then apply the subtract_factor
    if apply_subtract_factor:
        estimated_n -= subtract_n_factor

    if verbose:
        pass

    return int(estimated_n), ','.join(method)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Estimates cutoffs for cellbender remove background from total UMI
            counts. NOTE: different cutoff estimators work better on different
            tissues/protocols that have slightly different UMI count plots.
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
        help='Path to 10x data.'
    )
    # Let the user set the expected number of cells / droplets
    parser.add_argument(
        '--expected_ncells',
        action='store',
        dest='expected_ncells',
        default=0,
        type=int,
        help='Expected number of cells. If != 0, this number of cells \
            will be used for expected_cells.txt. (default: %(default)s).'
    )
    # Let the user set the expected number of cells / droplets
    parser.add_argument(
        '--expected_nemptydroplets',
        action='store',
        dest='expected_ndroplets',
        default=0,
        type=int,
        help='Expected total number of empty droplets. If != 0, this number of \
            will be used for total_droplets_included.txt. \
            (default: %(default)s).'
    )
    # Let the user set expected number of droplets based on UMI counts
    parser.add_argument(
        '--expected_nemptydroplets_umi_cutoff',
        action='store',
        dest='expected_ndroplets_umicutoff',
        default=0,
        type=int,
        help='Set expected number of empty droplets to droplets with UMI \
            cutoffs below this value. \
            If != 0 and expected_total_ndroplets == 0, the corresponding UMI \
            barcode ranks will be used for total_droplets_included.txt. \
            (default: %(default)s).'
    )
    # Estimate the expected number of cells
    parser.add_argument(
        '--method_estimate_ncells',
        action='store',
        dest='method_ncells',
        default='dropletutils::barcoderanks::inflection',
        help='Method to use to estimate number of cells if \
            expected_ncells == 0. Valid options = [{}]. \
            (default: %(default)s)'.format(','.join(valid_methods))
    )
    parser.add_argument(
        '--method_estimate_nemptydroplets',
        action='store',
        dest='method_ndroplets',
        default='dropletutils::barcoderanks::knee',
        help='Method to use to estimate total number of droplets if \
            total_ndroplets == 0. Valid options = [{}]. \
            One can also calculate the difference between two metrics \
            (default: %(default)s)'.format(','.join(valid_methods))
    )
    parser.add_argument(
        '--lower_bound_umis_estimate_ncells',
        action='store',
        dest='lb_estimated_ncells',
        default=1000,
        type=int,
        help='The lower bound on the total UMI count at or below which all \
            barcodes are assumed to correspond to empty droplets (for \
            estimating expected ncells). \
            (default: %(default)s)'
    )
    parser.add_argument(
        '--upper_bound_umis_estimate_nemptydroplets',
        action='store',
        dest='ub_estimated_ndroplets',
        default=250,
        type=int,
        help='The upper bound on the total UMI count used to estimate cells \
            with ambient RNA that will be included in cellbender to \
            estimate the background. \
            (default: %(default)s)'
    )
    parser.add_argument(
        '--lower_bound_umis_estimate_nemptydroplets',
        action='store',
        dest='lb_estimated_ndroplets',
        default=10,
        type=int,
        help='The lower bound on the total UMI count used to estimate cells \
            with ambient RNA that will be included in cellbender to \
            estimate the background. \
            (default: %(default)s)'
    )
    parser.add_argument(
        '--estimate_nemptydroplets_add_umifactor',
        action='store',
        dest='add_umifactor_estimated_ndroplets',
        default=0,
        type=int,
        help='Add this number from the estimated UMI cutoff \
            to get the final total_ndroplets. Adding shifts the cutoff on the \
            UMI curve to the left. This option only does something if \
            expected_nemptydroplets and expected_nemptydroplets_umi_cutoff \
            are 0. If the final value is outside the range of \
            ub_estimated_ndroplets, then it is set to ub_estimated_ndroplets. \
            (default: %(default)s)'
    )
    parser.add_argument(
        '--estimate_nemptydroplets_subtract_dropletfactor',
        action='store',
        dest='subtract_dropletfactor_estimated_ndroplets',
        default=0,
        type=int,
        help='Subtract this number from the estimated count cutoff \
            to get the final total_ndroplets. Subtract shifts the cutoff on \
            the UMI curve to the left. This option only does something if \
            expected_nemptydroplets and expected_nemptydroplets_umi_cutoff \
            are 0. (default: %(default)s)'
    )
    parser.add_argument(
        '--estimate_nemptydroplets_min_nemptydroplets',
        action='store',
        dest='min_ndroplets',
        default=0,
        type=int,
        help='Minimum value for empty drops cutoff. \
            If estimate <= this value, then we set the threshold to this \
            value. (default: %(default)s)'
    )
    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata',
        help='Basename of output files, assuming output in current \
            working directory. \
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Get clean names
    method_ncells = options.method_ncells.split(',')
    method_ndroplets = options.method_ndroplets.split(',')
    expected_ndroplets_umicutoff = options.expected_ndroplets_umicutoff
    lb_estimated_ncells = options.lb_estimated_ncells
    lb_estimated_ndroplets = options.lb_estimated_ndroplets
    ub_estimated_ndroplets = options.ub_estimated_ndroplets
    add_umifactor_estimated_ndroplets = (
        options.add_umifactor_estimated_ndroplets
    )
    subtract_dropletfactor_estimated_ndroplets = (
        options.subtract_dropletfactor_estimated_ndroplets
    )
    min_ndroplets = options.min_ndroplets
    output_file = options.of
    verbose = True

    # Check for something wonky
    if (
        add_umifactor_estimated_ndroplets != 0
    ) and (
        subtract_dropletfactor_estimated_ndroplets != 0
    ):
        raise Exception('Error in computing estimate.')

    # Load single cell file
    file = options.txd
    adata = sc.read_10x_mtx(
        path=file,
        # var_names='gene_symbols',
        var_names='gene_ids',
        make_unique=False
    )

    # Make the umi dataframe
    df = pd.DataFrame({
        'umi_counts': np.sort(
            (np.array(adata.X.sum(axis=1))).flatten()
        )[::-1]
    })
    df['barcode'] = df.index + 1

    # Make an analysis df using lower bound
    df_analysis = df.loc[df['umi_counts'] >= lb_estimated_ncells, :]
    df_analysis = df_analysis.reset_index(drop=True)
    # Get dropletutils estimates
    right_edge, left_edge, result_dropletutils = dropletutils_cutoff(
        df_analysis
    )
    # Get kneedle estimates
    result_kneedle = kneedle_cutoff(df_analysis, verbose=False)
    # Save the results to a list we will turn into a pandas matrix
    cell_estimate_outdict = []
    for i in result_dropletutils:
        cell_estimate_outdict.append(i)
    for i in result_kneedle:
        cell_estimate_outdict.append(i)

    # Set expected number of cells if we have some
    if options.expected_ncells != 0:
        # NOTE: with specific esimates, we should use the raw df
        method_ncells = 'expected'
        cell_estimate_outdict.append({
            'method': 'expected',
            'umi_counts_cutoff': df['umi_counts'][options.expected_ncells],
            'n_cells': options.expected_ncells
        })
        # If we have expected number of cells, then go ahead and
        # also estimate like CellRanger2.
        # Assumes a ~10-fold range of library sizes for real cells and
        # estimates this range using the expected number of cells.
        # https://scrnaseq-course.cog.sanger.ac.uk/website/processing-raw-scrna-seq-data.html
        # 99th percentile of top n_cells divided by 10
        #
        # NOTE: with specific esimates, we should use the raw df
        cellranger_expected_knee = df['umi_counts'][
            round(0.01*options.expected_ncells)
        ]/10
        cell_estimate_outdict.append({
            'method': 'cellrangerv2::expected',
            'umi_counts_cutoff': cellranger_expected_knee,
            'n_cells': (df['umi_counts'] >= cellranger_expected_knee).sum()
        })

    # Make a dataframe of all of the different knee calculations
    # NOTE: here knee is the number of cells that would be kept at that
    # threshold
    df_estimate_ncells = pd.DataFrame(cell_estimate_outdict)
    # Add the fit of the droplet utils model
    # df_fit = pd.DataFrame({
    #     'y': 10 ** dropletutils_fit(x),
    #     'x': 10 ** x
    # })

    # Get the final number of estimated cells
    estimated_ncells, estimated_ncells_method = get_final_estimates(
        df_estimates=df_estimate_ncells,
        method=method_ncells,
        subtract_n_factor=0,
        verbose=verbose
    )
    # Add the final estimate to the output dict
    cell_estimate_outdict.append({
        'method': 'estimated_ncells',
        'umi_counts_cutoff': 0,
        'n_cells': estimated_ncells
    })
    df_estimate_ncells = pd.DataFrame(cell_estimate_outdict)

    # Save all of our estimates
    df_estimate_ncells.to_csv(
        '{}-cell_estimate_cutoff.tsv.gz'.format(output_file),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )
    # Save the final estimates
    with open('{}-expected_cells.txt'.format(output_file), 'w') as f:
        f.write(str(int(estimated_ncells)))

    # Make a plot of the different cutoffs
    _ = estimate_cutoffs_plot(
        '{}-cell_estimate_cutoffs-zoomed'.format(output_file),
        df_analysis,
        df_estimate_ncells,
        # df_fit=df_fit,
        scale_x_log10=False
    )
    _ = estimate_cutoffs_plot(
        '{}-cell_estimate_cutoffs'.format(output_file),
        df,
        df_estimate_ncells,
        # df_fit=df_fit,
        scale_x_log10=True
    )

    # Run a similar proceedure but for total-droplets-included.
    # ...these will be the lower bound of the droplets used to estimate the
    # background signals
    # Make an analysis df using lower bound
    df_analysis = df.loc[df['umi_counts'] >= lb_estimated_ndroplets, :]
    df_analysis = df_analysis.loc[
        df_analysis['umi_counts'] <= ub_estimated_ndroplets,
        :
    ]
    df_analysis = df_analysis.reset_index(drop=True)
    # Get dropletutils estimates
    right_edge, left_edge, result_dropletutils = dropletutils_cutoff(
        df_analysis
    )
    # Get kneedle estimates
    result_kneedle = kneedle_cutoff(df_analysis, verbose=False)
    # Save the results to a list we will turn into a pandas matrix
    total_droplets_estimate_outdict = []
    for i in result_dropletutils:
        total_droplets_estimate_outdict.append(i)
    for i in result_kneedle:
        total_droplets_estimate_outdict.append(i)

    # If we have an expected total number of droplets use that.
    if options.expected_ndroplets != 0:
        # NOTE: with specific esimates, we should use the raw df
        method_ndroplets = 'expected'
        total_droplets_estimate_outdict.append({
            'method': 'expected',
            'umi_counts_cutoff': df['umi_counts'][options.expected_ndroplets],
            'n_cells': options.expected_ndroplets
        })
    # If we have expected total number of droplets based on total UMIs use that
    if expected_ndroplets_umicutoff != 0:
        # NOTE: with specific esimates, we should use the raw df
        total_droplets_estimate_outdict.append({
            'method': 'expected_umicutoff',
            'umi_counts_cutoff': expected_ndroplets_umicutoff,
            'n_cells': df.loc[
                df['umi_counts'] <= expected_ndroplets_umicutoff,
                'barcode'
            ].min()
        })
        # Never over-ride expected_ndroplets
        if options.expected_ndroplets == 0:
            method_ndroplets = 'expected_umicutoff'

    # Make a dataframe of all of the estimates
    df_estimate_ndroplets = pd.DataFrame(total_droplets_estimate_outdict)

    # Now, get the final number of estimated droplets
    estimated_ndroplets, estimated_ndroplets_method = get_final_estimates(
        df_estimates=df_estimate_ndroplets,
        method=method_ndroplets,
        subtract_n_factor=subtract_dropletfactor_estimated_ndroplets,
        verbose=verbose
    )
    # We have one last edit and that is if the user has told us to subtract
    # something based on umis
    if estimated_ndroplets_method not in ['expected', 'expected_umicutoff']:
        if add_umifactor_estimated_ndroplets:
            estimated_ndroplets_umicutoff = df.loc[
                df['barcode'] == estimated_ndroplets,
                'umi_counts'
            ].values[0]
            # if (
            #     estimated_ndroplets_umicutoff >= lb_estimated_ndroplets - 10
            # ) and (
            #     estimated_ndroplets_umicutoff <= lb_estimated_ndroplets + 10
            # ):
            estimated_ndroplets_umicutoff += add_umifactor_estimated_ndroplets
            estimated_ndroplets_umicutoff = max(
                estimated_ndroplets_umicutoff,
                ub_estimated_ndroplets
            )
            estimated_ndroplets = df.loc[
                df['umi_counts'] == estimated_ndroplets_umicutoff,
                'barcode'
            ].values[0]
    # If the estimated ndroplets is less than the user specified min, then set
    # to the user specified min.
    if (estimated_ndroplets != 0) and (estimated_ndroplets < min_ndroplets):
        estimated_ndroplets = min_ndroplets

    # Add the final estimate to the output
    total_droplets_estimate_outdict.append({
        'method': 'estimated_ndroplets',
        'umi_counts_cutoff': 0,
        'n_cells': estimated_ndroplets
    })
    df_estimate_ndroplets = pd.DataFrame(total_droplets_estimate_outdict)
    # Save all of our estimates
    df_estimate_ndroplets.to_csv(
        '{}-total_droplets_cutoff.tsv.gz'.format(output_file),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )
    # Save total_droplets_included
    with open('{}-total_droplets_included.txt'.format(output_file), 'w') as f:
        f.write(str(int(estimated_ndroplets)))

    # Make a plot of the different cutoffs
    _ = estimate_cutoffs_plot(
        '{}-total_drops_estimate_cutoffs-zoomed'.format(output_file),
        df_analysis,
        pd.DataFrame(total_droplets_estimate_outdict),
        # df_fit=df_fit,
        scale_x_log10=False,
        save_plot=True
    )
    _ = estimate_cutoffs_plot(
        '{}-total_drops_estimate_cutoffs'.format(output_file),
        df,
        pd.DataFrame(total_droplets_estimate_outdict),
        # df_fit=df_fit,
        scale_x_log10=True
    )

    # Plot the final estimates togher
    final_cutoffs = []
    final_cutoffs.append({
        'method': 'Estimated expected # cells',
        'umi_counts_cutoff': 0,
        'n_cells': estimated_ncells
    })
    final_cutoffs.append({
        'method': 'Estimated total # droplets',
        'umi_counts_cutoff': 0,
        'n_cells': estimated_ndroplets
    })
    # print(final_cutoffs)
    # Optinally zoom in the final estimates
    # df_plt = df.loc[df['barcode'] <= total_droplets_cutoff_ncells + 100, :]
    df_plt = df.loc[df['umi_counts'] >= 1, :]
    _ = estimate_cutoffs_plot(
        '{}-final_estimates'.format(output_file),
        df_plt,
        pd.DataFrame(final_cutoffs),
        # df_fit=df_fit,
        scale_x_log10=False,
        save_plot=True,
        add_text=True
    )
    _ = estimate_cutoffs_plot(
        '{}-final_estimates-scale_x_log10'.format(output_file),
        df_plt,
        pd.DataFrame(final_cutoffs),
        # df_fit=df_fit,
        scale_x_log10=True,
        save_plot=True,
        add_text=True
    )


if __name__ == '__main__':
    main()
