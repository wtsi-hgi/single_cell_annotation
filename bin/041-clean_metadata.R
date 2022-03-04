#!/usr/bin/env Rscript

library(optparse)


main <- function() {
    optionList <- list(
        optparse::make_option(c("-i", "--input_file"),
            type = "character",
            help = paste0(
                "Input metadata file"
            )
        ),

        optparse::make_option(c("--out_file"),
            type = "character",
            default = "",
            help = paste0(
                "Name (and possibly path) of output file. Will have tsv.gz",
                " appended to it. If '' then add '-cleaned.tsv.gz' to",
                " pca_file.",
                " [default: %default]"
            )
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Cleans metadata."
        )
    )

    # a hack to fix a bug in optparse that won"t let you use positional args
    # if you also have non-boolean optional args:
    getOptionStrings <- function(parserObj) {
        optionStrings <- character()
        for (item in parserObj@options) {
            optionStrings <- append(optionStrings,
                                    c(item@short_flag, item@long_flag))
        }
        optionStrings
    }
    optStrings <- getOptionStrings(parser)
    arguments <- optparse::parse_args(parser, positional_arguments = TRUE)

    # read in the parameters
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }

    input_file <- arguments$options$input_file
    output_file_basename <- arguments$options$out_file

    # Read in the sample sheet
    dat_meta <- read.csv(
        input_file,
        sep = "\t"
    )

    # Write the metadata file for OTAR
    # df_final$late_early_ratio = df_final$f2.f1_ratio
    cols_keep <- c(
        'date_of_sample',
        'patient_id',
        'sanger_sample_id',
        'biopsy_type',
        'disease_status',
        'inflammation_status',
        'sex',
        'age',
        'smoking_status',
        'experimentalist',
        'protocol',
        'enzyme_lot_BLP',
        'bead_version',
        'bead_lot',
        'chip_version',
        'chip_lot',
        'gem_lot',
        'id_run',
        'lane',
        'library_id',
        'season_sample_collected',
        'month_sample_collected',
        #'experiment_order',
        'time_to_chromium_processing',
        'hours_to_chromium_processing',
        'smoker_at_time_of_biopsy',
        'inflammation_status_grouped',
        'late_early_ratio',
        'late_early_ratio_float',
        'early_late_ratio',
        'early_late_ratio_float',
        'Genotyping_ID',
        'medications_details'
    )
    mising_cols <- cols_keep[which(!cols_keep %in% colnames(dat_meta))]
    if (length(mising_cols) > 0) {
        print(mising_cols)
        stop()
    }
    base <- output_file_basename
    if (output_file_basename == "") {
        base <- paste0(
            gsub(
                ".gz",
                "",
                gsub(".tsv", "", basename(input_file))
            ),
            "-cleaned"
        )
    }
    gzfh <- gzfile(paste0(base, ".tsv.gz"), "w", compression = 9)
    write.table(
        dat_meta[cols_keep],
        gzfh,
        row.names=F,
        col.names=T,
        quote=F,
        sep="\t",
        na=""
    )
    close(gzfh)

}


# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
