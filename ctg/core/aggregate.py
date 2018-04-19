import logging
import pandas as pd
import numpy as np
from functools import reduce


def aggregate(counts_files, samples=None, sep="\t", header=0, comment="#"):
    """ Merge a list of dataframes together

        Args:
            counts_files (list) - a list of file paths to the counts files to merge
            samples (list) - a list of sample names to use for the merged columns
            sep (str) - the column seperator in the counts files
            header (int) - the position of the header in the counts files, None if no header
            comment (str) - comment character to ignore in the counts file
        Returns:
            merged_df (pd.DataFrame) - the merged counts data frame
        Raises:
            RunTimeError if sample names are not congruent
    """
    sample_pos = -1
    
    if samples is not None:
        if len(samples)!=len(counts_files):
            logging.error("Number of sample names is not the same length as the number of counts files.")
            raise RunTimeError("")

    # read in all counts files
    counts_df = [pd.read_csv(file, sep=sep, header=header, comment=comment) for file in counts_files]

    # overwrite the sample names if provided
    if samples:
        for i in range(counts_df):
            counts_df[i].columns[sample_pos] = samples[i]
    else:
        # check sample names are all different
        sample_names_from_files = [df.columns[sample_pos] for df in counts_df]

        if (len(np.unique(sample_names_from_files))<len(counts_files)):
            logging.error("Sample names in counts files are not unique. Fix or provide a list of sample names to use.")
            raise RunTimeError()


    # merge the dataframes together
    merged_df = reduce(lambda x, y: pd.merge(x,y), counts_dfs)

    return merged_df