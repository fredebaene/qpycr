import numpy as np
import pandas as pd

def calculate_avg_cqs(
    cqs: pd.DataFrame
):
    """
    This function takes a data frame with at least the following 
    columns: (1) a column containg the raw Cq values, (2) a column 
    containing the sample identifiers, and (3) a column containing 
    the target names.

    For every sample-target combination, the average Cq value and 
    the standard deviation of the technical replicates are calculated.

    Args:
        cqs (pd.DataFrame): data frame containg sample identifiers, 
            target names, and raw Cq values.

    Returns:
        pd.DataFrame: data frame containg the mean and standard 
            deviation of the technical replicates.
    """
    # Check if the arguments are valid
    if not isinstance(cqs, pd.DataFrame):

        raise TypeError("The argument passed to `cqs` must be a `pd.DataFrame`")

    # Calculate mean Cq for each sample-target combination
    avg_cqs = cqs.groupby(
            by=["sample", "target"]
        ).agg(
            func={"cq": ["mean", "std"]}
        ).droplevel(
            level=0,
            axis=1
        ).reset_index()

    avg_cqs = cqs[
            ["sample", "target"]
        ].drop_duplicates().reset_index(
            drop=True
        ).merge(
            right=avg_cqs,
            how="left",
            on=["sample", "target"]
        )

    return avg_cqs


def calculate_d_cqs(
    cqs: pd.DataFrame,
    internal_controls: list
):
    """
    This function takes a data frame with at least the following 
    columns: (1) a column containg the raw Cq values, (2) a column 
    containing the sample identifiers, and (3) a column containing 
    the target names.

    The mean and standard deviation of the technical replicates are 
    calculated for every sample-target combination using the 
    <b>calculate_avg_cqs()</b> function. For each sample, the mean Cq 
    values of the targets in the <b>internal_controls</b> list are 
    averaged. This average is subtracted of the mean Ct of every 
    target of interest. The fold change is calculated as 2 to the 
    negative power of the delta Cq value.

    Args:
        cqs (pd.DataFrame): data frame containg sample identifiers, 
            target names, and raw Cq values.
        internal_controls (list): list specifying internal controls.

    Returns:
        pd.DataFrame: data frame containg the mean and standard 
            deviation of the technical replicates, the mean of 
            the internal controls, the delta Cq, and the fold change 
            of the target of interest relative to the internal 
            control(s).
    """
    # Check if the arguments are valid
    if not isinstance(cqs, pd.DataFrame):

        raise TypeError("The argument passed to `cqs` must be a `pd.DataFrame`")

    if not isinstance(internal_controls, list):

        raise TypeError("The argument passed to `internal_controls` must be a `list`")
    
    # Calculate mean Cq for each sample-target combination
    avg_cqs = calculate_avg_cqs(cqs)

    # Calculate mean of the mean Cq values for the internal controls
    avg_ic_cqs = avg_cqs[
            avg_cqs["target"].isin(internal_controls)
        ].groupby(
            by=["sample"]
        ).agg(
            func={"mean": "mean"}
        ).rename(
            mapper={"mean": "mean_ic"},
            axis=1
        )

    # Append mean Cq values of internal controls
    avg_cqs = avg_cqs.merge(
            right=avg_ic_cqs,
            how="left",
            on="sample"
        )
        
    # Calculate delta Cq and fold change
    avg_cqs["d_cq"] = np.where(
            ~avg_cqs["target"].isin(internal_controls),
            avg_cqs["mean"] - avg_cqs["mean_ic"], np.NaN
        )

    return avg_cqs

def calculate_dd_cqs(
    cqs: pd.DataFrame,
    internal_controls: list,
    calibrator: str
):

    # Check if the arguments are valid
    if not isinstance(cqs, pd.DataFrame):

        raise TypeError("The argument passed to `cqs` must be a `pd.DataFrame`")

    if not isinstance(internal_controls, list):

        raise TypeError("The argument passed to `internal_controls` must be a `list`")

    if not isinstance(calibrator, str):

        raise TypeError("The argument passed to `calibrator` must be a `str`")

    # Calculate the delta Cq values
    dd_cqs = q.calculate_d_cqs(cqs, internal_controls)

    # Extract the delta Cq values for each TOI in the calibrator
    d_cqs_ic = dd_cqs[
            (dd_cqs["sample"] == calibrator)
            & ~(dd_cqs["target"].isin(internal_controls))
        ][
            ["target", "d_cq"]
        ].rename(
            mapper={"d_cq": "d_cq_ic"},
            axis=1
        )

    dd_cqs = dd_cqs.merge(
            right=d_cqs_ic,
            how="left",
            on="target"
        )

    dd_cqs["dd_cq"] = dd_cqs["d_cq"] - dd_cqs["d_cq_ic"]
    dd_cqs["fc"] = 2**-dd_cqs["dd_cq"]

    return dd_cqs