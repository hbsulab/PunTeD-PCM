# Punctuated Tensor Decomposed Pair Clustering Method

This repository provides an implementation of the Punctuated Tensor Decomposed Pair Clustering Method (PunTed-PCM).

>requirements:\
>numpy>=1.26.0\
>pandas>=2.2.0\
>matplotlib>=3.10.0\
>seaborn>=0.13.2\
>scipy>=1.15.0\
>tensorly>=0.8.0\
>xlsxwriter>=3.2.0\
>scipy>=1.15.0

**Data Preparation**

Prepare your mutation data in the format: D614G;N501Y (multiple mutations separated by semicolons). Check https://github.com/hbsulab/L-index for more information of data preparation. Place your data under the path directory.

**Generate Composite Matrix**

>python Composite_matrix_TSVGISAID_HPC3.py

**Compute Covariance Networks**

Monthly-dependent covariance:

>python cov_compu.py

Covariance time series:

>python cov_time_series.py

**Calculate Hilbert Frequencies (MATLAB)**

Run the following MATLAB scripts:

>hilbert_compu.m


>save_fre.m

**Construct Coevolutionary Landscape**

>python seqeunce_space.py

**Skewness Calculation**

>python skew_sap_cal.py
