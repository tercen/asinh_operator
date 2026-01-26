# asinh operator

##### Description

`asinh` operator performs an inverse hyperbolic sine on values. Also denoted as arcsinh, arsinh, or argsinh.

##### Usage

Input projection|.
---|---
`y-axis` | numeric, values required to be transformed by the asinh operator
`row`    | channels, and scale values (optional, required for manual/auto methods)
`col`    | event, for example

Input parameters|.
---|---
`method` | string, transformation method: 'fixed' (default), 'manual', or 'auto'
`scale`  | numeric, the scaling factor to use before the asinh transformation (only used when method is 'fixed'), default is 5

Output relations|.
---|---
`asinh`| numeric, output transformation per `row` and `col`.

##### Details

Values are scaled first and then asinh transformation is performed. One data point per cell is required as input.
```r
asinh(value/scale)
```

##### Methods

**fixed** (default)
Uses a global scale parameter for all channels. A scale of 5 is recommended for CyTOF measurements and 150 for flow cytometry measurements.

**manual**
Reads per-channel scale values from the row factor. Requires a second factor in the `row` dimension that provides the scaling value for each channel.

**auto**
Automatically estimates optimal cofactors per channel using the flowVS variance stabilization algorithm. This method finds the cofactor that minimizes variance heterogeneity across cell populations identified in each channel.

The algorithm:
1. Identifies cell populations as density peaks using kernel density estimation
2. For each candidate cofactor, transforms data with `asinh(x/cofactor)`
3. Calculates Bartlett's statistic to measure variance homogeneity across populations
4. Uses golden section search over logarithmic cofactor space (10⁻¹ to 10¹⁰) to find the optimal cofactor that minimizes Bartlett's statistic

Requirements for auto method:
- Channel names must be provided in the row projection
- At least 100 data points per channel for reliable peak detection
- Falls back to percentile-based estimation if peak detection fails

##### References

- See the `base::asinh` [function of the R package for more information](https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions).
- Azad A, Rajwa B, Pothen A (2016). "flowVS: Channel-Specific Variance Stabilization in Flow Cytometry." BMC Bioinformatics, 17:291. doi:10.1186/s12859-016-1083-9
