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
Automatically estimates optimal cofactors per channel using percentile-based estimation. This method analyzes the data distribution for each channel and determines the best transformation parameter based on the 5th percentile of positive values.

Requirements for auto method:
- Channel names must be provided in the row projection
- Estimated cofactors are constrained between 1 and 1000

##### References

- See the `base::asinh` [function of the R package for more information](https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions).
