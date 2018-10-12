# asinh operator

#### Description
`asinh` operator performs an inverse hyperbolic sine of values.

##### Usage
Input projection|.
---|---
`y-axis` | values as input to the asinh operator

Input parameters|.
---|---
`scale` | numeric, the scaling factor to use after the asinh transformation, default is 5

Output relations|.
---|---
`asinh`| numeric, output transformation per data point

##### Details
An asinh of values is performed followed by a scaling.
```r
asinh(value)/scale
```
A scale of 5 is recommended for cytof measurement and 150 for flowcyto measurements.

#### References
see the `base::asinh` function of the R package for the documentation,
https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions


##### See Also

#### Examples
