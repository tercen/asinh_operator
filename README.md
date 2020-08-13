# asinh operator

#### Description
`asinh` operator performs an inverse hyperbolic sine on values

##### Usage
Input projection|.
---|---
`y-axis` | values as input to the asinh operator
`row`    | channels, and scale (optional)
`col`    | event, for example

Input parameters|.
---|---
`scale`             | numeric, the scaling factor to use after the asinh transformation, default is 5
`scale per channel` | boolean, indicates if a different scaling factor per channel is performed, default is false

Output relations|.
---|---
`asinh`| numeric, output transformation per `row` and `col`.

##### Details
An asinh of values is performed followed by a scaling. One data point per cell is required as input

```r
asinh(value/scale)
```
A scale of 5 is recommended for cytof measurement and 150 for flowcyto measurements.

If used in the `scale per channel` mode then an extra factor indicating the scaling value is required in the cross-tab projection on the `row` dimension in the second position.

#### References
see the `base::asinh` function of the R package for the documentation,
https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions


##### See Also

#### Examples
