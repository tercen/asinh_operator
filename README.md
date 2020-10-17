# asinh operator

##### Description

`asinh` operator performs an inverse hyperbolic sine on values

##### Usage

Input projection|.
---|---
`y-axis` | values required to be transformed by the asinh operator
`row`    | channels, and scale values (optional)
`col`    | event, for example

Input parameters|.
---|---
`scale`  | numeric, the scaling factor to use before the asinh transformation, a NULL value indicates different scaling values per channel, default is 5

Output relations|.
---|---
`asinh`| numeric, output transformation per `row` and `col`.

##### Details

Values are scaled first and then asinh transformation is performed. One data point per cell is required as input. 
```r
asinh(value/scale)
```
A scale of 5 is recommended for cytof measurement and 150 for flowcyto measurements.

If a NULL scaling value is used then the scaling factor is give by the input cross-tab projection, the second factor defined in the `row` dimension indicates the scaling values to use per channel.

##### References

See the `base::asinh` [function of the R package for more information](https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions).
