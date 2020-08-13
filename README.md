# asinh operator

#### Description
`asinh` operator performs an inverse hyperbolic sine on values

##### Usage
Input projection|.
---|---
`y-axis` | values required to be transformed by the asinh operator
`row`    | channels, and scale values (optional)
`col`    | event, for example

Input parameters|.
---|---
`scale`             | numeric, the scaling factor to use before the asinh transformation, default is 5

Output relations|.
---|---
`asinh`| numeric, output transformation per `row` and `col`.

##### Details
Values are scaled first and then asinh tranfrmation is performed. One data point per cell is required as input. 
```r
asinh(value/scale)
```
A scale of 5 is recommended for cytof measurement and 150 for flowcyto measurements.

if a NULL scaling value is used then it scaling factor is found in the cross-tab projection, the second factor in the `row` dimension.


#### References
see the `base::asinh` function of the R package for the documentation,
https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions


##### See Also

#### Examples
