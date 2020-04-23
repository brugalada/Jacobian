# Error Propagation

This two pieces of code are a translation into Python of a Fortran routine implemented by L. Lindegren in 1995.

## Jacobian
Returns the Jacobian of the transformation from ICRS to l,b,plx,U,V,W (6x6 matrix).   

    Input: x (iterable)
        x[0]: right ascention (equatorial) -> degrees
        x[1]: declination (equatorial) -> degrees
        x[2]: parallax -> mas
        x[3]: proper motion (mualpha*) -> mas/yr
        x[4]: proper motion (mudelta) -> mas/yr
        x[5]: radial velocity -> km/s
    Returns:
        Jacobian: ndarray[6,6]
    Original Fortran Code: @L Lindegren, Lund Observatory, 15 March 1995
    Adaptation to Python 3: @Pau Ramos, ICC-UB, 2017

## Jacobian_4d
Returns the Jacobian of the transformation from ICRS to l,b,plx,U,V,W but 
    ignoring positional errors (4x4 matrix). By doing so, we obtain a much
    faster routine. However, the 6x6 Cov. Matrix can be recovered by uncommenting
    the lines where the elements of the Jacobian are calculated, as well as 
    changing the sizes of such Jacobians at inicialization and including 'jac1' in
    the final matrix product chain.
    
    Input: x (iterable)
        x[0]: right ascention (equatorial) -> degrees
        x[1]: declination (equatorial) -> degrees
        x[2]: parallax -> mas
        x[3]: proper motion (mualpha*) -> mas/yr
        x[4]: proper motion (mudelta) -> mas/yr
        x[5]: radial velocity -> km/s
    Returns:
        Jacobian: ndarray[4,4]
    Original Fortran Code: @L Lindegren, Lund Observatory, 15 March 1995
    Adaptation to Python 3: @Pau Ramos, ICC-UB, 2017
    
## Jacobian_4d_tan
Returns the Jacobian of the transformation from ICRS to l,b,plx,vl,vb,vlos but 
    ignoring positional errors (4x4 matrix). By doing so, we obtain a much
    faster routine. However, the 6x6 Cov. Matrix can be recovered by uncommenting
    the lines where the elements of the Jacobian are calculated, as well as 
    changing the sizes of such Jacobians at inicialization and including 'jac1' in
    the final matrix product chain.
    
    Input: x (iterable)
        x[0]: right ascention (equatorial) -> degrees
        x[1]: declination (equatorial) -> degrees
        x[2]: parallax -> mas
        x[3]: proper motion (mualpha*) -> mas/yr
        x[4]: proper motion (mudelta) -> mas/yr
        x[5]: radial velocity -> km/s
    Returns:
        Jacobian: ndarray[4,4]
    Original Fortran Code: @L Lindegren, Lund Observatory, 15 March 1995
    Adaptation to Python 3: @Pau Ramos, ICC-UB, 2020
