# KleinianGenus2.jl

Introduction
------------

`KleinianGenus2.jl` is a library for computing [Kleinian hyperelliptic functions](https://arxiv.org/abs/solv-int/9603005) associated with curves of genus 2 written in [Julia](http://julialang.org/).

This is an implementation in Julia of the algorithm presented in the paper
[Computation of genus 2 Kleinian hyperelliptic functions via Richelot isogenies](https://arxiv.org/abs/2603.23188),
by M. Smirnov.

Usage
-----

The package provides types `Genus2WCurve` and `Genus2WCurveDivisor` along with constructors for them.

The type `Genus2WCurveDivisor` represents a divisor on a complex curve of genus 2 given in Weierstrass form (i.e. the curve with algebraic equation `y^2 = F(x) = 4x^5 + a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0`). It is a parametric type depending on a floating point type `T`. It has two fields `P1` and `P2`, which are both of the type `Array{Complex{T}}`. Each of these fields should be either an empty array, or an array `[x,y]`, where `x` and `y` satisfy the equation above. The empty array is used to denote that the corresponding point of the curve is at infinity. An object `D` of the type `Genus2WCurveDivisor` corresponds to the divisor `D.P1 + D.P2 - 2infty`.

The type `Genu2WCurve` represents an complex curve of genus 2 given in Weierstrass form (i.e. the curve with algebraic equation `y^2 = F(x) = 4x^5 + a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0`). It is a parametric type depending on an floating point type `T`. `Genus2WCurve` has the following fields:
* `F` -- the polynomial, defining equation of the curve (type -- `Polynomial{Complex{T}}`).
* `Roots` -- the roots of the polynomial `F` (type -- `Array{Complex{T}}`).
* `Periods` -- some of vectors that belong to the period lattice of the curve (type -- `Matrix{Complex{T}}`). It is a matrix of the size `2 x 2` or `2 x 4` containing a basis of a Lagrange subgroup in the period lattice or a whole basis in the period lattice respectively (see more detailed explanation below, in the description of the optional fields of the constructor -- `is_real` and `triple_of_discs`).
* `Eta` -- contains vectors that determine Weierstrass zeta function monodromy, corresponding to basis elements contained in `Periods` (type -- `Matrix{Complex{T}}`).
* `KleinianFunc` -- a function, that calculates the values of classical Kleinian functions at a given point. That is, `WeierstrassFunc(z)` is the tuple `(sigma(z), zeta_1(z), zeta_2(z), p_22(z), p_12(z), p_11(z))`.
* `KleinianFuncSigmaSquared` -- a function that calculates the tuple `(sigma^2(z), zeta_1(z), zeta_2(z), p_22(z), p_12(z), p_11(z))`. It is slightly faster than the previous and should be used, when there is no need in the value `sigma(z)`.
* `KleinianFuncWeight2` -- a function that calculates values of the [Kleinian functions of weight 2](https://link.springer.com/article/10.1134/S1995080225615334) `T(z) = (S(z), S_22(z), S_12(z), S_11(z))` and their first derivatives with respect to `z_1` and `z_2`. This is the core function, on which the functions `KleinianFunc` and `KleinianFuncSigmaSquared` are based. The function `T` is an entire analytic function, so `KleinianFuncWeight2` has a more stable behaviour near sigma-divisor (where functions zeta and p have poles). If the output of `KleinianFuncWeight2` suffices to one's purposes, it should be used instead of the foregoing functions.
* `JacobiInversion` -- a function that given a point `z` calculates the divisor `D` on the curve, which is mapped to `z` by the Abel map.
* `AbelMap` -- a function that given a divisor `D` on the curve calculates its image with respect to the Abel map.

The mandatory argument for the constructor is
* `data`, which determines the algebraic equation of a curve `y^2 = F(x)`. More precisely, `data` should be either a vector of roots of the polynomial `F`, or the vector of its coefficients, or the polynomial itself (i.e. an object of the type `Polynomial{...}`).

It should be noted that the type `T`, on which `Genus2WCurve` depends parametrically, is defined to be the type of the real part of the sum of the roots.

Also the constructor has optional arguments:
* `source`, which specifies what is contained in `data`. If `source == "Coefficients"`, then the constructor treats `data` as the vector of coefficients (listed in the standard order, from the 0th to 5th), if `source == "Roots"`, then the constructor treats `data` as the vector of roots (with length 5), and if `source == "Polynomial"`, then the constructor treats `data` as `F`. By default `source` is set to `"Coefficients"`.
* `n`, the upper bound for the number of iterations of Landen's transform. By default `n` is set to `15`.
* `e`, the threshold for stopping the iterations of constructing a Richelot isogenous curve: if the maximal distance in pairs of roots, which converge to a common limit, is smaller than `e`, then the iterations are stopped. By default `e` is set to `10*eps(T)`, where `T` is the parameter floating point type (see the definition of `T` above).
* `is_real`, the flag, which can be used to signify that all roots of `F` are real. In this case the constructor finds a whole basis in the period lattice (with corresponding values of `Eta`). The function `AbelMap`, given that `is_real = true` normalizes the output to lie in the (2 dimensional) fundamental paralellogram of the period lattice.
* `triple_of_discs`, the optional argument, that contains a tuple `(C1,C2,C3)` of discs on the Riemann sphere. More precisely, `Cj` is a tuple `(cj, rj)`, where  `cj` is a center of a disc, and `rj` is its radius. If `rj > 0`, then the disc is understood in the usual way, and if `rj < 0`, then `Cj` denotes the outer disc, i.e. the set `{z: |z - cj| > -rj}`. The discs are required to be disjoint, and each of the discs has to contain exactly two of the Weierstrass points of the curve (i.e. the roots of `F` and infinity). In other words, two of the discs should contain (each) two of the roots of `F` and be bounded discs in the complex plane, and the remaining disc should contain the remaining root of `F` and be unbounded. Such a triple is always required for the algorithm, and is constructed automatically using the algorithm presented in the paper [Any six points on the Riemann sphere can be split into three pairs by a triple of disjoint discs](https://arxiv.org/abs/2604.00351). This algorithm may behave unstable in certain configurations of roots of `F`, so we provided a way to specify the discs by hand. Moreover, it should be noted, that if the flag `is_real` is set to true, then the field `triple_of_discs` is ignored. Finally, if `is_real == false`, then the columns of the field `Periods` contains two vectors that are the integrals of the standard basis of holomorphic 1-forms on the curve with respect to cycles given by boundaries of the first two discs in the triple.

Examples
--------

Examples of constructing a `WCurve`:

```
g2 = 4.0
g3 = 0.0
C = WCurve([g2,g3]) #adding 'source = "Invariants"' does not affect the result
```
```
g2 = 4.0
g3 = 0.0
C = WCurve([g2,g3], rectangular = true) #imposing the rectangular normalization of the basis in period lattice
```
```
g2 = BigFloat(4.0)
g3 = 0.0
C = WCurve([g2,g3]) #returns a curve with functions that perform all calculations in BigFloat.
```
```
setprecision(5000)
g2 = BigFloat(4.0)
g3 = 0.0
C = WCurve([g2,g3], n = 30) #due to very high precision of the floating point arithmetic
#it may require more iterations of Landen's tranform to achieve machine precision
```
```
r = [-1,0,1]
C = WCurve(r, source = "Roots") #the same curve as in the above examples
```

Weierstrass functions calling examples:
```
#given that C is an object of WCurve type
z = 6 + 7im
WF = C.WeierstrassFunc(z)
w = C.AbelMap(WF[3], WF[4]) # should be equal to z modulo period lattice
WF = C.WeierstrassFunc(z, use_periodicity = true) # the same as the previous call, may be more stable if |z| is large
WF = C.WeierstrassFuncSigmaSquared(z) # the same as the previous two, but the first value (i.e. sigma) is squared
WF = C.WeierstrassFuncRaw(z) # calculates auxiliary functions S = sigma^2, R = sigma^2 p, and their derivatives
```
