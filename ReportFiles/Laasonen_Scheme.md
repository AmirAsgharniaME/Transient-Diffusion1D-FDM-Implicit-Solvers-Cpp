# 4.2.4 The Implicit Euler (Laasonen) Scheme

The Implicit Euler scheme, commonly known as the Laasonen scheme, is an implicit finite difference method for solving the one-dimensional transient diffusion equation.

Unlike the FTCS and Richardson schemes, which evaluate the spatial derivative using only previously known values, the Laasonen scheme evaluates the spatial derivative at the new time level $n+1$. Consequently, the unknown values at the new time level appear simultaneously in the discretized equation, and a system of linear algebraic equations must be solved at every time step.

The governing one-dimensional diffusion equation is given by [Equation (8)](#EQ08-1D-DIFFUSION-GOVERNING). Evaluating this equation at the grid point $(y_j,t^{n+1})$ gives:

<a id="EQ67-LAASONEN-DISCRETE-NOTATION"></a>

$$
\left(\frac{\partial u}{\partial t}\right)_j^{n+1}
=
\nu
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
\tag{67}
$$

The Laasonen scheme uses:

- a first-order backward difference in time;
- a second-order central difference in space;
- an implicit evaluation of the diffusion term at time level $n+1$;
- a tridiagonal linear system at each time step.

## A. Discretization of the First-Order Temporal Derivative

To approximate the first-order temporal derivative at time level $n+1$, the value $u_j^n$ is expanded about the point $(y_j,t^{n+1})$:

<a id="EQ68-LAASONEN-BACKWARD-TEMPORAL-TAYLOR"></a>

$$
u_j^n
=
u_j^{n+1}
-
\left(\frac{\partial u}{\partial t}\right)_j^{n+1}\Delta t
+
\frac{1}{2!}
\left(\frac{\partial^2 u}{\partial t^2}\right)_j^{n+1}
(\Delta t)^2
-
\frac{1}{3!}
\left(\frac{\partial^3 u}{\partial t^3}\right)_j^{n+1}
(\Delta t)^3
+
\mathcal{O}\left((\Delta t)^4\right)
\tag{68}
$$

Rearranging [Equation (68)](#EQ68-LAASONEN-BACKWARD-TEMPORAL-TAYLOR) to isolate the first-order temporal derivative gives:

<a id="EQ69-LAASONEN-BACKWARD-TIME"></a>

$$
\left(\frac{\partial u}{\partial t}\right)_j^{n+1}
=
\frac{
u_j^{n+1}
-
u_j^n
}{
\Delta t
}
+
\mathcal{O}(\Delta t)
\tag{69}
$$

Therefore, the backward difference approximation is first-order accurate in time.

## B. Discretization of the Second-Order Spatial Derivative

The second-order spatial derivative is evaluated at time level $n+1$. Expanding the neighboring values $u_{j+1}^{n+1}$ and $u_{j-1}^{n+1}$ about $(y_j,t^{n+1})$ gives:

<a id="EQ70-LAASONEN-FORWARD-SPATIAL-TAYLOR"></a>

$$
u_{j+1}^{n+1}
=
u_j^{n+1}
+
\left(\frac{\partial u}{\partial y}\right)_j^{n+1}\Delta y
+
\frac{1}{2!}
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
(\Delta y)^2
+
\frac{1}{3!}
\left(\frac{\partial^3 u}{\partial y^3}\right)_j^{n+1}
(\Delta y)^3
+
\frac{1}{4!}
\left(\frac{\partial^4 u}{\partial y^4}\right)_j^{n+1}
(\Delta y)^4
+
\cdots
\tag{70}
$$

<a id="EQ71-LAASONEN-BACKWARD-SPATIAL-TAYLOR"></a>

$$
u_{j-1}^{n+1}
=
u_j^{n+1}
-
\left(\frac{\partial u}{\partial y}\right)_j^{n+1}\Delta y
+
\frac{1}{2!}
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
(\Delta y)^2
-
\frac{1}{3!}
\left(\frac{\partial^3 u}{\partial y^3}\right)_j^{n+1}
(\Delta y)^3
+
\frac{1}{4!}
\left(\frac{\partial^4 u}{\partial y^4}\right)_j^{n+1}
(\Delta y)^4
+
\cdots
\tag{71}
$$

Adding [Equations (70)](#EQ70-LAASONEN-FORWARD-SPATIAL-TAYLOR) and [(71)](#EQ71-LAASONEN-BACKWARD-SPATIAL-TAYLOR) eliminates the odd-order spatial derivative terms:

<a id="EQ72-LAASONEN-CENTRAL-TAYLOR-SUM"></a>

$$
u_{j+1}^{n+1}
+
u_{j-1}^{n+1}
=
2u_j^{n+1}
+
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
(\Delta y)^2
+
\frac{1}{12}
\left(\frac{\partial^4 u}{\partial y^4}\right)_j^{n+1}
(\Delta y)^4
+
\cdots
\tag{72}
$$

Rearranging [Equation (72)](#EQ72-LAASONEN-CENTRAL-TAYLOR-SUM) gives the second-order central difference approximation:

<a id="EQ73-LAASONEN-CENTRAL-SPACE"></a>

$$
\left(\frac{\partial^2 u}{\partial y^2}\right)_j^{n+1}
=
\frac{
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
}{
(\Delta y)^2
}
+
\mathcal{O}\left((\Delta y)^2\right)
\tag{73}
$$

Thus, the spatial discretization is second-order accurate.

## C. Discrete Algebraic Equation

Substituting the temporal approximation from [Equation (69)](#EQ69-LAASONEN-BACKWARD-TIME) and the spatial approximation from [Equation (73)](#EQ73-LAASONEN-CENTRAL-SPACE) into the governing equation in [Equation (67)](#EQ67-LAASONEN-DISCRETE-NOTATION) gives:

<a id="EQ74-LAASONEN-DISCRETE-SUBSTITUTION"></a>

$$
\frac{
u_j^{n+1}
-
u_j^n
}{
\Delta t
}
=
\nu
\left(
\frac{
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
}{
(\Delta y)^2
}
\right)
\tag{74}
$$

Define the dimensionless diffusion number as:

<a id="EQ75-LAASONEN-DIFFUSION-NUMBER"></a>

$$
d
=
\frac{\nu\Delta t}{(\Delta y)^2}
\tag{75}
$$

Multiplying [Equation (74)](#EQ74-LAASONEN-DISCRETE-SUBSTITUTION) by $\Delta t$ and using [Equation (75)](#EQ75-LAASONEN-DIFFUSION-NUMBER) yields:

<a id="EQ76-LAASONEN-DIMENSIONLESS-FORM"></a>

$$
u_j^{n+1}
-
u_j^n
=
d
\left(
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
\right)
\tag{76}
$$

Rearranging [Equation (76)](#EQ76-LAASONEN-DIMENSIONLESS-FORM) so that all unknown values at time level $n+1$ are placed on the left-hand side gives:

<a id="EQ77-LAASONEN-TRIDIAGONAL-RELATION"></a>

$$
-d\,u_{j-1}^{n+1}
+
(1+2d)\,u_j^{n+1}
-
d\,u_{j+1}^{n+1}
=
u_j^n
\tag{77}
$$

[Equation (77)](#EQ77-LAASONEN-TRIDIAGONAL-RELATION) is the discrete algebraic form of the Laasonen scheme. It must be solved simultaneously for all spatial nodes at every time step.

## D. Accuracy of the Laasonen Scheme

The temporal truncation error is obtained from [Equation (69)](#EQ69-LAASONEN-BACKWARD-TIME):

$$
\text{Temporal truncation error}
=
\mathcal{O}(\Delta t)
$$

The spatial truncation error is obtained from [Equation (73)](#EQ73-LAASONEN-CENTRAL-SPACE):

$$
\text{Spatial truncation error}
=
\mathcal{O}\left((\Delta y)^2\right)
$$

Therefore, the Laasonen scheme is:

- first-order accurate in time;
- second-order accurate in space.

The combined local truncation error can be written as:

<a id="EQ78-LAASONEN-LOCAL-TRUNCATION-ERROR"></a>

$$
\mathcal{T}_{\mathrm{Laasonen}}
=
\mathcal{O}(\Delta t)
+
\mathcal{O}\left((\Delta y)^2\right)
\tag{78}
$$

The spatial accuracy of Laasonen is the same as that of FTCS, while its temporal accuracy is the same order as FTCS. The main difference between the two methods is their stability behavior and algebraic structure.




## G. Stability Analysis of the Laasonen Scheme

The Laasonen scheme is unconditionally stable for the one-dimensional diffusion equation. This means that its stability does not require a restriction of the form $d\leq\frac{1}{2}$, unlike the explicit FTCS scheme.

To demonstrate this property, consider a Fourier error mode:

<a id="EQ82-LAASONEN-FOURIER-ERROR-MODE"></a>

$$
\epsilon_j^n
=
G^n e^{ij\theta}
\tag{82}
$$

where $G$ is the amplification factor, $\theta$ is the phase angle, and $i=\sqrt{-1}$.

Substituting the Fourier mode from [Equation (82)](#EQ82-LAASONEN-FOURIER-ERROR-MODE) into the Laasonen discretization in [Equation (77)](#EQ77-LAASONEN-TRIDIAGONAL-RELATION) gives:

<a id="EQ83-LAASONEN-FOURIER-SUBSTITUTION"></a>

$$
\frac{
G e^{ij\theta}
-
e^{ij\theta}
}{
\Delta t
}
=
\nu
\frac{
G e^{i(j+1)\theta}
-
2G e^{ij\theta}
+
G e^{i(j-1)\theta}
}{
(\Delta y)^2
}
\tag{83}
$$

Dividing [Equation (83)](#EQ83-LAASONEN-FOURIER-SUBSTITUTION) by $e^{ij\theta}$ and using the diffusion number defined in [Equation (75)](#EQ75-LAASONEN-DIFFUSION-NUMBER) gives:

<a id="EQ84-LAASONEN-AMPLIFICATION-RELATION"></a>

$$
G-1
=
dG
\left(
e^{i\theta}
-
2
+
e^{-i\theta}
\right)
\tag{84}
$$

Using the identities $e^{i\theta}+e^{-i\theta}=2\cos\theta$ and $\cos\theta-1=-2\sin^2\left(\frac{\theta}{2}\right)$, [Equation (84)](#EQ84-LAASONEN-AMPLIFICATION-RELATION) becomes:

<a id="EQ85-LAASONEN-AMPLIFICATION-EQUATION"></a>

$$
G-1
=
-4dG\sin^2\left(\frac{\theta}{2}\right)
\tag{85}
$$

Rearranging [Equation (85)](#EQ85-LAASONEN-AMPLIFICATION-EQUATION) to isolate the amplification factor gives:

<a id="EQ86-LAASONEN-AMPLIFICATION-FACTOR"></a>

$$
G
=
\frac{
1
}{
1+4d\sin^2\left(\frac{\theta}{2}\right)
}
\tag{86}
$$

For every $d\geq0$ and every Fourier mode, the denominator in [Equation (86)](#EQ86-LAASONEN-AMPLIFICATION-FACTOR) is greater than or equal to one. Therefore:

<a id="EQ87-LAASONEN-STABILITY-CONDITION"></a>

$$
0<G\leq1
\tag{87}
$$

Consequently, all Fourier error modes remain bounded, and the Laasonen scheme is unconditionally stable.

This result contrasts with the explicit methods considered previously:

- FTCS is stable only when $d\leq\frac{1}{2}$.
- Richardson is unconditionally unstable.
- DuFort-Frankel is unconditionally stable.
- Laasonen is unconditionally stable.

Although Laasonen is unconditionally stable, unconditional stability does not imply unlimited accuracy. An excessively large time step may still produce significant temporal discretization errors because the Laasonen scheme is only first-order accurate in time.

## H. Comparison with the Explicit Schemes

The main characteristics of the four numerical schemes considered so far are summarized below:

| Scheme | Time discretization | Spatial discretization | Algebraic system | Temporal accuracy | Stability |
|---|---|---|---|---|---|
| FTCS | Forward explicit | Central | Not required | First-order | Conditionally stable, $d\leq\frac{1}{2}$ |
| Richardson | Central explicit | Central | Not required | Second-order | Unconditionally unstable |
| DuFort-Frankel | Central explicit | Modified central | Not required | Nominally second-order | Unconditionally stable |
| Laasonen | Backward implicit | Central | Tridiagonal system | First-order | Unconditionally stable |

The Laasonen scheme has the same formal temporal accuracy as the FTCS scheme, but it does not impose the FTCS stability restriction. Its main computational requirement is the solution of a tridiagonal linear system at every time step.

In contrast:

- FTCS is simple to implement because it requires only a direct explicit update.
- Richardson has second-order formal accuracy in both time and space, but it is unusable because of unconditional instability.
- DuFort-Frankel is explicit and unconditionally stable, but it uses three time levels and requires a starting procedure.
- Laasonen is implicit and unconditionally stable, but its solution requires a tridiagonal-system solver.

For the present transient Plane Couette Flow problem, Laasonen is a suitable method when stability and robustness are more important than minimizing the computational work per time step. Since its coefficient matrix is tridiagonal, the resulting system can be solved efficiently using a specialized tridiagonal solver, which will be described separately in a later section.

## I. Summary of the Laasonen Scheme

The main properties of the Implicit Euler (Laasonen) scheme are:

- It uses a first-order backward difference in time.
- It uses a second-order central difference in space.
- It is implicit because the unknown values at time level $n+1$ appear in the diffusion term.
- It produces a tridiagonal linear system at each time step.
- It is first-order accurate in time and second-order accurate in space.
- It is unconditionally stable for the one-dimensional diffusion equation.
- It does not require a special initialization step beyond the prescribed initial condition.
- It permits larger time steps than FTCS without violating the FTCS stability condition.
- Its linear system can be solved efficiently using a specialized tridiagonal solver.

Therefore, the Laasonen scheme provides a robust implicit alternative to the explicit FTCS, Richardson, and DuFort-Frankel schemes. The detailed solution procedure for the resulting tridiagonal system, including its implementation in the present C++ solver, is presented separately in the section dedicated to the Thomas Algorithm.
