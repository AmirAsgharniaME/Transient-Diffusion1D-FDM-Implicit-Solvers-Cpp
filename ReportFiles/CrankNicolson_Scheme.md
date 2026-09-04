#### 4.2.5 The Crank-Nicolson Scheme

The Crank-Nicolson scheme is an implicit finite difference method that achieves second-order accuracy in time by evaluating the diffusion term at the arithmetic average of the current and future time levels. In contrast to the Laasonen scheme, which uses the diffusion term only at time level $n+1$, the Crank-Nicolson method blends the contributions from time levels $n$ and $n+1$.

For the one-dimensional transient diffusion equation, the Crank-Nicolson formulation can be written as the average of the spatial diffusion term at the two adjacent time levels:

<a id="EQ89-CRANK-NICOLSON-GOVERNING-FORM"></a>

$$
\left( \frac{\partial u}{\partial t} \right)_j^{n+\frac{1}{2}}
=
\frac{\nu}{2}
\left[
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n+1}
+
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n}
\right]
\tag{89}
$$

This midpoint form is the key distinction between Crank-Nicolson and Laasonen. Both methods are implicit, but Crank-Nicolson is more accurate in time.

## A. Time Derivative Discretization

To derive the finite difference approximation for the time derivative, expand $u_j^{n+1}$ and $u_j^n$ about the midpoint time level $t^{n+\frac{1}{2}}$.

Expanding $u_j^{n+1}$ about $t^{n+\frac{1}{2}}$ gives:

<a id="EQ90-CRANK-NICOLSON-TIME-TAYLOR-FORWARD"></a>

$$
u_j^{n+1}
=
u_j^{n+\frac{1}{2}}
+
\frac{\Delta t}{2}
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
+
\frac{(\Delta t)^2}{8}
\left(
\frac{\partial^2 u}{\partial t^2}
\right)_j^{n+\frac{1}{2}}
+
\mathcal{O}((\Delta t)^3)
\tag{90}
$$

Similarly, expanding $u_j^n$ about $t^{n+\frac{1}{2}}$ gives:

<a id="EQ91-CRANK-NICOLSON-TIME-TAYLOR-BACKWARD"></a>

$$
u_j^{n}
=
u_j^{n+\frac{1}{2}}
-
\frac{\Delta t}{2}
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
+
\frac{(\Delta t)^2}{8}
\left(
\frac{\partial^2 u}{\partial t^2}
\right)_j^{n+\frac{1}{2}}
-
\mathcal{O}((\Delta t)^3)
\tag{91}
$$

Subtracting [Equation (91)](#EQ91-CRANK-NICOLSON-TIME-TAYLOR-BACKWARD) from [Equation (90)](#EQ90-CRANK-NICOLSON-TIME-TAYLOR-FORWARD) yields:

<a id="EQ92-CRANK-NICOLSON-TIME-SUBTRACTION"></a>

$$
u_j^{n+1} - u_j^n
=
\Delta t
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
+
\mathcal{O}((\Delta t)^3)
\tag{92}
$$

Therefore, the midpoint time derivative can be approximated as:

<a id="EQ93-CRANK-NICOLSON-TIME-DERIVATIVE"></a>

$$
\left(
\frac{\partial u}{\partial t}
\right)_j^{n+\frac{1}{2}}
=
\frac{u_j^{n+1} - u_j^n}{\Delta t}
+
\mathcal{O}((\Delta t)^2)
\tag{93}
$$

Thus, the Crank-Nicolson scheme is second-order accurate in time, whereas Laasonen is only first-order accurate in time.

## B. Spatial Derivative Discretization

For the diffusion term, the second derivative is approximated at both time levels $n$ and $n+1$ using the central difference formula.

At time level $n+1$:

<a id="EQ94-CRANK-NICOLSON-SPATIAL-NPLUS1"></a>

$$
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n+1}
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
\mathcal{O}((\Delta y)^2)
\tag{94}
$$

At time level $n$:

<a id="EQ95-CRANK-NICOLSON-SPATIAL-N"></a>

$$
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n}
=
\frac{
u_{j+1}^{n}
-
2u_j^{n}
+
u_{j-1}^{n}
}{
(\Delta y)^2
}
+
\mathcal{O}((\Delta y)^2)
\tag{95}
$$

Taking the arithmetic average of [Equations (94)](#EQ94-CRANK-NICOLSON-SPATIAL-NPLUS1) and [(95)](#EQ95-CRANK-NICOLSON-SPATIAL-N) gives:

<a id="EQ96-CRANK-NICOLSON-SPATIAL-AVERAGE"></a>

$$
\frac{1}{2}
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n+1}
+
\frac{1}{2}
\left(
\frac{\partial^2 u}{\partial y^2}
\right)_j^{n}
=
\frac{1}{2(\Delta y)^2}
\left[
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
+
u_{j+1}^{n}
-
2u_j^{n}
+
u_{j-1}^{n}
\right]
+
\mathcal{O}((\Delta y)^2)
\tag{96}
$$

This confirms that the Crank-Nicolson scheme retains second-order spatial accuracy, the same as Laasonen and FTCS.

## C. Discrete Algebraic Equation

Substituting [Equation (93)](#EQ93-CRANK-NICOLSON-TIME-DERIVATIVE) and [Equation (96)](#EQ96-CRANK-NICOLSON-SPATIAL-AVERAGE) into [Equation (89)](#EQ89-CRANK-NICOLSON-GOVERNING-FORM) gives:

<a id="EQ97-CRANK-NICOLSON-DISCRETE-FORM"></a>

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
\frac{1}{2(\Delta y)^2}
\left[
u_{j+1}^{n+1}
-
2u_j^{n+1}
+
u_{j-1}^{n+1}
+
u_{j+1}^{n}
-
2u_j^{n}
+
u_{j-1}^{n}
\right]
\tag{97}
$$

Defining the diffusion number as:

<a id="EQ98-CRANK-NICOLSON-DIFFUSION-NUMBER"></a>

$$
d
=
\frac{\nu \Delta t}{(\Delta y)^2}
\tag{98}
$$

the Crank-Nicolson scheme can be rearranged into the standard implicit form:

<a id="EQ99-CRANK-NICOLSON-FINAL-FORM"></a>

$$
-\frac{d}{2}u_{j-1}^{n+1}
+
(1+d)u_j^{n+1}
-
\frac{d}{2}u_{j+1}^{n+1}
=
\frac{d}{2}u_{j-1}^{n}
+
(1-d)u_j^{n}
+
\frac{d}{2}u_{j+1}^{n}
\tag{99}
$$

Equation (99) is the discrete algebraic equation for Crank-Nicolson. It has the same tridiagonal coupling pattern as Laasonen, but the right-hand side now also contains information from time level $n$.

## D. Tridiagonal Linear System

For an interior node $j$, [Equation (99)](#EQ99-CRANK-NICOLSON-FINAL-FORM) couples only three adjacent unknowns at time level $n+1$. Therefore, the coefficient matrix is tridiagonal.

For a domain with interior nodes $j=1,2,\dots,N-2$, the system can be written as:

<a id="EQ100-CRANK-NICOLSON-TRIDIAGONAL-SYSTEM"></a>

$$
\begin{bmatrix}
1+d & -\frac{d}{2} & 0 & \cdots & 0 \\
-\frac{d}{2} & 1+d & -\frac{d}{2} & \ddots & \vdots \\
0 & -\frac{d}{2} & 1+d & \ddots & 0 \\
\vdots & \ddots & \ddots & \ddots & -\frac{d}{2} \\
0 & \cdots & 0 & -\frac{d}{2} & 1+d
\end{bmatrix}
\begin{bmatrix}
u_1^{n+1} \\
u_2^{n+1} \\
u_3^{n+1} \\
\vdots \\
u_{N-2}^{n+1}
\end{bmatrix}
=
\begin{bmatrix}
\text{RHS}_1 \\
\text{RHS}_2 \\
\text{RHS}_3 \\
\vdots \\
\text{RHS}_{N-2}
\end{bmatrix}
\tag{100}
$$

where the right-hand side contains the known values from time level $n$ and the contributions from the boundary conditions.

As in the Laasonen method, the resulting matrix is tridiagonal and therefore suitable for solution by a tridiagonal solver. The Thomas algorithm is the standard choice, and its implementation can be described in a separate file dedicated to the linear solver.

## E. Incorporation of Dirichlet Boundary Conditions

For a domain with $N$ grid points, the boundary nodes are $j=0$ and $j=N-1$, while the interior nodes are $j=1,\ldots,N-2$.

There are two equivalent ways to impose Dirichlet boundary conditions.

### Approach A: Interior System

In the reduced interior formulation, the boundary values are excluded from the unknown vector and their contributions are moved to the right-hand side. This produces an $(N-2)\times(N-2)$ tridiagonal system.

This is mathematically compact and is often convenient for analysis.

### Approach B: Global $N\times N$ System

In the global formulation, all grid nodes are retained in the solution vector, including the boundary nodes. The boundary rows are replaced by identity equations:

<a id="EQ101-CRANK-NICOLSON-BOUNDARY-CONDITIONS"></a>

$$
u_0^{n+1}=U_{\mathrm{bottom}},
\qquad
u_{N-1}^{n+1}=U_{\mathrm{top}}
\tag{101}
$$

The corresponding global system is:

<a id="EQ102-CRANK-NICOLSON-GLOBAL-SYSTEM"></a>

$$
\begin{bmatrix}
1 & 0 & 0 & \cdots & 0 \\
-\frac{d}{2} & 1+d & -\frac{d}{2} & \cdots & 0 \\
0 & -\frac{d}{2} & 1+d & \ddots & \vdots \\
\vdots & \ddots & \ddots & \ddots & -\frac{d}{2} \\
0 & \cdots & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
u_0^{n+1} \\
u_1^{n+1} \\
u_2^{n+1} \\
\vdots \\
u_{N-1}^{n+1}
\end{bmatrix}
=
\begin{bmatrix}
U_{\mathrm{bottom}} \\
\text{RHS}_1 \\
\text{RHS}_2 \\
\vdots \\
U_{\mathrm{top}}
\end{bmatrix}
\tag{102}
$$

The omitted entries follow the same tridiagonal pattern. This global formulation is convenient for code implementation because the algebraic index and the physical grid index remain identical.

For the present C++ solver, the global $N\times N$ representation is preferred because it keeps the indexing uniform and reduces the chance of boundary-handling errors.

## F. Accuracy of the Crank-Nicolson Scheme

The temporal truncation error is:

$$
\mathcal{O}((\Delta t)^2)
$$

The spatial truncation error is:

$$
\mathcal{O}((\Delta y)^2)
$$

Therefore, the Crank-Nicolson scheme is:

- second-order accurate in time;
- second-order accurate in space.

The combined local truncation error is:

<a id="EQ103-CRANK-NICOLSON-LOCAL-ERROR"></a>

$$
\mathcal{T}_{\mathrm{CN}}
=
\mathcal{O}((\Delta t)^2)
+
\mathcal{O}((\Delta y)^2)
\tag{103}
$$

This is the main accuracy advantage of Crank-Nicolson over Laasonen, which is only first-order accurate in time.

## G. Stability Analysis

The Crank-Nicolson scheme is unconditionally stable for the one-dimensional diffusion equation under standard Fourier analysis. However, unlike Laasonen, Crank-Nicolson is only marginally dissipative and may generate oscillatory behavior for large time steps or sharp initial gradients.

To analyze stability, consider a Fourier error mode:

<a id="EQ104-CRANK-NICOLSON-FOURIER-MODE"></a>

$$
\epsilon_j^n = G^n e^{ij\theta}
\tag{104}
$$

Substituting this mode into the Crank-Nicolson discretization yields the amplification factor:

<a id="EQ105-CRANK-NICOLSON-AMPLIFICATION-FACTOR"></a>

$$
G
=
\frac{
1 - 2d\sin^2\left(\frac{\theta}{2}\right)
}{
1 + 2d\sin^2\left(\frac{\theta}{2}\right)
}
\tag{105}
$$

Since the denominator is always positive and the magnitude of $G$ satisfies:

<a id="EQ106-CRANK-NICOLSON-STABILITY-BOUND"></a>

$$
|G|\leq 1
\tag{106}
$$

the scheme is unconditionally stable.

However, when $d$ is large, $G$ can become negative for some wave numbers, which causes nonphysical oscillations in the numerical solution. This is why Crank-Nicolson, although unconditionally stable, is not always the best choice when monotonic decay is required.

## H. Comparison Between Crank-Nicolson and Laasonen

The two implicit schemes used in this report are Laasonen and Crank-Nicolson. Both methods lead to tridiagonal linear systems and both are unconditionally stable for the one-dimensional diffusion equation. Their principal difference lies in how the diffusion term is treated in time.

### Laasonen
- first-order accurate in time;
- second-order accurate in space;
- fully implicit at time level $n+1$;
- more strongly damping;
- numerically robust and monotone for diffusion-dominated problems.

### Crank-Nicolson
- second-order accurate in time;
- second-order accurate in space;
- averages the diffusion term between $n$ and $n+1$;
- less numerical damping than Laasonen;
- may show oscillations for large $d$.

In summary, Crank-Nicolson is more accurate than Laasonen, but Laasonen is often more diffusive and therefore more forgiving for strongly dissipative transient problems. If the goal is higher temporal accuracy, Crank-Nicolson is preferable. If the goal is stronger numerical smoothing and simpler behavior, Laasonen can be advantageous.

## I. Comparison of Implicit and Explicit Schemes

The numerical schemes discussed in this report can be divided into explicit and implicit classes.

### Explicit schemes
- FTCS
- Richardson
- DuFort-Frankel

### Implicit schemes
- Laasonen
- Crank-Nicolson

The main differences are:

| Property | Explicit schemes | Implicit schemes |
|---|---|---|
| Unknowns at new time level | Computed directly | Found by solving a linear system |
| Stability | Usually restricted or problematic | Typically unconditionally stable for diffusion |
| Computational cost per step | Low | Higher |
| Time-step flexibility | Limited by stability | Much larger |
| Implementation | Simpler | More involved |
| Accuracy control | Often lower or unstable | Can be higher and more robust |

For the present diffusion problem:

- FTCS is simple but conditionally stable.
- Richardson is formally attractive but unusable because it is unconditionally unstable.
- DuFort-Frankel is explicit and unconditionally stable, but it uses three time levels and requires a startup procedure.
- Laasonen is robust and unconditionally stable, but only first-order accurate in time.
- Crank-Nicolson is the most accurate among the schemes considered here, with second-order accuracy in both time and space, but it may require care when large time steps create oscillatory behavior.

Therefore, the implicit schemes provide better stability properties, while the explicit schemes provide simpler updates but weaker stability guarantees.

## J. Summary of the Crank-Nicolson Scheme

The main properties of the Crank-Nicolson scheme are:

- It is an implicit finite difference method.
- It evaluates the diffusion term at the average of time levels $n$ and $n+1$.
- It is second-order accurate in time.
- It is second-order accurate in space.
- It produces a tridiagonal linear system at every time step.
- It is unconditionally stable for the one-dimensional diffusion equation.
- It may produce oscillations for large diffusion numbers.
- It is more accurate in time than Laasonen.
- It requires a tridiagonal solver such as the Thomas algorithm.

Therefore, Crank-Nicolson provides a high-accuracy implicit alternative for the transient Plane Couette Flow problem, especially when accuracy is more important than the extra damping provided by Laasonen.
