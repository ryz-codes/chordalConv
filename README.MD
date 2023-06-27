This code is a MATLAB implementation of *chordal conversion* that uses
[MOSEK](https://www.mosek.com/) to solve a standard-form semidefinite
program (SDP) over a single *n* × *n* matrix variable and subject to *m*
linear constraints:

min<sub>*X* ≽ 0</sub>{⟨*C*,*X*⟩:lb<sub>*i*</sub>≤⟨*A*<sub>*i*</sub>,*X*⟩≤ub<sub>*i*</sub> for all *i*∈{1,2,…,*m*}},

where *X* ≽ 0 indicates that the matrix variable *X* should be symmetric
positive semidefinite.

The problem data are the objective matrix *C* ∈ 𝕊<sup>*n*</sup>, the
list of *m* constraint matrices
*A*<sub>1</sub>, *A*<sub>2</sub>, …, *A*<sub>*m*</sub> ∈ 𝕊<sup>*n*</sup>
and associated lower-bounds
lb<sub>1</sub>, lb<sub>2</sub>, …, lb<sub>*m*</sub> ∈ ℝ ∪ { − ∞} and
upper-bounds
ub<sub>1</sub>, ub<sub>2</sub>, …, ub<sub>*m*</sub> ∈ ℝ ∪ { + ∞}, where
𝕊<sup>*n*</sup> is the set of *n* × *n* real symmetric matrices. In this
case, the matrix inner product is explicitly written
⟨*C*, *X*⟩ = trace(*C*<sup>*T*</sup>*X*).

**Main feature.** Under suitable *small treewidth* assumptions (see
below), the code will reduce MOSEK’s per-iteration cost to *O*(*m*+*n*)
time. Given that MOSEK usually converges to *ϵ* = 10<sup>−8</sup> in
less than ten iterations, the total cost of solving the SDP is also
essentially *O*(*m*+*n*) time.

To give some benchmark numbers, the code is able to solve a Lovasz Theta
problem with *n* = 10<sup>6</sup> and *m* = 1.5 × 10<sup>6</sup> in 30
minutes, and an AC optimal power flow relaxation with
*n* = 8.2 × 10<sup>4</sup> and *m* ≈ 2.5 × 10<sup>5</sup> in less than 4
hours, on a modest workstation with a Xeon 3.3 GHz quad-core CPU and 32
GB of RAM.

**Complex data.** The code also natively supports complex-valued data
*C*, *A*<sub>1</sub>, *A*<sub>2</sub>, …, *A*<sub>*m*</sub> ∈ ℍ<sup>*n*</sup>,
where ℍ<sup>*n*</sup> is the set of *n* × *n* complex Hermitian
matrices. The matrix inner product is now
⟨*C*, *X*⟩ = Re{trace(*C*<sup>*H*</sup>*X*)}, in which the superscript
*H* indicates the conjugate transpose. (The bounds lb<sub>*i*</sub> and
ub<sub>*i*</sub> must still be real-valued.)

**Citation.** If you find this code helpful or use it in your projects,
please cite:

# Installation and usage notes

This code works with the commerical [MOSEK](https://www.mosek.com/)
solver ([a free academic license is available
here](https://www.mosek.com/products/academic-licenses/)). The base
command `mosekopt` must be added to the search path.

Once unzipped, the base calling sequence is the following.

    [U,v,info] = solveChordalConv(c,At,lb,ub)

The input is provided in SeDuMi format:

-   `c` is an *n*<sup>2</sup> × 1 sparse vector representing the
    column-wise vectorization vec(*C*). In MATLAB, this kind of
    vectorization can be implemented by calling `c=C(:). `

-   `At` is an *n*<sup>2</sup> × *m* sparse matrix, whose *j*-th column
    represents the column-wise vectorization vec(*A*<sub>*i*</sub>) of
    the *i*-th constraint matrix. The *n* × *n* matrices
    *C*, *A*<sub>1</sub>, …, *A*<sub>*m*</sub> can be either real
    symmetric or complex Hermitian.

-   `lb` and `ub` are dense *m* × 1 vectors, whose *i*-th elements
    represent lb<sub>*i*</sub> ∈ ℝ ∪ { − ∞} and
    ub<sub>*i*</sub> ∈ ℝ ∪ { + ∞} respectively. All elements must
    satisfy lb<sub>*i*</sub> ≤ ub<sub>*i*</sub>.

The output is given as follows

-   `U` is the *n* × *r* dense factor matrix *U*, which recovers the SDP
    solution *X* via the Hermitian outer product `X=UU’`. The matrix `U`
    may be real or complex. **Warning:** Explicitly forming *X* will
    ruin the *O*(*m*+*n*) time complexity of the algorithm. In almost
    all cases, operations with *X* can be directly replaced by
    operations with *U*.

-   `v` is an *m* × 2 dense matrix, containing the dual Lagrange
    multipliers. 

-   `info` is a struct containing various timing information.

# The small treewidth assumption

A necessary condition for the code to be fast is that the aggregate
sparsity graph *G* = (*V*,*E*) defined below should have a *small
treewidth*:

*V* = {1, 2, …, *n*},  *E* = spar(*C*) ∪ spar(*A*<sub>1</sub>) ∪ ⋯ ∪ spar(*A*<sub>*m*</sub>),

where
spar(*C*) ≡ {(*i*,*j*) : *C*\[*i*,*j*\] ≠ 0 for *i* \> *j*}.
In other words, if *G* does *not* have a small treewidth, then the code
is guaranteed *not* to be fast. (But the code may be slow even if *G*
does have a small treewidth.)

In MATLAB, the adjacency matrix for the graph *G* is easily computed.

    Adj = reshape(any([c,At],2),n,n); 

We can compute an upper-bound on tw(*G*) by performing symbolic Cholesky
factorization on the adjacency matrix.

    p = amd(Adj); % fill-reducing permutation
    twub = max(symbfact(Adj(p,p))-1; % symbolic Cholesky factor

Conversely, if *G* has a small treewidth, then the code may be fast or
it may be slow. For example, if we set
*C*, *A*<sub>1</sub>, …, *A*<sub>*m*</sub> to all be diagonal matrices
with dense diagonals, then the code will be slow.

A sufficient condition for the code to be fast is that the *extended*
aggregate sparsity graph *G*<sup>ext</sup> = (*V*,*E*<sup>ext</sup>)
defined below also needs to have a small treewidth:

*V* = {1, 2, …, *n*},  *E*<sup>ext</sup> = spar(*C*) ∪ clique(*A*<sub>1</sub>) ∪ ⋯ ∪ clique(*A*<sub>*m*</sub>)

where
clique(*A*) = {(*i*,*j*) : *A*\[*i*,*k*\] ≠ 0 or *A*\[*k*,*j*\] ≠ 0 for some *k*}.
In other words, if *G*<sup>ext</sup> has a small treewidth, then the
code is guaranteed to be fast. (But the code may be fast even if
*G*<sup>ext</sup> does not have a small treewidth.)

It turns out to be surprisingly difficult to efficiently compute the
adjancecy matrix of *G*<sup>ext</sup> once *n* and *m* become large. We
provide an efficient implementation below.

    [~,AdjExt] = sparseGraph(c,At);

We can then proceed to upper-bound tw(*G*<sup>ext</sup>) using the same
idea as above:

    p = amd(AdjExt); % fill-reducing permutation
    twub = max(symbfact(AdjExt(p,p))-1; % symbolic Cholesky factor

The rigorous guarantee is useful for applications where worst-case
runtime guarantee is required.

# Options structure

The code accepts an options structure for some advanced tuning.

    opt = struct('verbose',0);
    [U,v,info] = solveChordalConv(c,At,lb,ub,opt)

The fields of the options structure `opt` are:

-   ’verbose’ is set to 0 to turn off printouts.

-   ’norecover’ is set to TRUE if we only solve the SDP without
    recovering the minimizer *X*.

-   ’perm’ allows a custom elimination ordering to be provided in
    computing the tree decomposition.

-   ’natural’ forces MOSEK to use the natural elimination ordering,
    instead of the graph partition ordering.

-   ’pfeas’ overrides the MSK_DPAR_INTPNT_CO_TOL_PFEAS parameter in
    MOSEK.

# Known issues

MOSEK can sometimes terminate with “STATUS UNKNOWN” despite seemingly
achieving reasonable accuracy in its terminal output. This can usually
be fixed by lowering the primal feasibility tolerance, albeit at the
cost of a less accurate solution.

    opt = struct('pfeas',1e-3); % Reduce MOSEK pfeas tolerance
    [U,v] = solveChordalConv(c,At,lb,ub,opt);

The issue seems to occur at extremely large scales, i.e. with *n* on the
order of 1 million or more, or with highly ill-conditioned data. (Full
credit to Martin S. Andersen for discovering this fix.)

# Acknowledgements

This code is developed by Richard Y. Zhang (ryz@illinois.edu). The
author gratefully acknowledges Martin S. Andersen for discussions and
advice, and for early numerical experiments.
