This code is a MATLAB implementation of *chordal conversion* that uses
[MOSEK](https://www.mosek.com/) to solve a standard-form semidefinite
program (SDP) over a single $n\times n$ matrix variable and subject to
$m$ linear constraints:
$$\min_{X\succeq0}\left\{ \langle C,X\rangle:\mathrm{lb}_{i}\le\langle A_{i},X\rangle\le\mathrm{ub}_{i}\text{ for all }i\in\{1,2,\dots,m\}\right\} ,$$
where $X\succeq0$ indicates that the matrix variable $X$ should be
symmetric positive semidefinite.

The problem data are the objective matrix $C\in\mathbb{S}^{n}$, the list
of $m$ constraint matrices $A_{1},A_{2},\dots,A_{m}\in\mathbb{S}^{n}$
and associated lower-bounds
$\mathrm{lb}_{1},\mathrm{lb}_{2},\dots,\mathrm{lb}_{m}\in\mathbb{R}\cup\{-\infty\}$
and upper-bounds
$\mathrm{ub}_{1},\mathrm{ub}_{2},\dots,\mathrm{ub}_{m}\in\mathbb{R}\cup\{+\infty\}$,
where $\mathbb{S}^{n}$ is the set of $n\times n$ real symmetric
matrices. In this case, the matrix inner product is explicitly written
$\langle C,X\rangle=\mathrm{trace}(C^{T}X)$.

**Main feature.** Under suitable *small treewidth* assumptions (see
below), the code will reduce MOSEK's per-iteration cost to $O(m+n)$
time. Given that MOSEK usually converges to $\epsilon=10^{-8}$ in less
than ten iterations, the total cost of solving the SDP is also
essentially $O(m+n)$ time.

To give some benchmark numbers, the code is able to solve a Lovasz Theta
problem with $n=10^{6}$ and $m=1.5\times10^{6}$ in 30 minutes, and an AC
optimal power flow relaxation with $n=8.2\times10^{4}$ and
$m\approx2.5\times10^{5}$ in less than 4 hours, on a modest workstation
with a Xeon 3.3 GHz quad-core CPU and 32 GB of RAM.

**Complex data.** The code also natively supports complex-valued data
$C,A_{1},A_{2},\dots,A_{m}\in\mathbb{H}^{n}$, where $\mathbb{H}^{n}$ is
the set of $n\times n$ complex Hermitian matrices. The matrix inner
product is now
$\langle C,X\rangle=\mathrm{Re}\{\mathrm{trace}(\overline{C}^{T}X)\}$,
in which the overline indicates elementwise complex conjugation. (The
bounds $\mathrm{lb}_{i}$ and $\mathrm{ub}_{i}$ must still be
real-valued.)

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

-   `c` is an $n^{2}\times1$ sparse vector representing the column-wise
    vectorization $\mathrm{vec}(C)$. In MATLAB, this kind of
    vectorization can be implemented by calling `c=C(:). `

-   `At` is an $n^{2}\times m$ sparse matrix, whose $j$-th column
    represents the column-wise vectorization $\mathrm{vec}(A_{i})$ of
    the $i$-th constraint matrix. The $n\times n$ matrices
    $C,A_{1},\dots,A_{m}$ can be either real symmetric or complex
    Hermitian.

-   `lb` and `ub` are dense $m\times1$ vectors, whose $i$-th elements
    represent $\mathrm{lb}_{i}\in\mathbb{R}\cup\{-\infty\}$ and
    $\mathrm{ub}_{i}\in\mathbb{R}\cup\{+\infty\}$ respectively. All
    elements must satisfy $\mathrm{lb}_{i}\le\mathrm{ub}_{i}$.

The output is given as follows

-   `U` is the $n\times r$ dense factor matrix $U$, which recovers the
    SDP solution $X$ via the Hermitian outer product `X=UU’`. The matrix
    `U` may be real or complex. **Warning:** Explicitly forming $X$ will
    ruin the $O(m+n)$ time complexity of the algorithm. In almost all
    cases, operations with $X$ can be directly replaced by operations
    with $U$.

-   `v` is an $m\times2$ dense matrix, containing the dual Lagrange
    multipliers. These solve the following problem
    $$\max_{v^{+},v^{-}\ge0}\left\{ \langle\mathrm{lb},v^{-}\rangle+\langle\mathrm{ub},v^{+}\rangle:\sum_{i=1}^{m}(v_{i}^{+}-v_{i}^{-})A_{i}\preceq C\right\} .$$

-   `info` is a struct containing various timing information.

# The small treewidth assumption

A necessary condition for the code to be fast is that the aggregate
sparsity graph $G=(V,E)$ defined below should have a *small treewidth*:
$$\begin{gathered}
V=\{1,2,\dots,n\},\quad E=\mathrm{spar}(C)\cup\mathrm{spar}(A_{1})\cup\cdots\cup\mathrm{spar}(A_{m}),\\
\text{where }\mathrm{spar}(C)\equiv\{(i,j):C[i,j]\ne0\text{ for }i>j\}.
\end{gathered}$$ In other words, if $G$ does *not* have a small
treewidth, then the code is guaranteed *not* to be fast. (But the code
may be slow even if $G$ does have a small treewidth.)

In MATLAB, the adjacency matrix for the graph $G$ is easily computed.

    Adj = reshape(any([c,At],2),n,n); 

We can compute an upper-bound on $\mathrm{tw}(G)$ by performing symbolic
Cholesky factorization on the adjacency matrix.

    p = amd(Adj); % fill-reducing permutation
    twub = max(symbfact(Adj(p,p))-1; % symbolic Cholesky factor

Conversely, if $G$ has a small treewidth, then the code may be fast or
it may be slow. For example, if we set $C,A_{1},\dots,A_{m}$ to all be
diagonal matrices with dense diagonals, then the code will be slow.

A sufficient condition for the code to be fast is that the *extended*
aggregate sparsity graph $\overline{G}=(V,\overline{E})$ defined below
also needs to have a small treewidth: $$\begin{gathered}
V=\{1,2,\dots,n\},\quad\overline{E}=\mathrm{spar}(C)\cup\mathrm{clique}(A_{1})\cup\cdots\cup\mathrm{clique}(A_{m})\\
\text{where }\mathrm{clique}(A)=\{(i,j):A[i,k]\ne0\text{ or }A[k,j]\ne0\text{ for some }k\}.
\end{gathered}$$ In other words, if $\overline{G}$ has a small
treewidth, then the code is guaranteed to be fast. (But the code may be
fast even if $\overline{G}$ does not have a small treewidth.)

It turns out to be surprisingly difficult to efficiently compute the
adjancecy matrix of $\overline{G}$ once $n$ and $m$ become large. We
provide an efficient implementation below.

    [~,AdjExt] = sparseGraph(c,At);

We can then proceed to upper-bound $\mathrm{tw}(\overline{G})$ using the
same idea as above:

    p = amd(AdjExt); % fill-reducing permutation
    twub = max(symbfact(AdjExt(p,p))-1; % symbolic Cholesky factor

The rigorous guarantee is useful for applications where worst-case
runtime guarantee is required.

# Options structure

The code accepts an options structure for some advanced tuning.

    opt = struct('verbose',0);
    [U,v,info] = solveChordalConv(c,At,lb,ub,opt)

The fields of the options structure `opt` are:

-   'verbose' is set to 0 to turn off printouts.

-   'norecover' is set to TRUE if we only solve the SDP without
    recovering the minimizer $X$.

-   'perm' allows a custom elimination ordering to be provided in
    computing the tree decomposition.

-   'natural' forces MOSEK to use the natural elimination ordering,
    instead of the graph partition ordering.

-   'pfeas' overrides the MSK_DPAR_INTPNT_CO_TOL_PFEAS parameter in
    MOSEK.

# Known issues

MOSEK can sometimes terminate with "STATUS UNKNOWN" despite seemingly
achieving reasonable accuracy in its terminal output. This can usually
be fixed by lowering the primal feasibility tolerance, albeit at the
cost of a less accurate solution.

    opt = struct('pfeas',1e-3); % Reduce MOSEK pfeas tolerance
    [U,v] = solveChordalConv(c,At,lb,ub,opt);

The issue seems to occur at extremely large scales, i.e. with $n$ on the
order of 1 million or more, or with highly ill-conditioned data. (Full
credit to Martin S. Andersen for discovering this fix.)

# Acknowledgements

This code is developed by Richard Y. Zhang (ryz@illinois.edu). 
The author gratefully acknowledges Martin S. Andersen for
discussions and advice, and for early numerical experiments.