# MA106 Linear Algebra — Completed Solutions
**IIT Bombay, Spring 2021**
Solutions for all questions whose answers were missing or incomplete.

---

## Tutorial 1

### 1.6(iii) — General Solution of Ax = b (completed)

After row-reducing the augmented matrix for part (iii):

```
[0  2  -2  1 | 2]      swap R1↔R3      [1  3  0  1 | 8]
[2  -8  14 -5 | 2]    ─────────────►   [0 -14 14 -7 | -14]
[1   3   0  1 | 8]    then R2→R2-2R1   [0   0  0  0 |  0]
                       then 7R3+R2
```

From **row 2**: −14x₂ + 14x₃ − 7x₄ = −14 → x₂ = 1 + x₃ − x₄/2  
From **row 1**: x₁ = 8 − 3x₂ − x₄ = 5 − 3x₃ + x₄/2

Setting x₃ = s, x₄ = 2t (free variables):

**General solution:**
```
x = [5, 1, 0, 0]ᵀ  +  s[-3, 1, 1, 0]ᵀ  +  t[1, -1, 0, 2]ᵀ
```

---

## Tutorial 2

### 2.2 — Inverse of A by Gauss-Jordan Method

**A = [[1,0,0],[1,1,0],[1,1,1]]**

Set up [A | I] and perform row operations:

```
[1  0  0 | 1  0  0]
[1  1  0 | 0  1  0]   R2 → R2 − R1
[1  1  1 | 0  0  1]   R3 → R3 − R1

[1  0  0 | 1   0  0]
[0  1  0 | -1  1  0]   R3 → R3 − R2
[0  1  1 | -1  0  1]

[1  0  0 |  1   0  0]
[0  1  0 | -1   1  0]
[0  0  1 |  0  -1  1]
```

**A⁻¹ = [[1,0,0],[-1,1,0],[0,-1,1]]**

---

### 2.5 — Linear Independence

**(i)** {[1,−1,1], [3,5,2], [1,2,1], [1,1,1]} ⊂ ℝ¹ˣ³

We have **4 vectors in a 3-dimensional space**. By the dimension theorem, any set of more than 3 vectors in ℝ³ is automatically linearly **dependent**.

**(ii)** {[1,9,9,8], [2,0,0,3], [2,0,0,8]} ⊂ ℝ¹ˣ⁴

Row-reduce the matrix formed by these row vectors:

```
[1  9  9  8]            [1  9  9  8]
[2  0  0  3]  → R2-2R1 → [0 -18 -18 -13]
[2  0  0  8]    R3-2R1   [0 -18 -18  -8]
                → R3-R2: [0   0   0   5]
```

Rank = 3 = number of vectors, so they are **linearly independent**.

**(iii)** {[1,−1,0]ᵀ, [3,−5,2]ᵀ, [1,−2,1]ᵀ} ⊂ ℝ³ˣ¹

Row-reduce the matrix with these as columns:

```
[1   3   1]             [1  3  1]
[-1 -5  -2]  → R2+R1 →  [0 -2 -1]
[0   2   1]    (R3+R2): [0  0  0]
```

Rank = 2 < 3, so the set is **linearly dependent**.

---

### 2.7 — Ranks

**(i) A = [[8,−4],[−2,1],[6,−3]]**

Notice: Row 2 = (−1/4)·Row 1 and Row 3 = (3/4)·Row 1. All rows are multiples of [8,−4]. **Rank = 1**.

**(ii) A = [[0,8,−1],[1,2,0],[0,0,3],[0,4,5]]**

Swap R1 ↔ R2, then R4 → 2R4 − R2:

```
[1  2  0]
[0  8 -1]
[0  0  3]
[0  0 11]  → R4 - (11/3)R3 → [0  0  0]
```

Three pivot columns. **Rank = 3**.

---

### 2.8 — Subspaces of ℝ³ˣ¹

**(i) W = {[x₁,x₂,x₃]ᵀ : x₁+x₂+x₃ = 0}** — **YES, a subspace.**

- **Zero vector**: 0+0+0=0 ✓  
- **Closure under addition**: (x₁+y₁)+(x₂+y₂)+(x₃+y₃)=0 ✓  
- **Closure under scalar multiplication**: α(x₁+x₂+x₃)=0 ✓  
- **Basis**: {[1,−1,0]ᵀ, [1,0,−1]ᵀ}, **dim = 2**

**(ii) W = {[x₁+x₂+x₃, x₂+x₃, x₃]ᵀ : x₁,x₂,x₃ ∈ ℝ}** — **YES, a subspace.**

This is the column space of the invertible matrix [[1,1,1],[0,1,1],[0,0,1]] (upper triangular with det=1), so it equals all of ℝ³. **Basis**: {e₁,e₂,e₃}, **dim = 3**.

**(iii) W = {[x₁,x₂,x₃]ᵀ : x₁x₂x₃ = 0}** — **NOT a subspace.**

**Counterexample**: [1,1,0]ᵀ and [0,0,1]ᵀ are both in W, but their sum [1,1,1]ᵀ has x₁x₂x₃ = 1 ≠ 0. Not closed under addition.

**(iv) W = {[x₁,x₂,x₃]ᵀ : |x₁|,|x₂|,|x₃| ≤ 1}** — **NOT a subspace.**

**Counterexample**: [1,1,1]ᵀ ∈ W, but 2·[1,1,1]ᵀ = [2,2,2]ᵀ ∉ W. Not closed under scalar multiplication.

---

### 2.9 — Subspaces of ℝ, ℝ²ˣ¹, ℝ³ˣ¹, ℝ⁴ˣ¹

- **ℝ**: {0} and ℝ itself (dimension 0 or 1).
- **ℝ²ˣ¹**: {0} (a point), any **line through the origin** (dimension 1), ℝ² itself. *Geometrically*: a point, a line, or the full plane.
- **ℝ³ˣ¹**: {0}, **lines through the origin** (dim 1), **planes through the origin** (dim 2), ℝ³ itself. *Geometrically*: a point, a line, a plane, or 3-space.
- **ℝ⁴ˣ¹**: {0}, 1-dim subspaces, 2-dim subspaces, 3-dim subspaces (hyperplanes), ℝ⁴.

---

## Tutorial 3

### 3.1 — Spanning Sets and Dimension

Let V be a subspace with dim V = r, span S = V, |S| = s.

**(i) s ≥ r:** Since span S = V, we can extract a basis of V from S. A basis has exactly r elements. Since we can only extract, not add, we need s ≥ r.

**(ii) If s = r, S is a basis:** S spans V and has exactly r = dim V elements. A spanning set of size equal to dim V is necessarily linearly independent (otherwise we could reduce it further, contradicting span S = V with fewer vectors). Hence S is a basis.

**(iii) If s > r, S contains a basis:** Extract a maximal linearly independent subset T ⊆ S. Since span S = V, this T also spans V (any element of V reachable from S is reachable from T by removing redundant vectors). T has exactly r elements, so T is a basis for V.

---

### 3.2 — Pivotal Columns of REF Form a Basis for Column Space

Let A₀ be in REF with pivot columns at positions j₁ < j₂ < ... < jᵣ.

**Linear independence:** Each pivot column has a leading 1 in a row where all other pivot columns have 0. Hence no pivot column is a linear combination of the others.

**Spanning C(A₀):** Every non-pivot column of A₀ can be expressed as a linear combination of the pivot columns to its left (the entries above the pivot rows directly give the coefficients). Hence every column of A₀ lies in the span of the pivot columns.

---

### 3.3 — Row Space Equals Rank

**R(A) is a subspace:** 0 = 0·r₁+···+0·rₘ ∈ R(A). If u,v ∈ R(A), then u+v and αu are again linear combinations of rows, so R(A) is closed. ✓

**EROs preserve R(A):** Each ERO replaces a row by a linear combination of rows of A, so R(A') ⊆ R(A). Since EROs are invertible (each has an inverse ERO), R(A) ⊆ R(A'). Hence R(A') = R(A).

**dim R(A) = rank A:** The nonzero rows of the REF of A are linearly independent (each has a leading entry in a column where all others are 0) and they span R(A) (by the above). Their count equals the number of pivots = rank A.

---

### 3.4 — rank(AB) ≤ min{rank A, rank B}

**rank(AB) ≤ rank A:** Every column of AB is A times a column of B, so C(AB) ⊆ C(A). Hence rank(AB) = dim C(AB) ≤ dim C(A) = rank A.

**rank(AB) ≤ rank B:** Every row of AB is a row of A times B, so R(AB) ⊆ R(B). Hence rank(AB) = dim R(AB) ≤ dim R(B) = rank B.

---

### 3.5 — Rank, Nullity, and General Solution

**A = [[0,0,0,−2,1],[0,2,−2,14,−1],[0,2,3,13,1]]**

Row reduce: R3 → R3 − R2, then reorder pivots:

```
[0  2  -2  14  -1]   (pivot col 2)
[0  0   5  -1   2]   (pivot col 3)
[0  0   0  -2   1]   (pivot col 4)
```

**Rank = 3**, **Nullity = 5 − 3 = 2** (free variables: x₁, x₅).

**Null space (Ax = 0):** Free variables x₁ = s, x₅ = 10t:

From row 3: x₄ = x₅/2 = 5t  
From row 2: x₃ = −3x₅/10 = −3t  
From row 1: x₂ = −33x₅/10 = −33t

**Null space basis:** {[1,0,0,0,0]ᵀ, [0,−33,−3,5,10]ᵀ}

**General solution of Ax = b, b = [2,2,3]ᵀ:**

Augmented matrix reduces to:

```
[0  2  -2  14  -1 | 2]
[0  0   5  -1   2 | 1]
[0  0   0  -2   1 | 2]
```

From row 3: x₄ = (x₅−2)/2; From row 2: x₃ = −3x₅/10; From row 1: x₂ = 8 − 33x₅/10.

**Particular solution** (x₁=0, x₅=0): **x_p = [0, 8, 0, −1, 0]ᵀ**

**General solution:**
```
x = [0, 8, 0, -1, 0]ᵀ + s[1,0,0,0,0]ᵀ + t[0,-33,-3,5,10]ᵀ
```

---

### 3.9 — Rank via Determinants

**(i) A = [[0,2,−3],[2,0,5],[−3,5,0]]**

det(A) = 0(0−25) − 2(0+15) + (−3)(10−0) = 0 − 30 − 30 = **−60 ≠ 0**

Since det(A) ≠ 0, **rank = 3**.

*Verification by REF:* R1↔R2, then eliminate to get 3 pivots. ✓

**(ii) A = [[4,3],[−8,−6],[16,12]]**

Row 2 = −2·Row 1; Row 3 = 4·Row 1. Matrix has rank 1.

All 2×2 minors: det([[4,3],[−8,−6]])=−24+24=0, etc. All zero.

Since A ≠ O, **rank = 1**.

---

## Tutorial 4

### 4.1 — Cramer's Rule

System: x+2y+3z=20, x+3y+z=13, x+6y+αz=α.

**det(A)** = det([[1,2,3],[1,3,1],[1,6,α]])  
= 1(3α−6)−2(α−1)+3(6−3) = 3α−6−2α+2+9 = **α+5**

**Cramer's rule is applicable when α+5 ≠ 0, i.e., α ≠ −5.**

**When α = −5:**
Augmented matrix reduces (R2−R1, R3−R1, R3−4R2):

```
[1  2   3 | 20]
[0  1  -2 | -7]
[0  0   0 |  3]
```

Last row: 0 = 3. **No solution** when α = −5.

**When α ≠ −5** (Cramer's rule applies):

- det(Aₓ) [col 1 replaced by b] = 27α+114
- det(Aᵧ) [col 2 replaced by b] = −5α−19  
- det(A_z) [col 3 replaced by b] = α+8

**Unique solution:**
```
x = (27α+114)/(α+5),   y = (−5α−19)/(α+5),   z = (α+8)/(α+5)
```

---

### 4.2 — Cofactor Matrices and Inverses

**(i) A = [[a,b],[c,d]]**

Cofactor matrix: C = [[d,−c],[−b,a]]  
det(A) = ad−bc

**A⁻¹ = (1/(ad−bc))·[[d,−b],[−c,a]]**

**(ii) A = [[0,9,5],[2,0,0],[0,2,0]]**

det(A) = 5·det([[2,0],[0,2]]) = 5·4 = **20**

Cofactors:

| C₁₁=0 | C₁₂=0 | C₁₃=4 |
|--------|--------|--------|
| C₂₁=10 | C₂₂=0 | C₂₃=0 |
| C₃₁=0 | C₃₂=10 | C₃₃=−18 |

**A⁻¹ = (Cᵀ)/20 = [[0,1/2,0],[0,0,1/2],[1/5,0,−9/10]]**

**(iii) A = 3×3 Hilbert matrix [[1, 1/2, 1/3],[1/2, 1/3, 1/4],[1/3, 1/4, 1/5]]**

det(A) = **1/2160**

The cofactor matrix (A is symmetric so C = Cᵀ):

```
C₁₁ = 1/240,   C₁₂ = -1/60,   C₁₃ = 1/72
C₂₁ = -1/60,   C₂₂ = 4/45,    C₂₃ = -1/12
C₃₁ = 1/72,    C₃₂ = -1/12,   C₃₃ = 1/12
```

**A⁻¹ = 2160·C =**
```
[  9    -36    30 ]
[ -36   192  -180 ]
[  30  -180   180 ]
```

---

## Tutorial 5

### 5.1 — Eigenvalues and Diagonalizability

**(i) A = [[5,−1],[1,3]]**

Characteristic polynomial: det(A−λI) = (5−λ)(3−λ)+1 = λ²−8λ+16 = **(λ−4)²**

Eigenvalue λ=4, **algebraic multiplicity = 2**.

(A−4I) = [[1,−1],[1,−1]], rank=1, nullity=1. **Geometric multiplicity = 1**.

Since AM ≠ GM, **A is NOT diagonalizable**.

---

**(ii) A = [[3,2,1,0],[0,1,0,1],[0,2,−1,0],[0,0,0,1/2]]**

Expand det(A−λI) along column 1 (only top-left entry is nonzero there):

det = (3−λ)·det([[1−λ,0,1],[2,−1−λ,0],[0,0,1/2−λ]])

Inner det (expand along 3rd column): = (1/2−λ)·(1−λ)(−1−λ) = (1/2−λ)(λ²−1)

**Characteristic polynomial = (3−λ)(1/2−λ)(λ−1)(λ+1)**

**Eigenvalues: λ = 3, 1, −1, 1/2** — all distinct, so **A IS diagonalizable**.

Eigenvectors (solve (A−λI)x=0 for each):
- **λ=3**: v₁ = [1,0,0,0]ᵀ
- **λ=1**: v₂ = [−3,2,2,0]ᵀ
- **λ=−1**: v₃ = [1,0,−4,0]ᵀ
- **λ=1/2**: v₄ = [8,−6,−8,3]ᵀ

**P = [v₁|v₂|v₃|v₄]**, so P⁻¹AP = diag(3, 1, −1, 1/2).

---

**(iii) A = [[2,1,0],[0,2,1],[0,0,2]]**

Characteristic polynomial: **(2−λ)³** → eigenvalue λ=2, **AM=3**.

(A−2I) = [[0,1,0],[0,0,1],[0,0,0]], rank=2, nullity=1. **GM=1**.

Since AM ≠ GM, **A is NOT diagonalizable** (Jordan form has two Jordan blocks, but it is not diagonal).

---

## Tutorial 6

### 6.1 — Gram-Schmidt Orthonormalization

**(i) (e₁, e₁+e₂, e₁+e₂+e₃, e₁+e₂+e₃+e₄)**

Apply Gram-Schmidt:

- **u₁** = e₁ = [1,0,0,0]ᵀ
- x₂ = e₁+e₂; ⟨x₂,u₁⟩=1; x₂−u₁ = [0,1,0,0]ᵀ → **u₂** = [0,1,0,0]ᵀ = e₂
- x₃ = e₁+e₂+e₃; projections = 1,1; x₃−u₁−u₂ = [0,0,1,0]ᵀ → **u₃** = e₃
- x₄ = e₁+e₂+e₃+e₄; projections = 1,1,1; x₄−u₁−u₂−u₃ = [0,0,0,1]ᵀ → **u₄** = e₄

The Gram-Schmidt process on this set produces the **standard basis (e₁, e₂, e₃, e₄)**.

---

**(ii) (e₁+e₂+e₃+e₄, −e₁+e₂, −e₁+e₃, −e₁+e₄)**

Let x₁=[1,1,1,1]ᵀ, x₂=[−1,1,0,0]ᵀ, x₃=[−1,0,1,0]ᵀ, x₄=[−1,0,0,1]ᵀ.

**u₁** = x₁/‖x₁‖ = [1,1,1,1]ᵀ/2

⟨x₂,u₁⟩ = (−1/2+1/2+0+0) = 0 → x₂ already orthogonal.  
**u₂** = x₂/‖x₂‖ = [−1,1,0,0]ᵀ/√2

⟨x₃,u₁⟩ = 0; ⟨x₃,u₂⟩ = 1/√2  
x₃ − (1/√2)u₂ = [−1,0,1,0]ᵀ − [1/2,−1/2,0,0]ᵀ = [−1/2,−1/2,1,0]ᵀ  
‖ ‖= √(1/4+1/4+1) = √(3/2)  
**u₃** = [−1,−1,2,0]ᵀ/√6

⟨x₄,u₁⟩ = 0; ⟨x₄,u₂⟩ = 1/√2; ⟨x₄,u₃⟩ = 1/√6  
x₄ − (1/√2)u₂ − (1/√6)u₃ = [−1,0,0,1]ᵀ − [1/2,−1/2,0,0]ᵀ − [1/6,1/6,−2/6,0]ᵀ  
= [−1/3, −1/3, −1/3, 1]ᵀ, ‖ ‖ = 2/√3  
**u₄** = [−1,−1,−1,3]ᵀ/(2√3)

---

### 6.2 — Extending to u₄ (completing the solution)

We need u₄ = [α₁,α₂,α₃,α₄]ᵀ orthogonal to u₁,u₂,u₃ with ‖u₄‖=1.

Orthogonality to x₁=[1,−1,2,0]ᵀ: α₁−α₂+2α₃=0  
Orthogonality to x₂=[1,1,2,0]ᵀ: α₁+α₂+2α₃=0  
Orthogonality to x₃=[3,0,0,1]ᵀ: 3α₁+α₄=0

From equations 1&2: 2α₁+4α₃=0 → α₁=−2α₃ and −2α₂+4α₃=0... wait:  
Adding: 2α₁+4α₃=0 → α₁=−2α₃. Subtracting: −2α₂+0=0 → α₂=0... hmm, let me redo.

Equations: α₁−α₂+2α₃=0 and α₁+α₂+2α₃=0.  
Sum: 2α₁+4α₃=0 → α₁=−2α₃. Difference: −2α₂=0 → α₂=0.  
From eq 3: 3(−2α₃)+α₄=0 → α₄=6α₃. Set α₃=1: α=[−2,0,1,6]ᵀ.

**u₄ = [−2,0,1,6]ᵀ/√41**

**Express v=[1,−1,1,−1]ᵀ in this ONB:**

⟨u₁,v⟩ = (1+1+2+0)/√6 = 4/√6  
⟨u₂,v⟩ = (1−5+2+0)/√30 = −2/√30  
⟨u₃,v⟩ = (12+0+6(−1)+5(−1))... wait, using orthonormal vectors:  
u₁=[1,−1,2,0]ᵀ/√6: ⟨u₁,v⟩=(1+1+2+0)/√6=4/√6  
u₂=[1,5,2,0]ᵀ/√30: ⟨u₂,v⟩=(1−5+2+0)/√30=−2/√30  
u₃=[12,0,−6,5]ᵀ/√205: ⟨u₃,v⟩=(12+0−6−5)/√205=1/√205  
⟨u₄,v⟩=(−2+0+1−6)/√41=−7/√41

**v = (4/√6)u₁ − (2/√30)u₂ + (1/√205)u₃ − (7/√41)u₄**

---

### 6.3 — Unitary ↔ Rows Orthonormal

**A is unitary iff AA* = I.**

(AA*)_{ij} = Σₖ A_{ik}·Ā_{jk} = ⟨row_i of A, row_j of A⟩ (inner product in ℝ¹ˣⁿ)

This equals δ_{ij} for all i,j iff the rows form an orthonormal set. ∎

---

### 6.4 — M_F^E(I) is Unitary

The matrix M_F^E(I) has as its j-th column the coordinate vector of u_j with respect to E=(e₁,...,eₙ). Since uⱼ = Σₖ(uⱼ)ₖeₖ, the j-th column is simply the vector uⱼ in standard coordinates. Since {u₁,...,uₙ} is an ONB, the columns are orthonormal, so M_F^E(I) is unitary. ∎

---

### 6.5 — p(λ) is Eigenvalue of p(A)

If Ax = λx (x ≠ 0), then by induction Aᵏx = λᵏx for all k ≥ 0.

For any polynomial p(t) = Σ aₖtᵏ:
```
p(A)x = Σ aₖAᵏx = Σ aₖλᵏx = p(λ)·x
```
Hence p(λ) is an eigenvalue of p(A). ∎

---

### 6.6 — Eigenvalues of A given A³−6A²+11A = 6I

Consider p(λ) = λ³−6λ²+11λ−6 = (λ−1)(λ−2)(λ−3).

By the given equation, p(A)+6I−6I = 0 → **A³−6A²+11A−6I = O**, i.e., p(A) = O.

By Cayley-Hamilton (or minimal polynomial argument), every eigenvalue λ must satisfy p(λ)=0, so **λ ∈ {1, 2, 3}**.

det(A) = product of eigenvalues. Given 5 ≤ det(A) ≤ 7, and the only product of elements from {1,2,3} (with 3 eigenvalues) lying in [5,7] is **1·2·3 = 6**.

**Eigenvalues: λ₁=1, λ₂=2, λ₃=3.** Since they are distinct, **A is diagonalizable**.

---

### 6.7 — Spectral Decomposition A = Σλⱼuⱼuⱼ*

Since {u₁,...,uₙ} is an orthonormal basis, any x ∈ Kⁿ can be written x = Σⱼ⟨uⱼ,x⟩uⱼ.

Then:
```
Ax = Σⱼ ⟨uⱼ,x⟩ Auⱼ = Σⱼ ⟨uⱼ,x⟩ λⱼuⱼ = Σⱼ λⱼ(uⱼuⱼ*)x
```
Since this holds for all x: **A = Σⱼ λⱼuⱼuⱼ*** ∎

---

### 6.8 — Unitary and Skew Self-Adjoint Eigenvalues

**(i)** λ is an eigenvalue of A iff λ̄ is an eigenvalue of A*.
(Already proved in 5.4 using det(B)=0 ↔ det(B*)=0 with B=A−λI.)

**(ii)** A unitary: ‖Ax‖²=⟨Ax,Ax⟩=⟨x,A*Ax⟩=⟨x,Ix⟩=‖x‖². If Ax=λx: |λ|‖x‖=‖Ax‖=‖x‖ → **|λ|=1**.

**(iii)** A skew self-adjoint (A*=−A). If Ax=λx, then:
```
⟨Ax,x⟩ = λ‖x‖² and ⟨Ax,x⟩ = ⟨x,A*x⟩ = −⟨x,Ax⟩ = −λ̄‖x‖²
```
So λ = −λ̄, i.e., λ+λ̄=0, meaning Re(λ)=0. Writing λ=iy (y∈ℝ), we get **iλ = i(iy) = −y ∈ ℝ**. ∎

---

### 6.9 — Normal iff Σ|aⱼₖ|² = Σ|λⱼ|²

Note Σ|aⱼₖ|² = tr(A*A) and Σ|λⱼ|² = tr(D*D) where D=diag(λ₁,...,λₙ).

By Schur's decomposition: A = UTU* with U unitary and T upper triangular (diagonal entries = eigenvalues). Then tr(A*A) = tr(T*T) = Σⱼ|tⱼⱼ|² + Σⱼ≠ₖ|tⱼₖ|² = Σ|λⱼ|² + Σⱼ<ₖ|tⱼₖ|².

Equality holds **iff all off-diagonal entries of T are zero, iff T is diagonal, iff A is unitarily diagonalizable, iff A is normal**. ∎

---

### 6.10 — Nilpotent Matrices

**Part 1:** If A is n×n upper triangular with 0 diagonal, then (Aᵏ)_{ij}=0 whenever j−i<k (entries below the k-th superdiagonal are zero). For k=n, all entries have j−i<n, so **Aⁿ=O**.

**Part 2:** A∈Cⁿˣⁿ is nilpotent iff 0 is the only eigenvalue.

- (⟹) If Aᵐ=O and Ax=λx, then 0=Aᵐx=λᵐx, so λ=0.
- (⟸) If 0 is the only eigenvalue, the characteristic polynomial is (−λ)ⁿ=(−1)ⁿλⁿ. By Cayley-Hamilton, (−1)ⁿAⁿ=O, so **Aⁿ=O**. ∎

---

## Tutorial 7

### 7.1 — A is Self-Adjoint iff Normal + Real Eigenvalues

**(⟹)** If A=A*, then AA*=AA=A*A → **A is normal**. For any eigenpair (λ,x): λ⟨x,x⟩=⟨Ax,x⟩=⟨x,A*x⟩=⟨x,Ax⟩=λ̄⟨x,x⟩, so **λ=λ̄ ∈ ℝ**.

**(⟸)** By spectral theorem for normal matrices, A=UDU* with D diagonal. If all diagonal entries of D are real, then D*=D̄=D, hence **A*=(UDU*)*=UD*U*=UDU*=A**. ∎

---

### 7.2 — Spectral Theorem for Skew Self-Adjoint

**Theorem:** If A∈Cⁿˣⁿ with A*=−A (skew self-adjoint), then:
1. A is unitarily diagonalizable.
2. All eigenvalues of A are purely imaginary.

**Proof:** Let B=iA. Then B*=(iA)*=−iA*=−i(−A)=iA=B. So B is self-adjoint. By the spectral theorem for self-adjoint matrices, B=UDU* where U is unitary and D=diag(d₁,...,dₙ) with each dⱼ∈ℝ.

Hence A=−iB=U(−iD)U*, which is a unitary diagonalization of A with diagonal entries **λⱼ=−idⱼ** (purely imaginary). ∎

---

### 7.3 — Orthonormal Eigenvectors and Spectral Representation

**A = diag(−1,−1) ⊕ [[−1,−4],[−4,−1]]**

**Top-left 2×2 block:** Eigenvalue λ=−1 (multiplicity 2), eigenvectors u₁=e₁, u₂=e₂.

**Bottom-right 2×2 block B = [[−1,−4],[−4,−1]]:**

char poly: (−1−λ)²−16=0 → λ²+2λ−15=(λ+5)(λ−3)=0

- **λ₃=−5:** (B+5I)x=0 → 4x₁−4x₂=0 → x₁=x₂. Eigenvector: **u₃=[0,0,1,1]ᵀ/√2**
- **λ₄=3:** (B−3I)x=0 → −4x₁−4x₂=0 → x₁=−x₂. Eigenvector: **u₄=[0,0,1,−1]ᵀ/√2**

**Eigenvalues:** −1 (×2), −5, 3.

**Spectral representation:**
```
A = (−1)(u₁u₁* + u₂u₂*) + (−5)u₃u₃* + 3u₄u₄*
```

**Computing A⁷x for x=[1,2,3,4]ᵀ:**

⟨u₁,x⟩=1, ⟨u₂,x⟩=2, ⟨u₃,x⟩=7/√2, ⟨u₄,x⟩=−1/√2

A⁷ = (−1)⁷(u₁u₁*+u₂u₂*) + (−5)⁷u₃u₃* + 3⁷u₄u₄*
   = −(u₁u₁*+u₂u₂*) − 78125·u₃u₃* + 2187·u₄u₄*

A⁷x = −1·u₁ − 2·u₂ + (−78125)(7/√2)·u₃ + 2187·(−1/√2)·u₄

**First two components:** [−1,−2,0,0]ᵀ

**Last two components** (using B⁷[3,4]ᵀ):
```
(−5)⁷·(7/2)[1,1]ᵀ + 3⁷·(−1/2)[1,−1]ᵀ
= −78125·(7/2)[1,1]ᵀ − 2187·(1/2)[1,−1]ᵀ
= [−274531, −272344]ᵀ
```

**A⁷x = [−1, −2, −274531, −272344]ᵀ**

---

### 7.4 — Positive Definite iff All Eigenvalues Positive

By spectral theorem: A=UDU* with D=diag(λ₁,...,λₙ), U unitary.

For any x≠0, let y=U*x≠0. Then:
```
⟨Ax,x⟩ = ⟨UDU*x,x⟩ = ⟨Dy,y⟩ = Σⱼλⱼ|yⱼ|²
```

**(⟸)** If all λⱼ>0 and y≠0, then Σλⱼ|yⱼ|²>0 since at least one |yⱼ|²>0. So A is positive definite.

**(⟹)** If A positive definite, take x=uⱼ (eigenvector, unit): ⟨Auⱼ,uⱼ⟩=λⱼ>0. ∎

---

### 7.5 — Corner Replacement Problem

Let α=[α₁,α₂,α₃,α₄]ᵀ. The replacement rule gives **β=Aα** where:
```
A = [[0,1/2,0,1/2],[1/2,0,1/2,0],[0,1/2,0,1/2],[1/2,0,1/2,0]]
```
After k steps: αₖ = Aᵏα.

**Finding orthonormal eigenvectors of A:**

- **λ₁=1**: u₁=[1,1,1,1]ᵀ/2
- **λ₂=0**: u₂=[1,0,−1,0]ᵀ/√2
- **λ₃=−1**: u₃=[1,−1,1,−1]ᵀ/2
- **λ₄=0**: u₄=[0,1,0,−1]ᵀ/√2

**Spectral representation:** A = 1·u₁u₁* + (−1)·u₃u₃*

By spectral theorem:
```
Aᵏ = 1ᵏ·u₁u₁* + (-1)ᵏ·u₃u₃*
```

Let S=α₁+α₂+α₃+α₄ and D=α₁−α₂+α₃−α₄. Then:

**After k steps:**
```
αₖ = (S/4)[1,1,1,1]ᵀ + (-1)ᵏ(D/4)[1,-1,1,-1]ᵀ
```

Explicitly:
- (αₖ)₁ = (αₖ)₃ = (S + (−1)ᵏD)/4
- (αₖ)₂ = (αₖ)₄ = (S − (−1)ᵏD)/4

As k→∞, all values converge to the **average (α₁+α₂+α₃+α₄)/4**.

---

### 7.6 — Lagrange Multiplier and Quadratic Forms

**Setup:** Minimize/maximize Q(x)=xᵀAx subject to g(x)=‖x‖²−1=0.

By the Lagrange multiplier condition: ∇Q = λ₀∇g.

Since Q(x) = xᵀAx with A symmetric: ∇Q = 2Ax. And ∇g = 2x.

So **2Ax = 2λ₀x**, i.e., **Ax₀ = λ₀x₀**. Hence x₀ is a unit eigenvector of A. ∎

Moreover, Q(x₀) = x₀ᵀ(Ax₀) = x₀ᵀ(λ₀x₀) = λ₀‖x₀‖² = λ₀.

So the **constrained extrema of Q equal the eigenvalues of A**:
- Maximum of Q on unit sphere = largest eigenvalue
- Minimum of Q on unit sphere = smallest eigenvalue ∎

---

## Tutorial 8

### 8.2 — Subspaces of Polynomial Space V

**(i) W₁ = {p∈V : p(0)=0}** — **YES, a subspace.**

p(0)=a₀=0. So W₁ = span{t, t², ..., tⁿ}. **Basis: {t, t², ..., tⁿ}, dim = n.**

**(ii) W₂ = {p∈V : p'(0)=0 = p''(0)}** — **YES, a subspace.**

p'(0)=a₁=0 and p''(0)=2a₂=0, so a₁=a₂=0.  
W₂ = {a₀ + a₃t³ + a₄t⁴ + ··· + aₙtⁿ}.  
**Basis: {1, t³, t⁴, ..., tⁿ}, dim = n−1.**

**(iii) W₃ = {p∈V : p odd, i.e., p(−t)=−p(t)}** — **YES, a subspace.**

An odd polynomial has only odd-degree terms.  
**Basis: {t, t³, t⁵, ...} (up to degree n), dim = ⌊n/2⌋.**

---

### 8.3 — Linear Independence of Trigonometric Functions

**S₁ = {1, cos, sin}:** Suppose α₁·1 + α₂·cos(t) + α₃·sin(t) = 0 for all t.

- t=0: α₁ + α₂ = 0
- t=π: α₁ − α₂ = 0  
- From these two: α₁=0, α₂=0.
- t=π/2: α₁ + α₃ = 0 → α₃=0.

All coefficients zero, so **S₁ is linearly independent**.

**S₂ = {1, cos², sin²}:** Note that **1·cos²(t) + 1·sin²(t) = 1** for all t, i.e.:

1·(1) − 1·cos² − 1·sin² = 0 (coefficient of 1 is 1, others are −1).

This is a nontrivial relation, so **S₂ is linearly dependent**.

---

### 8.4 — Multiple Representations in R¹ˣ²

{v₁,v₂,v₃} = {[1,0],[1,1],[1,−1]} contains 3 vectors in ℝ² (dim=2), so they are **linearly dependent**. Indeed: 2v₁ = v₂+v₃.

Since there is a nontrivial null combination (2v₁−v₂−v₃=0), if 4v₁+16v₂+4v₃ = (24,12), we can add any multiple c(2v₁−v₂−v₃)=0:

(4+2c)v₁ + (16−c)v₂ + (4−c)v₃ = (24,12) for any c.

- c=0: 4v₁+16v₂+4v₃ ✓
- c=−2: 6v₁+15v₂+3v₃... let me check: (4−4)=0? No: c is added differently.

More directly: both representations are equal to (24,12) by direct computation, and they differ because the set is linearly dependent, permitting infinitely many representations.

---

### 8.5 — Dimensions of Matrix Subspaces (n×n real matrices)

| Subspace | Description | Dimension |
|----------|-------------|-----------|
| Diagonal matrices W₁ | n free diagonal entries | **n** |
| Upper triangular W₂ | diagonal + upper off-diagonal | **n(n+1)/2** |
| Symmetric W₃ | aᵢⱼ=aⱼᵢ, determined by upper triangle | **n(n+1)/2** |
| Skew-symmetric W₄ | aᵢⱼ=−aⱼᵢ, aᵢᵢ=0 | **n(n−1)/2** |

---

### 8.6 — V×W as a Vector Space

**Addition:** (v₁,w₁)+(v₂,w₂) = (v₁+v₂, w₁+w₂). Inherits commutativity, associativity from V and W. Zero element: (0_V, 0_W). Negatives: −(v,w)=(−v,−w).

**Scalar multiplication:** α(v,w)=(αv,αw). Distributivity and scalar axioms follow from V and W.

All eight vector space axioms are verified component-wise. ✓

**Dimension:** Take any basis {v₁,...,vₙ} of V and {w₁,...,wₘ} of W. The set {(v₁,0),...,(vₙ,0),(0,w₁),...,(0,wₘ)} is a basis for V×W.

**dim(V×W) = n+m = dim V + dim W**

---

### 8.7 — Matrix of T: K²ˣ² → K²ˣ²

The transformation T maps a 2×2 matrix X to the 2×2 matrix Y by **vec(Y) = A·vec(X)**, where vec stacks entries as [x₁₁,x₁₂,x₂₁,x₂₂]ᵀ.

The matrix of T with respect to the ordered basis {E₁₁,E₁₂,E₂₁,E₂₂} is computed by applying T to each basis element:

T(E₁₁): vec = [1,0,0,0]ᵀ → A[1,0,0,0]ᵀ = col 1 of A = [a₁₁,a₂₁,a₃₁,a₄₁]ᵀ

In general, the j-th column of the matrix of T is the j-th column of A. Therefore:

**The matrix of T with respect to this basis is A itself.**

---

### 8.8 — Matrix of T: P₂ → K²ˣ¹

T(α₀+α₁t+α₂t²) = [α₀+α₁, α₁+α₂]ᵀ

**With E=(1,t,t²) and F=(e₁,e₂):**

T(1) = [1,0]ᵀ → coefficients in F: [1,0]  
T(t) = [1,1]ᵀ → [1,1]  
T(t²) = [0,1]ᵀ → [0,1]

**M_F^E = [[1,1,0],[0,1,1]]**

**With E'=(1, 1+t, (1+t)²) and F'=(e₁, e₁+e₂):**

Express T of each basis element in terms of F'=(e₁, e₁+e₂): if T gives [a,b]ᵀ, write ae₁+b(e₁+e₂)... wait: αe₁+β(e₁+e₂)=[α+β,β]ᵀ. So α=a−b, β=b.

T(1) = [1,0]ᵀ: α=1,β=0 → [1,0] in F'  
T(1+t) = [2,1]ᵀ: α=1,β=1 → [1,1] in F'  
T((1+t)²) = T(1+2t+t²) = [1+2, 2+1]ᵀ=[3,3]ᵀ: α=0,β=3 → [0,3] in F'

**M_{F'}^{E'} = [[1,1,0],[0,1,3]]**

---

### 8.9 — Parallelogram Law in Inner Product Space

**Proof:**

‖v+w‖² = ⟨v+w, v+w⟩ = ‖v‖² + ⟨v,w⟩ + ⟨w,v⟩ + ‖w‖²

‖v−w‖² = ⟨v−w, v−w⟩ = ‖v‖² − ⟨v,w⟩ − ⟨w,v⟩ + ‖w‖²

Adding:  
**‖v+w‖² + ‖v−w‖² = 2‖v‖² + 2‖w‖²** ∎

---

### 8.10 — ⟨A,B⟩ = tr(A*B) is an Inner Product on Kᵐˣⁿ

**(i) Conjugate symmetry:** 
⟨B,A⟩ = tr(B*A). Since tr(B*A) = Σⱼ,ₖ B̄_{kj}A_{kj} = conj(Σⱼ,ₖ A*_{jk}B_{jk}) ... 

More directly: tr(B*A) = tr((A*B)*) (property of trace and conjugate transpose) = conj(tr(A*B)) = **⟨A,B⟩̄** ✓

**(ii) Linearity in second argument:**
⟨A, αB+γC⟩ = tr(A*(αB+γC)) = α·tr(A*B)+γ·tr(A*C) = α⟨A,B⟩+γ⟨A,C⟩ ✓

**(iii) Positive definiteness:**
⟨A,A⟩ = tr(A*A) = Σⱼ,ₖ|A_{jk}|² ≥ 0, with equality iff all A_{jk}=0 iff A=O. ✓

---

### 8.11 — Orthonormal Fourier Set in C([−π,π])

Use inner product ⟨f,g⟩ = ∫₋π^π f(t)g̅(t) dt.

**Normalization:**
- ∫₋π^π (1/2π)dt = 1, so ‖1/√(2π)‖=1 ✓  
- ∫₋π^π cos²(nt)/π dt = (1/π)·π = 1 ✓ (using ∫cos²(nt)dt=π for n≥1)  
- ∫₋π^π sin²(nt)/π dt = 1 ✓

**Orthogonality:** For m≠n (or between different types):
- ∫₋π^π cos(mt)cos(nt)dt = 0  
- ∫₋π^π sin(mt)sin(nt)dt = 0  
- ∫₋π^π cos(mt)sin(nt)dt = 0 (odd integrand, symmetric interval)  
- ∫₋π^π 1·cos(nt)dt = 0, ∫₋π^π 1·sin(nt)dt = 0

These are the **standard Fourier orthogonality relations**, establishing that the given set is orthonormal. ∎ (This is the foundation of Fourier series theory.)

---

### 8.12 — Properties of Hermitian Operators

Let T be Hermitian: T* = T (i.e., ⟨T(v),w⟩=⟨v,T(w)⟩ for all v,w).

**(i) ⟨T(v),v⟩ ∈ ℝ:**

⟨T(v),v⟩ = ⟨v,T*(v)⟩ = ⟨v,T(v)⟩ = **⟨T(v),v⟩̄**

A complex number equal to its own conjugate must be real. ∎

**(ii) Every eigenvalue is real:**

If T(v)=λv (v≠0), then λ‖v‖² = ⟨T(v),v⟩ ∈ ℝ (by part i). Since ‖v‖²>0, **λ∈ℝ**. ∎

**(iii) Eigenvectors for distinct eigenvalues are orthogonal:**

If T(v)=λv and T(w)=μw with λ≠μ:
```
λ⟨v,w⟩ = ⟨T(v),w⟩ = ⟨v,T(w)⟩ = μ̄⟨v,w⟩ = μ⟨v,w⟩
```
(using μ∈ℝ from part ii). So (λ−μ)⟨v,w⟩=0. Since λ≠μ, **⟨v,w⟩=0**. ∎

**(iv) T(W⁺) ⊆ W⊥ if T(W) ⊆ W:**

Let u ∈ W⊥. For any w ∈ W:
```
⟨T(u), w⟩ = ⟨u, T*(w)⟩ = ⟨u, T(w)⟩ = 0
```
The last equality holds because T(w) ∈ W (given) and u ∈ W⊥. Since ⟨T(u),w⟩=0 for all w∈W, we conclude **T(u) ∈ W⊥**. ∎

---
*End of Completed Solutions*
