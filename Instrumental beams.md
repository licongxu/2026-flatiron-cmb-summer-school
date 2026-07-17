# CMB Beam Systematics — Mathematical Derivation

This document derives, step by step, the mathematics behind the
`CMB_Beam_Systematics.ipynb` notebook. The notebook studies how
non-ideal telescope beams contaminate CMB maps: **T→P leakage**
(temperature leaking into polarization) and **E→B leakage** (E-mode
mixing into B-mode).

---

## 0. Set-up: maps, beams, and convolution

A CMB map is a 2D field on the flat sky. The temperature map is $T(\mathbf{x})$;
the Stokes polarization maps are $Q(\mathbf{x})$ and $U(\mathbf{x})$.
A detector does not measure the sky directly — it measures the sky
**convolved with its beam** $B(\mathbf{x})$:

$$
d(\mathbf{x}) = (B * T)(\mathbf{x}) = \int d^2\mathbf{x}'\, B(\mathbf{x}-\mathbf{x}')\, T(\mathbf{x}').
$$

In Fourier space a convolution becomes a product:

$$
d(\ell) = B(\ell)\, T(\ell),
$$

where $B(\ell)$ is the **beam transfer function** (the 2D
Fourier transform of the real-space beam). The data-reduction code assumes a
*model* beam $B_{\mathrm{assumed}}$; the real instrument has $B_{\mathrm{real}}$.
Everything that follows comes from the gap between these two.

---

## 1. Beam shapes and their transfer functions

### Circular Gaussian

$$
B_G(\mathbf{x}) = \frac{1}{2\pi\sigma^2}\exp\!\left(-\frac{|\mathbf{x}|^2}{2\sigma^2}\right),
\qquad \sigma = \frac{\mathrm{FWHM}}{2\sqrt{2\ln 2}}.
$$

Its Fourier transform is also a Gaussian:

$$
B_G(\boldsymbol{\ell}) = \exp\!\left(-\tfrac{1}{2}\sigma^2 \ell^2\right),
$$

a low-pass filter: it damps small scales (high $\ell$). This is the
"assumed" beam in most of the notebook.

### Elliptical Gaussian

$$
B_E(\mathbf{x}) \propto \exp\!\left[-\tfrac{1}{2}\left(
\frac{X_r^2}{\sigma_x^2} + \frac{Y_r^2}{\sigma_y^2}\right)\right],
$$

with $(X_r,Y_r)$ the coordinates rotated by an angle $\theta$, and
ellipticity

$$
\varepsilon = \frac{\mathrm{fwhm}_x - \mathrm{fwhm}_y}{\mathrm{fwhm}_x + \mathrm{fwhm}_y}.
$$

$B_E$ breaks azimuthal symmetry: an $x$-cut and a $y$-cut differ. Its
$B(\boldsymbol{\ell})$ is an *anisotropic* Gaussian in Fourier space — it
mixes different $\boldsymbol{\ell}$ directions. This anisotropy is the seed
of the beam-asymmetry leakage below.

### Airy disk

$$
B_A(r) \propto \left[\frac{2 J_1(\pi r / r_0)}{\pi r / r_0}\right]^2,
$$

the diffraction pattern of a circular aperture. It has **rings (sidelobes)**
that a Gaussian model misses entirely — the seed of the sidelobe leakage below.

---

## 2. T→P leakage: the general mechanism

A detector at scan angle $\psi$ obeys the (single-detector) data model

$$
d(\psi) = T + Q\cos 2\psi + U\sin 2\psi.
$$

The map-maker solves the per-pixel linear system
$\hat m = (A^\top A)^{-1} A^\top d$ with rows $A_i = [1,\cos 2\psi_i,\sin 2\psi_i]$.
For uniform angle coverage this reduces to projecting the timestream onto
$1,\cos 2\psi,\sin 2\psi$:

$$
\hat T = \frac{1}{n}\sum_i d_i,\qquad
\hat Q = \frac{2}{n}\sum_i d_i\cos 2\psi_i,\qquad
\hat U = \frac{2}{n}\sum_i d_i\sin 2\psi_i.
$$

With a **perfect circular beam** the detector really does see
$d_i = T + Q\cos 2\psi_i + U\sin 2\psi_i$, and the projection recovers
$(T,Q,U)$ exactly. Leakage appears only when the *real* beam differs from
the *assumed* beam, because then the timestream contains a temperature-shaped
residual the solver cannot place in $T$ (it has already fit $T$ with the
modelled beam), so it is absorbed into the $Q/U$ solution:

$$
\delta T \equiv (B_{\rm real} - B_{\rm assumed}) * T
\quad\Longrightarrow\quad \text{spurious } \hat Q,\hat U.
$$

Because $T \gg P$ (temperature fluctuations are $\sim 100\times$ larger than
polarization), even a sub-percent beam error acts on a signal $\sim 100\times$
bigger than the target, so the fake polarization can rival the real $E/B$.
The leakage power spectrum inherits the shape of $C_\ell^{TT}$, so it
contaminates a **broad band** of $\ell$, not a single scale.

---

## 3. T→P leakage from sidelobes

### 3.1 Real vs assumed beam

A real beam has low-amplitude **sidelobes** far from the pointing centre
that the assumed Gaussian model does not describe. Write

$$
B_{\mathrm{real}} = B_{\mathrm{assumed}} + \Delta B,
$$

where $\Delta B$ is the unmodelled residual beam (main-beam mismatch + rings).
Convolution is linear, so the detector timestream (ignoring true $Q,U$ for a
moment) is

$$
d = B_{\mathrm{real}} * T = B_{\mathrm{assumed}} * T + \Delta B * T.
$$

The map-maker / data-reduction code only knows $B_{\mathrm{assumed}}$. It
attributes the first term to temperature and has **no free parameter left**
for the second term. That residual

$$
\delta T \equiv \Delta B * T = (B_{\mathrm{real}} - B_{\mathrm{assumed}}) * T
$$

is temperature-shaped, but the temperature solution has already been fixed
by $B_{\mathrm{assumed}}*T$. In the linear Stokes fit

$$
d(\psi) = T + Q\cos 2\psi + U\sin 2\psi,
$$

any leftover angle-independent (or poorly modelled) piece of $\delta T$ is
absorbed into the $Q/U$ solution. That is **T→P leakage**.

In Fourier space the residual is multiplicative:

$$
\delta T(\boldsymbol{\ell}) = \Delta B(\boldsymbol{\ell})\, T(\boldsymbol{\ell}),
$$

so its power spectrum is

$$
C_\ell^{\mathrm{leak}}
=
\bigl|\Delta B(\ell)\bigr|^2\, C_\ell^{TT}
\quad\text{(azimuthally averaged)}.
$$

Two immediate consequences:

1. The leakage **tracks $C_\ell^{TT}$** (broad-band contamination, not a
   single multipole).
2. Because $T_{\mathrm{rms}}\gg P_{\mathrm{rms}}$ ($\sim 100\times$), even a
   sub-percent $\Delta B$ acting on $T$ can produce a fake polarization
   competitive with true $E$ or $B$.

### 3.2 Toy model used in the notebook

The notebook builds a concrete $\Delta B$ as a **Gaussian main beam + ring
sidelobe**:

$$
B_{\mathrm{real}}(\mathbf{x})
=
B_G(\mathbf{x};\,\mathrm{fwhm})
+
f_{\mathrm{sl}}\,
\frac{B_G^{\mathrm{peak}}}{R^{\mathrm{peak}}}\,
R(\mathbf{x}),
$$

then renormalises so $\int B_{\mathrm{real}}=1$. Here $f_{\mathrm{sl}}$ is the
sidelobe peak amplitude relative to the main-beam peak, and the ring is a
radial Gaussian shell at radius $r_{\mathrm{sl}} \sim 4\times\mathrm{fwhm}$:

$$
R(\mathbf{x})
=
\exp\!\left(
-\frac{1}{2}
\frac{\bigl(|\mathbf{x}|-r_{\mathrm{sl}}\bigr)^2}{(\mathrm{fwhm}/4)^2}
\right).
$$

(This is a caricature of Airy rings / spillover, not a diffraction calculation.)
For small $f_{\mathrm{sl}}$ the difference beam is approximately linear in
$f_{\mathrm{sl}}$:

$$
\Delta B \;\approx\; f_{\mathrm{sl}}\,\Delta B_1
\qquad\Rightarrow\qquad
\delta T \;\approx\; f_{\mathrm{sl}}\,(\Delta B_1 * T).
$$

Hence the map-level leakage rms scales as

$$
\mathrm{rms}\!\left(Q_{\mathrm{leak}}^{\mathrm{(sl)}}\right)
\sim
f_{\mathrm{sl}}\, T_{\mathrm{rms}},
$$

a **straight line through the origin**, vanishing at $f_{\mathrm{sl}}=0$.
(The notebook measures this at $f_{\mathrm{sl}}=1\%,5\%,10\%$.)

### 3.3 Super-simplified projection onto $Q$

A full polarised map-maker would distribute $\delta T$ between $\hat Q$ and
$\hat U$ according to the scan angles. The notebook uses the *simplest*
proxy: put **all** of the residual into one Stokes component,

$$
Q_{\mathrm{leak}}^{\mathrm{(sl)}} = \delta T
\qquad\bigl(\text{equivalent to }\psi=0:\ 
Q_{\mathrm{leak}}=\delta T\cos 2\psi,\ 
U_{\mathrm{leak}}=\delta T\sin 2\psi\bigr).
$$

This is enough to read off the **amplitude and spectrum** of the
contamination; it is not a claim that real mapmaking always dumps everything
into $Q$.

### 3.4 Pair differencing and common-mode cancellation

For an orthogonal detector pair $(A,B)$ the half-difference is the usual
polarisation estimator. If the *far* sidelobes of $A$ and $B$ are nearly
identical (common mode),

$$
\frac{d_A-d_B}{2}
=
\frac{(B_A-B_B)*T}{2}
+
\text{(true }Q/U\text{ terms)},
$$

and $\Delta B_A \approx \Delta B_B$ for far rings $\Rightarrow$ the
sidelobe contribution **cancels**. That is why the notebook emphasises:
pair-differencing mitigates far-sidelobe T→P relative to single-detector
solving. (Near sidelobes / main-beam mismatches that differ between $A$ and
$B$ do *not* cancel — that is the differential-beam case of §5.)

---

## 4. T→P leakage from beam asymmetry (ellipticity)

### 4.1 Ideal data model vs rotated elliptical beam

For a single detector at polarisation / scan angle $\psi$, the ideal model is

$$
d(\psi) = T + Q\cos 2\psi + U\sin 2\psi.
$$

With an **asymmetric** beam the temperature that actually enters the
timestream is not the scalar $T$, but $T$ convolved with the beam **rotated
to $\psi$**:

$$
d_i = (T * B_{\psi_i}) + Q\cos 2\psi_i + U\sin 2\psi_i + \cdots.
$$

In the pure-temperature toy problem used in the notebook ($Q=U=0$ on the
sky),

$$
d_i = T * B_{\psi_i}.
$$

The elliptical Gaussian at angle $\psi$ is

$$
B_\psi(\mathbf{x})
\propto
\exp\!\left[
-\tfrac{1}{2}
\left(
\frac{X_\psi^2}{\sigma_x^2}
+
\frac{Y_\psi^2}{\sigma_y^2}
\right)
\right],
$$

with $(X_\psi,Y_\psi)$ the coordinates rotated by $\psi$, and ellipticity

$$
\varepsilon
=
\frac{\mathrm{fwhm}_x-\mathrm{fwhm}_y}{\mathrm{fwhm}_x+\mathrm{fwhm}_y}
=
\frac{\sigma_x-\sigma_y}{\sigma_x+\sigma_y}.
$$

(The notebook keeps the mean FWHM fixed: $\mathrm{fwhm}_{x,y}=\mathrm{fwhm}_{\mathrm{nom}}(1\pm\varepsilon)$.)

### 4.2 Map-maker projection

With $n$ discrete angles $\psi_i\in[0,\pi)$ and uniform coverage, the
least-squares solution $\hat m=(A^\top A)^{-1}A^\top d$ with rows
$A_i=[1,\cos 2\psi_i,\sin 2\psi_i]$ reduces to the discrete projections

$$
\begin{aligned}
\hat T &= \frac{1}{n}\sum_i d_i,\\
\hat Q &= \frac{2}{n}\sum_i d_i\cos 2\psi_i,\\
\hat U &= \frac{2}{n}\sum_i d_i\sin 2\psi_i.
\end{aligned}
$$

(The factor of $2$ is $\bigl\langle\cos^2 2\psi\bigr\rangle^{-1}=2$ over a
full period.) Substituting $d_i=T*B_{\psi_i}$ gives exactly the notebook
estimator:

$$
\hat Q
=
\frac{2}{n}\sum_i (T*B_{\psi_i})\cos 2\psi_i,
\qquad
\hat U
=
\frac{2}{n}\sum_i (T*B_{\psi_i})\sin 2\psi_i.
$$

### 4.3 Azimuthal harmonic expansion of the beam

This subsection is the idea that usually confuses people. The short version:

> Fix a sky pixel $\mathbf{x}$. As you rotate the telescope (or the beam's major
> axis) through angle $\psi$, the beam *value at that pixel*,
> $B_\psi(\mathbf{x})$, is just a **periodic function of one angle $\psi$**.
> Any periodic function of an angle can be written as a Fourier series.
> That series is the "harmonic expansion of the beam."

It is **not** a spherical-harmonic $Y_{\ell m}$ expansion of the sky. It is a
Fourier expansion **in the scan / beam-orientation angle $\psi$**, done
separately at every real-space (or Fourier) location.

#### 4.3.1 Warm-up: Fourier series of a periodic function

Any smooth function of an angle $\psi$ that is $2\pi$-periodic can be written

$$
f(\psi)
=
a_0
+
\sum_{m=1}^{\infty}
\bigl[
a_m\cos(m\psi) + b_m\sin(m\psi)
\bigr].
$$

The integer $m$ is the **azimuthal harmonic order**:

| $m$ | name | how many times it wiggles per full $2\pi$ turn |
|-----|------|-----------------------------------------------|
| $0$ | monopole | constant — no $\psi$ dependence |
| $1$ | dipole | one peak, one trough |
| $2$ | quadrupole | two peaks, two troughs |
| $4$ | hexadecapole | four of each, … |

The coefficients are ordinary Fourier integrals, e.g.

$$
a_m = \frac{1}{\pi}\int_0^{2\pi} f(\psi)\,\cos(m\psi)\,d\psi
\quad(m\ge 1),\qquad
a_0 = \frac{1}{2\pi}\int_0^{2\pi} f(\psi)\,d\psi.
$$

#### 4.3.2 Apply that idea to the *beam*, not the sky

The rotated beam $B_\psi(\mathbf{x})$ depends on:

- **where** you look on the focal plane / map: $\mathbf{x}=(x,y)$,
- **how** the beam is oriented: $\psi$.

Hold $\mathbf{x}$ fixed and only vary $\psi$. Then

$$
\psi \;\longmapsto\; B_\psi(\mathbf{x})
$$

is a number-valued function of one angle. Expand *that* function:

$$
B_\psi(\mathbf{x})
=
\sum_{m=0}^{\infty}
\Bigl[
B_m^{(c)}(\mathbf{x})\,\cos(m\psi)
+
B_m^{(s)}(\mathbf{x})\,\sin(m\psi)
\Bigr].
$$

The coefficients $B_m^{(c)}(\mathbf{x})$, $B_m^{(s)}(\mathbf{x})$ are themselves
maps (or beams) on the plane. People often package the $m$-th sector as
one complex (or two real) beam multipole map(s), and write schematically

$$
B_\psi = B_0 + B_1(\psi) + B_2(\psi) + B_3(\psi) + \cdots,
$$

where $B_m(\psi)$ collects everything proportional to $\cos(m\psi)$ and
$\sin(m\psi)$.

**Picture (one fixed pixel $\mathbf{x}$ off the beam centre):**

- Circular Gaussian: rotating does nothing $\Rightarrow$
  $B_\psi(\mathbf{x})=B_0(\mathbf{x})$ for all $\psi$. Only $m=0$.
- Slightly elliptical: as you rotate, that pixel sometimes sits more under
  the major axis, sometimes more under the minor axis $\Rightarrow$
  $B_\psi(\mathbf{x})$ oscillates twice per full $360^\circ$ turn (because an
  ellipse looks the same after $180^\circ$) $\Rightarrow$ a large $m=2$ piece.

#### 4.3.3 Why only *even* $m$ for an elliptical beam?

An elliptical (or any mirror-symmetric) beam satisfies

$$
B_{\psi+\pi}(\mathbf{x}) = B_\psi(\mathbf{x}):
$$

rotate the major axis by $180^\circ$ and the beam pattern is unchanged. In the
Fourier series that forces all **odd** harmonics to vanish
($m=1,3,5,\ldots$), because those change sign under $\psi\to\psi+\pi$.
Only

$$
m = 0,2,4,6,\ldots
$$

survive.

(If the beam had a lopsided "coma" / pointing offset that prefers one side,
you would get $m=1$. That is a different systematic.)

#### 4.3.4 Small-ellipticity expansion: $B_0 + \varepsilon B_2 + \cdots$

For a weakly elliptical Gaussian one can expand in the small parameter
$\varepsilon=(\sigma_x-\sigma_y)/(\sigma_x+\sigma_y)$. The result has the
structure

$$
B_\psi(\mathbf{x})
=
\underbrace{B_0(\mathbf{x})}_{m=0}
+
\varepsilon\,
\underbrace{B_2(\psi;\mathbf{x})}_{m=2}
+
\varepsilon^2\,
\underbrace{\bigl(B_0^{(2)} + B_4(\psi)\bigr)}_{m=0\ \mathrm{and}\ m=4}
+
\cdots
$$

Interpretation of the first two terms:

1. **$B_0$ (monopole)**  
   The circular part — average of $B_\psi$ over all orientations:

   $$
   B_0(\mathbf{x})
   =
   \frac{1}{\pi}\int_0^\pi B_\psi(\mathbf{x})\,d\psi.
   $$

   For small $\varepsilon$ this is almost the circular Gaussian with the
   mean FWHM. It does **not** depend on $\psi$.

2. **$\varepsilon B_2$ (quadrupole)**  
   The leading *anisotropy*. Its $\psi$-dependence is pure $m=2$:

   $$
   B_2(\psi;\mathbf{x})
   =
   B_2^{(c)}(\mathbf{x})\,\cos 2\psi
   +
   B_2^{(s)}(\mathbf{x})\,\sin 2\psi.
   $$

   The maps $B_2^{(c)}$, $B_2^{(s)}$ have a characteristic four-lobed
   (quadrupole) spatial pattern: positive along one axis, negative along
   the perpendicular axis. That is the mathematical shape of "stretched in
   $x$, squeezed in $y$."

Higher even multipoles ($m=4,\ldots$) appear only at $\mathcal{O}(\varepsilon^2)$
and higher for a pure ellipse; for small $\varepsilon$ the $m=2$ term
dominates the leakage.

#### 4.3.5 Concrete sketch for a weakly elliptical Gaussian

Write $\sigma_x=\sigma(1+\varepsilon)$, $\sigma_y=\sigma(1-\varepsilon)$ with
$\varepsilon\ll 1$ and mean width $\sigma$. Then

$$
\frac{X_\psi^2}{\sigma_x^2}+\frac{Y_\psi^2}{\sigma_y^2}
=
\frac{r^2}{\sigma^2}
-
2\varepsilon\,\frac{X_\psi^2 - Y_\psi^2}{\sigma^2}
+
\mathcal{O}(\varepsilon^2),
$$

where $r^2=X_\psi^2+Y_\psi^2=x^2+y^2$ is rotation-invariant. The circular
Gaussian is $B_G\propto e^{-r^2/(2\sigma^2)}$. Expanding the exponential,

$$
B_\psi
\approx
B_G
\Biggl[
1
+
\varepsilon\,\frac{X_\psi^2 - Y_\psi^2}{\sigma^2}
+
\mathcal{O}(\varepsilon^2)
\Biggr].
$$

But $X_\psi^2-Y_\psi^2$ is exactly a quadrupole in the *sky* coordinates,
and under a rotation of the axes by $\psi$ it transforms as

$$
X_\psi^2 - Y_\psi^2
=
(x^2-y^2)\cos 2\psi
+
(2xy)\sin 2\psi.
$$

So

$$
B_\psi
\approx
\underbrace{B_G}_{B_0}
+
\varepsilon\,
\underbrace{B_G\cdot
\frac{(x^2-y^2)\cos 2\psi + 2xy\sin 2\psi}{\sigma^2}
}_{B_2(\psi)}
+
\cdots
$$

This is the explicit form of "$B_0 + \varepsilon B_2$": the correction is
the circular beam times a quadrupole angular factor that spins with
$2\psi$.

#### 4.3.6 Why $m=2$ is dangerous for polarisation

Recall the map-maker templates:

$$
\hat Q \propto \sum_i d_i \cos 2\psi_i,
\qquad
\hat U \propto \sum_i d_i \sin 2\psi_i.
$$

Those templates are themselves pure **$m=2$** harmonics of $\psi$. Fourier
modes are orthogonal, so:

| beam multipole in $d_i=T*B_{\psi_i}$ | projects onto $\hat Q,\hat U$? |
|-------------------------------------|-------------------------------|
| $m=0$ (monopole / circular) | **no** — averages into $\hat T$ only |
| $m=1$ (dipole) | no (orthogonal to $\cos 2\psi$) |
| **$m=2$ (quadrupole / ellipticity)** | **yes — fully** |
| $m=4$ | no (orthogonal to $\cos 2\psi$) |
| $\ldots$ | only $m=2$ couples at linear order |

So ellipticity is special **because its $\psi$-dependence is the same
$2\psi$ dependence the polarisation map-maker is designed to extract**. The
pipeline cannot tell "true $Q$ on the sky" from "temperature modulated by a
spinning elliptical beam."

A one-line slogan:

$$
\underbrace{\text{beam ellipticity}}_{m=2\ \text{in }\psi}
\;\longleftrightarrow\;
\underbrace{\text{Stokes }Q/U}_{m=2\ \text{in }\psi}.
$$

#### 4.3.7 How this connects to the rest of §4

Once you accept

$$
B_\psi = B_0 + \varepsilon\, B_2(\psi) + \cdots,
$$

the rest is linear algebra:

1. $d_i = T * B_{\psi_i} = T*B_0 + \varepsilon\,(T*B_2(\psi_i))+\cdots$
2. Project onto $\cos 2\psi$: the $B_0$ piece dies; the $B_2$ piece survives
   (§4.4).
3. Continuum average of that projection $\Rightarrow$ one fixed kernel
   $K_Q$ with $\hat Q=\varepsilon\,(T*K_Q)$ (§4.5).

You do **not** need the full infinite series in practice: for small
$\varepsilon$, keeping $B_0$ and $\varepsilon B_2$ is enough to predict
the scaling $\mathrm{rms}(\hat Q)\propto\varepsilon$ and the fact that it
does not average down with more scan angles.

#### 4.3.8 Common confusions

- **"Is this $\ell,m$ of the CMB?"**  
  No. Here $m$ labels Fourier modes in the *instrument orientation*
  $\psi$, not sky multipoles $(\ell,m)$.

- **"Is $B_2$ a map of the sky?"**  
  No. $B_2(\psi;\mathbf{x})$ is a property of the *beam* as a function of
  focal-plane / map coordinate $\mathbf{x}$ and orientation $\psi$. The sky
  only enters when you convolve: $T*B_2$.

- **"Why expand at all? The notebook just sums over angles."**  
  The notebook's discrete sum
  $\hat Q=\frac{2}{n}\sum_i (T*B_{\psi_i})\cos 2\psi_i$
  *is* the full calculation. The harmonic expansion is the *explanation*
  of why that sum is nonzero for ellipses, linear in $\varepsilon$, and
  independent of $n$ in the large-$n$ limit.

### 4.4 Monopole projects out; quadrupole survives

Plug the expansion into $\hat Q$:

$$
\hat Q
=
\frac{2}{n}\sum_i \bigl(T*B_0\bigr)\cos 2\psi_i
+
\varepsilon\,\frac{2}{n}\sum_i \bigl(T*B_2(\psi_i)\bigr)\cos 2\psi_i
+
\cdots
$$

Over a uniform set of angles on $[0,\pi)$,

$$
\sum_i \cos 2\psi_i = 0
\quad\Rightarrow\quad
\text{the monopole term vanishes}.
$$

(Physically: a circular beam gives the same $d_i$ at every $\psi$, so it is
absorbed entirely into $\hat T$ and never looks polarised.) Only the
quadrupole (and higher even multipoles) survive:

$$
\hat Q
=
\varepsilon\,\frac{2}{n}\sum_i \bigl(T*B_2(\psi_i)\bigr)\cos 2\psi_i
+
\mathcal{O}(\varepsilon^2).
$$

### 4.5 Continuum limit: a single fixed kernel

As $n\to\infty$ the discrete sum becomes an integral. Because convolution
is linear and $T$ is independent of $\psi$,

$$
\hat Q
\;\longrightarrow\;
\varepsilon\,(T * K_Q),
$$

with the **quadrupole moment of the beam**

$$
K_Q(\mathbf{x})
=
\frac{2}{\pi}
\int_0^\pi
B_2(\psi;\mathbf{x})\,\cos 2\psi\,d\psi.
$$

Similarly

$$
\hat U
\;\longrightarrow\;
\varepsilon\,(T * K_U),
\qquad
K_U(\mathbf{x})
=
\frac{2}{\pi}
\int_0^\pi
B_2(\psi;\mathbf{x})\,\sin 2\psi\,d\psi.
$$

$K_Q$ and $K_U$ are **fixed real-space kernels** (independent of the sky and
of $n$). In Fourier space

$$
\hat Q(\boldsymbol{\ell})
=
\varepsilon\, K_Q(\boldsymbol{\ell})\, T(\boldsymbol{\ell}),
$$

so the leakage spectrum is

$$
C_\ell^{\mathrm{leak}}
\propto
\varepsilon^2\,
|K_Q(\ell)|^2\,
C_\ell^{TT}.
$$

Again the contamination is **broad-band** and tracks $C_\ell^{TT}$, filtered
by the smooth transfer $|K_Q(\ell)|^2$ of the beam quadrupole.

### 4.6 Consequence 1 — linear in $\varepsilon$ (Ex.2)

To leading order

$$
\hat Q \propto \varepsilon
\quad\Rightarrow\quad
\mathrm{rms}(\hat Q) \propto \varepsilon,
$$

a straight line through the origin, exactly zero for a circular beam
($\varepsilon=0$). This is what the notebook's Ex.2 measures when scanning
$\varepsilon\in\{0,0.05,0.10,0.15,0.20\}$.

### 4.7 Consequence 2 — does **not** average down with $n$ (Ex.1)

Write the finite-$n$ estimator as a Riemann sum for the continuum kernel:

$$
\hat Q_n
=
\varepsilon\,(T*K_Q)
+
\mathcal{O}\!\left(\frac{1}{n}\right)
\quad\text{(discretisation error)}.
$$

As $n$ increases, $\hat Q_n$ converges to the **fixed, non-zero** map
$\varepsilon(T*K_Q)$ — it does **not** converge to zero. Adding scan angles
only shrinks the $\mathcal{O}(1/n)$ sampling error, not the leakage itself.

Why? The leakage is **coherent** with the map-maker's polarisation
templates $\cos 2\psi$ and $\sin 2\psi$: it is the projection of a real
physical $\psi$-dependence of the beam onto those same templates. Random
detector noise is uncorrelated with $\cos 2\psi$ and would fall as
$1/\sqrt{n}$; beam asymmetry is *correlated*, so the projection plateaus.

This is Ex.1 of the asymmetry cell: plot $\mathrm{rms}(\hat Q)$ vs
$n_{\mathrm{angles}}$ and compare to a $1/\sqrt{n}$ reference — the leakage
stays flat while noise would drop.

**Key physics:** uniform angle coverage improves conditioning of the
$(T,Q,U)$ fit but **does not kill** asymmetry-induced T→P. One needs
**deprojection** (fit out the beam-quadrupole templates) or
**beam-matched mapmaking** (build $B_\psi$ into the pointing matrix $A$).

### 4.8 Contrast with sidelobes

| | Sidelobes (§3) | Asymmetry (§4) |
|---|---|---|
| Residual | $\Delta B * T$ (fixed beam mismatch) | $(B_\psi - B_0)*T$ (angle-dependent) |
| How it enters $Q$ | dumped by imperfect $T$ model / scan | *projected* by $\cos 2\psi$ on purpose |
| Scales as | $f_{\mathrm{sl}}$ (linear) | $\varepsilon$ (linear) |
| vs more angles | N/A in the simplified proxy | **does not** $\to 0$ |
| Pair-diff helps? | yes, if far sidelobes are common-mode | no (each detector still has $B_\psi$) |
| Spectrum | $\|\Delta B\|^2 C_\ell^{TT}$ | $\varepsilon^2\|K_Q\|^2 C_\ell^{TT}$ |

---

## 5. T→P leakage from differential beam shapes (pair differencing)

This section explains the notebook cell
`T->P leakage: Differential beam shapes (pair differencing)`.

### 5.1 Why pair differencing?

A polarisation-sensitive pixel is usually a **pair** of detectors $A$ and $B$
with orthogonal polarisation axes (e.g. $A$ sensitive to $+Q$, $B$ to $-Q$).
The ideal timestreams (no beam mismatch, pure Stokes sky) look like

$$
\begin{aligned}
d_A &= T + Q,\\
d_B &= T - Q,
\end{aligned}
$$

so the **half-difference** isolates polarisation and kills temperature:

$$
\frac{d_A - d_B}{2} = Q.
$$

(The half-sum $(d_A+d_B)/2$ recovers $T$.) This is why pair differencing is
the standard first step in many CMB polarimeters: common-mode $T$ (and many
common systematics) cancel.

**But** that cancellation assumes $A$ and $B$ see the *same* temperature
through the *same* beam. If their beams differ, temperature no longer cancels
perfectly — that residual is **differential-beam T→P leakage**.

### 5.2 Mathematics: the differential beam

Let detector $A$ have beam $B_A$ and detector $B$ have beam $B_B$. Ignoring
true polarisation for a moment (pure $T$ sky, as in the notebook),

$$
d_A = T * B_A,
\qquad
d_B = T * B_B.
$$

The half-difference is

$$
\frac{d_A - d_B}{2}
=
\frac{T*B_A - T*B_B}{2}
=
T * \underbrace{\frac{B_A - B_B}{2}}_{\delta B}.
$$

Define the **differential beam**

$$
\delta B \equiv \frac{B_A - B_B}{2}.
$$

Then, in the notebook's super-simplified proxy (all residual dumped into $Q$),

$$
Q_{\mathrm{leak}}^{\mathrm{(diff)}}
=
T * \delta B.
$$

This is exactly the code:

```python
dbeam = (B_A - B_B) / 2.0
P_leak = convolve_map_with_beam(cmb_T, dbeam)   # = T * δB
```

In Fourier space,

$$
Q_{\mathrm{leak}}(\boldsymbol{\ell})
=
\delta B(\boldsymbol{\ell})\, T(\boldsymbol{\ell}),
$$

so the leakage power spectrum is

$$
C_\ell^{\mathrm{leak}}
=
\bigl|\delta B(\ell)\bigr|^2\, C_\ell^{TT}.
$$

Same story as sidelobes: contamination **tracks $C_\ell^{TT}$** (broad-band),
filtered by the transfer function of $\delta B$.

### 5.3 Parameterisation used in the notebook

Both beams are circular Gaussians; only the FWHM differs. With a fractional
mismatch $\delta$ (called `delta_fwhm_fraction` in the code),

$$
\begin{aligned}
\mathrm{fwhm}_A &= \mathrm{fwhm}_{\mathrm{nom}}\,(1+\delta),\\
\mathrm{fwhm}_B &= \mathrm{fwhm}_{\mathrm{nom}}\,(1-\delta).
\end{aligned}
$$

Properties of this choice:

- Mean FWHM stays $\mathrm{fwhm}_{\mathrm{nom}}$:
  $(\mathrm{fwhm}_A+\mathrm{fwhm}_B)/2 = \mathrm{fwhm}_{\mathrm{nom}}$.
- Half-difference of widths is
  $(\mathrm{fwhm}_A-\mathrm{fwhm}_B)/2 = \delta\,\mathrm{fwhm}_{\mathrm{nom}}$.
- $\delta=0$ $\Rightarrow$ $B_A=B_B$ $\Rightarrow$ $\delta B=0$ $\Rightarrow$
  no leakage.
- The notebook scans $\delta \in \{0.01, 0.05, 0.10\}$ (1%, 5%, 10%).

```python
def diff_beam(delta):
    fwhm_A = fwhm_nom * (1 + delta)
    fwhm_B = fwhm_nom * (1 - delta)
    B_A = make_2d_gaussian_beam(N, pix_size, fwhm_A)
    B_B = make_2d_gaussian_beam(N, pix_size, fwhm_B)
    return (B_A - B_B) / 2.0, B_A, B_B
```

### 5.4 Why the leakage is linear in $\delta$ (Ex.2)

Write the Gaussian beam as a function of width, $B(\sigma)$. For small
mismatch $\delta$,

$$
\begin{aligned}
B_A &= B\!\bigl(\sigma(1+\delta)\bigr)
     \approx B(\sigma) + \delta\,\sigma\,\partial_\sigma B,\\
B_B &= B\!\bigl(\sigma(1-\delta)\bigr)
     \approx B(\sigma) - \delta\,\sigma\,\partial_\sigma B.
\end{aligned}
$$

Subtract and divide by 2:

$$
\delta B
=
\frac{B_A-B_B}{2}
\approx
\delta\,\sigma\,\partial_\sigma B
\equiv
\delta\, \delta B_1.
$$

Therefore

$$
Q_{\mathrm{leak}}
\approx
\delta\, (T * \delta B_1)
\quad\Rightarrow\quad
\mathrm{rms}\!\left(Q_{\mathrm{leak}}\right)
\propto
\delta.
$$

A straight line through the origin — same pattern as $f_{\mathrm{sl}}$
(sidelobes) and $\varepsilon$ (ellipticity). The notebook's Ex.2 plots
exactly this: `rms_grid` vs `delta_grid`.

### 5.5 Shape of the differential beam (Ex.4)

A narrower Gaussian is **taller and thinner**; a wider one is **shorter and
fatter** (both normalised to unit integral). So at the centre

$$
B_{\mathrm{narrow}}(0) > B_{\mathrm{wide}}(0),
$$

while in the wings the wide beam wins. With
$B_A$ wider ($\delta>0$) and $B_B$ narrower,

$$
\delta B(0) = \frac{B_A(0)-B_B(0)}{2} < 0
$$

(negative at centre), then $\delta B$ **crosses zero** and becomes positive
in the wings. That sign-changing radial profile is what the notebook plots
in the middle bottom panel:

```text
δB / nominal peak
  ^
  |      ____          ← positive wings (wide beam)
  |     /    \
--+----/------\----→ r
  |   /        \
  |__/          \__    ← negative core (narrow beam taller)
```

Intuition: pair differencing is measuring "temperature seen by $A$ minus
temperature seen by $B$." Because $A$ blurs more, it sees a smoother $T$
than $B$; the difference is a high-pass-filtered version of $T$ — i.e.
$T$ convolved with a zero-mean kernel $\delta B$.

### 5.6 What each subplot of the notebook is doing

| Panel | Exercise | Math / meaning |
|-------|----------|----------------|
| Top row: 3 leakage maps | Ex.1 | $Q_{\mathrm{leak}}=T*\delta B$ for $\delta=1\%,5\%,10\%$; maps look like filtered $T$; rms grows with $\delta$ |
| Bottom-left: $C_\ell$ | Ex.3 | $C_\ell^{\mathrm{leak}}=\|\delta B\|^2 C_\ell^{TT}$ vs true $Q$ (EE) and BB ($r=0.01$) |
| Bottom-middle: $\delta B(r)$ | Ex.4 | radial cut of $\delta B$; sign flip core→wings |
| Bottom-right: rms vs $\delta$ | Ex.2 | $\mathrm{rms}\propto\delta$ (linear) |

### 5.7 Contrast with the other two T→P channels

| | Sidelobes (§3) | Asymmetry (§4) | **Differential FWHM (§5)** |
|---|---|---|---|
| Setup | one detector, unmodelled rings | one detector, elliptical $B_\psi$ | **pair** $A,B$ with different widths |
| Residual | $(B_{\mathrm{real}}-B_{\mathrm{assumed}})*T$ | projected $T*B_\psi$ | $T*\delta B$, $\delta B=(B_A-B_B)/2$ |
| Small parameter | $f_{\mathrm{sl}}$ | $\varepsilon$ | $\delta$ (FWHM fraction) |
| Pair-diff helps? | yes (common far rings cancel) | no | **this is the residual left after pair-diff** |
| rms scaling | $\propto f_{\mathrm{sl}}$ | $\propto\varepsilon$ | $\propto\delta$ |

Slogan for this section:

> Pair differencing cancels *common* temperature.  
> Whatever is *different* between the two beams — $\delta B$ — remains as fake $Q$.

### 5.8 Practical takeaway

Matching the two detectors of a pair (same FWHM, same ellipticity, same
pointing) is as important as knowing the absolute beam. A 1% FWHM mismatch
is not a 1% error on $T$; because it multiplies the *large* temperature field
and is read as polarisation, it can still dominate a small $B$-mode target.
That is why the notebook compares leakage spectra to BB at $r=0.01$.

---

## 5A. Case study: pure-$T$ sky, two general detectors

**Setup.** The sky has **only temperature**: $Q_{\mathrm{sky}}=U_{\mathrm{sky}}=0$.
Two detectors $A$ and $B$ have

- beams $B_A(\mathbf{x})$, $B_B(\mathbf{x})$ (arbitrary — not necessarily
  equal, not necessarily Gaussian),
- polariser angles $\psi_A$, $\psi_B$ (arbitrary — **not** necessarily
  orthogonal),
- optional gains $g_A$, $g_B$ (fold into the beams if you like:
  $B_i \to g_i B_i$).

**Question.** Where does T→P leakage appear, and what is it?

### 5A.1 What the detectors actually measure

A polarising detector responds to

$$
d = B * \bigl(T + Q\cos 2\psi + U\sin 2\psi\bigr)
$$

(flat-sky, single beam for $T$ and pol). On a **pure-$T$ sky** the $Q$ and
$U$ terms vanish, so the polariser angle **drops out of the true data**:

$$
\boxed{
d_A = B_A * T,
\qquad
d_B = B_B * T.
}
$$

Physical reading: unpolarised light through a polariser is just
(half) intensity. Orientation does not matter when $Q=U=0$. All that
survives is *how each detector spatially filters $T$*.

### 5A.2 What the map-maker *thinks* the data are

The pipeline assumes a common beam (or works with already
beam-smoothed maps) and the linear Stokes model

$$
\begin{pmatrix} d_A \\ d_B \end{pmatrix}
=
\underbrace{
\begin{pmatrix}
1 & \cos 2\psi_A & \sin 2\psi_A \\
1 & \cos 2\psi_B & \sin 2\psi_B
\end{pmatrix}
}_{A\ \text{(pointing / response matrix)}}
\begin{pmatrix} T \\ Q \\ U \end{pmatrix}.
$$

With only two measurements and three Stokes parameters the system is
underdetermined. In practice one either

1. forms **pair sum / difference** (classic), or
2. solves a reduced system (e.g. estimate $(T,Q)$ in a frame where $U=0$),
   or
3. accumulates many hits with many angles and solves
   $\hat m = (A^\top A)^{-1}A^\top d$.

The leakage structure is clearest from (1)–(2); (3) is the same idea with
more rows in $A$.

### 5A.3 Decomposition: common mode vs differential mode

Rewrite the true data as

$$
\begin{aligned}
d_A &= B_\star * T + \delta B * T,\\
d_B &= B_\star * T - \delta B * T,
\end{aligned}
$$

where

$$
B_\star \equiv \frac{B_A + B_B}{2}
\quad\text{(mean / common beam)},
\qquad
\delta B \equiv \frac{B_A - B_B}{2}
\quad\text{(differential beam)}.
$$

Equivalently

$$
\boxed{
\frac{d_A + d_B}{2} = B_\star * T,
\qquad
\frac{d_A - d_B}{2} = \delta B * T.
}
$$

- **Half-sum** = temperature through the *average* beam. This is not
  “polarisation leakage”; it is just $T$ with a slightly wrong beam if the
  pipeline assumed something else.
- **Half-difference** = pure temperature residual filtered by $\delta B$.
  This is the seed of **all** T→P leakage for a two-detector observation
  of a pure-$T$ sky.

**Key fact:** for pure $T$, the half-difference depends only on
$B_A-B_B$. It does **not** depend on $\psi_A,\psi_B$. The angles only
decide *how the pipeline labels* that residual ($Q$ vs $U$ vs a mix).

### 5A.4 Where the leakage sits: general angles

The map-maker interprets $(d_A,d_B)$ through $A$. Two useful projectors:

**Intensity-like combination** (weights $w$ with $w\cdot(1,1)$ style):

$$
\hat T \;\sim\; \text{mostly }\frac{d_A+d_B}{2} = B_\star * T.
$$

**Polarisation-like combination** (weights orthogonal to the intensity
column of $A$): anything proportional to the half-difference is read as
polarisation.

Solve the $2\times 2$ system for $(T,P)$ in the polarisation direction
spanned by the two detectors. Write

$$
\mathbf{c}_i = (\cos 2\psi_i,\ \sin 2\psi_i),
\qquad i=A,B
$$

(the polarisation response vectors in the $Q$–$U$ plane). Ideal model:

$$
d_i = T + \mathbf{c}_i\cdot \mathbf{P},
\qquad \mathbf{P}=(Q,U).
$$

On pure $T$ with unequal beams the data are $d_i = B_i * T$, so the
least-squares estimate of $\mathbf{P}$ is whatever best fits the *difference*
between $d_A$ and $d_B$. One obtains the structure

$$
\boxed{
\hat{\mathbf{P}}_{\mathrm{leak}}
=
\bigl(\delta B * T\bigr)
\;\times\;
\mathbf{n}(\psi_A,\psi_B),
}
$$

where $\mathbf{n}$ is a unit-ish direction in the $(Q,U)$ plane fixed by
the geometry of $\psi_A,\psi_B$ (normalisation depends on how you invert
$A$). In components:

$$
\begin{aligned}
\hat Q_{\mathrm{leak}}
&=
\alpha_Q(\psi_A,\psi_B)\,(\delta B * T),\\
\hat U_{\mathrm{leak}}
&=
\alpha_U(\psi_A,\psi_B)\,(\delta B * T).
\end{aligned}
$$

So:

| piece | where it goes | formula |
|-------|---------------|---------|
| mean beam $B_\star$ | $\hat T$ | $B_\star * T$ |
| differential beam $\delta B$ | $\hat Q$ and/or $\hat U$ | $\delta B * T$, mixed by angles |
| true sky $Q,U$ | — | zero by assumption |

**Leakage lives entirely in the polarisation maps**, and its *sky pattern*
is always $\delta B * T$. The angles only rotate that pattern between
$\hat Q$ and $\hat U$.

### 5A.5 Special cases

#### (i) Identical beams, any angles — no T→P leakage

$B_A=B_B=B$ $\Rightarrow$ $\delta B=0$ $\Rightarrow$
$d_A=d_B=B*T$. The two detectors agree. Any polarisation estimator built
from $d_A-d_B$ gives **exactly zero**. Angles do not matter.

$$
\hat Q_{\mathrm{leak}}=\hat U_{\mathrm{leak}}=0.
$$

#### (ii) Orthogonal pair $\psi_A=0$, $\psi_B=\pi/2$ (notebook case)

Response: $d_A^{\mathrm{model}}=T+Q$, $d_B^{\mathrm{model}}=T-Q$.

Estimators:

$$
\hat T = \frac{d_A+d_B}{2} = B_\star * T,
\qquad
\hat Q = \frac{d_A-d_B}{2} = \delta B * T,
\qquad
\hat U = 0
\ \text{(not measured by this pair)}.
$$

So **all** differential-beam leakage sits in $\hat Q$:

$$
Q_{\mathrm{leak}} = T * \delta B,
\qquad
U_{\mathrm{leak}} = 0.
$$

This is exactly the notebook cell
`Q_leak = T (conv) delta_B` with
$\delta B=(B_A-B_B)/2$.

#### (iii) General non-orthogonal pair

Example: $\psi_A=0$, $\psi_B=\pi/4$ (not $90^\circ$). Ideal model

$$
d_A = T + Q,
\qquad
d_B = T + \tfrac{\sqrt{2}}{2}(Q+U).
$$

The two polarisation columns of $A$ are linearly independent, so in
principle both $Q$ and $U$ are constrained (together with $T$, still
underdetermined with 2 data — you need a prior or a third measurement).
Schematically, after projecting out $T$, the residual $d_A-d_B$ is
attributed to a linear mix:

$$
\hat Q_{\mathrm{leak}} = \alpha\,(\delta B * T),
\qquad
\hat U_{\mathrm{leak}} = \beta\,(\delta B * T),
$$

with $(\alpha,\beta)$ set by $\psi_A,\psi_B$. Same sky template
$\delta B*T$, split across $Q$ and $U$.

#### (iv) Same FWHM, but different shapes (ellipticity, pointing, …)

$\delta B$ need not come from FWHM alone. Any mismatch counts:

$$
\delta B = \tfrac12\bigl(B_A - B_B\bigr)
\quad\text{with e.g.}\quad
B_A = B_{\varepsilon_A},\;
B_B = B_{\varepsilon_B},
$$

or a relative pointing offset $\mathbf{p}$:
$B_B(\mathbf{x})=B_A(\mathbf{x}-\mathbf{p})$, etc. The leakage formula is
unchanged: $Q/U$ get $\delta B * T$.

#### (v) Same beams but wrong assumed angles

If $B_A=B_B$ and sky is pure $T$, true $d_A=d_B$, so $\hat P=0$ even if
the pipeline uses wrong $\psi$. Angle errors matter when **true
polarisation** is present (mixing $Q\leftrightarrow U$, E→B). They do
**not** by themselves create T→P from pure $T$ with matched beams.

### 5A.6 Full Stokes least-squares view (many hits)

With many samples $i=1,\ldots,N_{\mathrm{hit}}$ (many pairs or many scan
angles), stack

$$
\mathbf{d} = A\,\mathbf{m} + \text{noise},
\qquad
\mathbf{m}=(T,Q,U)^\top.
$$

The true pure-$T$ data with detector-dependent beams can be written

$$
d_i = (B_{\star,i} * T) + s_i\,(\delta B_i * T),
$$

where $s_i=\pm 1$ labels which detector (or more generally the coefficient
of the differential mode). The LS solution is

$$
\hat{\mathbf{m}} = (A^\top A)^{-1} A^\top \mathbf{d}.
$$

Split $\mathbf{d} = \mathbf{d}_{\mathrm{common}} + \mathbf{d}_{\mathrm{diff}}$:

- $\mathbf{d}_{\mathrm{common}}$ lies mostly in the column space of the
  intensity column of $A$ $\rightarrow$ absorbed into $\hat T$.
- $\mathbf{d}_{\mathrm{diff}}$ is orthogonal (or partly orthogonal) to that
  column $\rightarrow$ absorbed into $(\hat Q,\hat U)$.

So again:

$$
\boxed{
\text{T→P leakage}
=
\text{the component of }\{\delta B_i * T\}
\text{ that is not parallel to the }T\text{-column of }A.
}
$$

If the differential mode is coherent with $\cos 2\psi$ / $\sin 2\psi$
(as in the elliptical single-detector case of §4), it is *maximally*
aligned with the $Q/U$ columns and does not average down. If it is a
static pair mismatch (§5), it is a fixed map $\delta B*T$ assigned to
whatever pol direction that pair measures.

### 5A.7 One-line summary of the case study

$$
\underbrace{d_A = B_A * T,\quad d_B = B_B * T}_{\text{true pure-}T\text{ data}}
\qquad\Longrightarrow\qquad
\underbrace{\hat T \leftarrow B_\star * T}_{\text{common mode}}
,\;
\underbrace{(\hat Q,\hat U) \leftarrow \delta B * T}_{\text{differential mode = leakage}}.
$$

- **Leakage location:** polarisation maps $(\hat Q,\hat U)$, never “new”
  physics in $T$ beyond using $B_\star$ instead of some nominal beam.
- **Leakage template on the sky:** always $T$ convolved with
  $\delta B=(B_A-B_B)/2$.
- **Role of angles:** only rotate leakage between $\hat Q$ and $\hat U$;
  they do not create the residual.
- **No leakage iff** $B_A=B_B$ (matched beams), for any $\psi_A,\psi_B$.

This unifies the notebook’s differential-FWHM exercise (orthogonal pair,
$\delta B$ from FWHM) with the general pair-differencing picture.

---

## 6. E→B leakage: cross-polarization

Polarization is encoded in Stokes $Q,U$. The E and B modes are recovered by
a flat-sky rotation in Fourier space (the `QU2EB` function):

$$
fE(\ell) = fQ\cos 2\varphi + fU\sin 2\varphi,
$$
$$
fB(\ell) = -fQ\sin 2\varphi + fU\cos 2\varphi,
$$

with $\varphi = \operatorname{atan}(\ell_y/\ell_x)$ the 2D polarization angle. A perfect
detector measures $Q,U$ cleanly. A detector with **cross-polarization**
fraction $\gamma$ swaps a fraction $\gamma$ of $Q\leftrightarrow U$:

$$
Q_{\mathrm{obs}} = (1-\gamma)Q + \gamma\,U,\qquad
U_{\mathrm{obs}} = (1-\gamma)U + \gamma\,Q.
$$

Substituting into the $fB$ rotation, the recovered $B$ picks up a term
$\propto\gamma$ coming from $fE$ (since $E\gg B$):

$$
\Delta B \;\propto\; \gamma\, E \quad\Rightarrow\quad
\Delta C_\ell^{BB} \propto \gamma^2\, C_\ell^{EE}.
$$

So E→B leakage from cross-pol is **quadratic in $\gamma$** (Ex. of the
cross-pol cell: "scaling with $\gamma^2$"), and it contaminates $B$ with the
much larger $E$ signal — dangerous for primordial $B$-mode searches.

---

## 7. Putting the three T→P systematics side by side

| systematic | residual that leaks | leakage $\hat Q$ | rms scaling |
|---|---|---|---|
| sidelobe | $(B_{\mathrm{real}}-B_{\mathrm{assumed}})*T$ | $\delta T$ | $\propto f_{\mathrm{sl}}$ (linear) |
| asymmetry | $T*B_\psi$ projected on $\cos2\psi$ | $\varepsilon\,(T*K_Q)$ | $\propto\varepsilon$ (linear), **does not average down** |
| differential FWHM | $T*\delta B$, $\delta B=(B_A-B_B)/2$ | $T*\delta B$ | $\propto\delta$ (linear) |

All three are **linear in their small beam parameter** ($f_{\rm sl}$,
$\varepsilon$, $\delta$) and all produce a leakage spectrum that **tracks
$C_\ell^{TT}$** (broad-band, not a single scale), because each is ultimately
"temperature convolved with a small beam-derived kernel." The closing
exercise asks: for a fixed leakage rms of, say, $10^{-4}\,\mu$K at
cosmological scales, what $f_{\rm sl}$, $\varepsilon$, $\delta$ give the same
contamination? Because each rms is linear in its parameter, the answer is a
direct read-off from the three rms-vs-parameter curves — they are all
interchangeable up to a single proportionality constant per systematic.

---

## 8. Key takeaways

- **T→P leakage = a temperature residual with nowhere to go but $Q/U$.** It is
  always of the form $(\text{beam error}) * T$, so its spectrum follows
  $C_\ell^{TT}$ — broad-band, and amplified by $T\gg P$.
- **Sidelobes** give a residual $(B_{\mathrm{real}}-B_{\mathrm{assumed}})*T$; pair
  differencing cancels the common far-sidelobe mode.
- **Beam asymmetry** couples through the $m=2$ (quadrupole) moment of the
  beam: $\hat Q = \varepsilon\,(T*K_Q)$. It is **coherent** with the
  $\cos2\psi$ template, so it does **not** average down with more scan
  angles — the most insidious of the three, requiring deprojection.
- **Differential beams** give $\hat Q = T*\delta B$ with
  $\delta B=(B_A-B_B)/2$; linear in the FWHM mismatch $\delta$.
- **E→B leakage** (cross-pol) mixes the large $E$ into $B$:
  $\Delta C_\ell^{BB}\propto\gamma^2 C_\ell^{EE}$ — quadratic in $\gamma$,
  and the leading threat for primordial $B$-mode detection.

The unifying picture: a beam systematic is a small *kernel* derived from the
beam error; convolving the (large) temperature or E-mode field with that
kernel gives the leakage, which is why all of these scale linearly (or
quadratically, for E→B) with the size of the beam imperfection.



