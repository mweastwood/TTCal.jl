# Beam Models

In TTCal, a beam model is a function that takes a frequency (measured in Hz), an azimuth,
and an elevation (both measured in radians) and returns a Jones matrix. This definition
is only workable for zenith-pointing instruments like the OVRO LWA (for which TTCal was
developed).

\[
    \text{beam model}: (\nu, {\rm az}, {\rm el}) \mapsto J
\]

Beam models are used during model visibility generation to account for attenuation,
polarization leakage, and other instrumental effects. This is accomplished by using a
[congruence transform](http://mathworld.wolfram.com/CongruenceTransformation.html).

\[
    V_\text{corrupted model} = J V_\text{pristine model} J^*,
\]

where $V$ is a 2x2 matrix containing the `xx`, `xy`, `yx`, and `yy` model visibilities.

The available beam models are listed below.

## Constant Beam

The "constant beam" simply assumes the Jones matrix is the identity at every frequency
and in every direction. It goes without saying that this beam model is not realistic
(especially away from zenith), but it can be useful.
For example if the flux of a source has already been attenuated by a beam model, you
should use the constant beam to prevent TTCal applying a second beam model.

Use `ttcal.jl ... --beam constant ...` to select the constant beam model from the
command line or use `beam = ConstantBeam()` from within a Julia script.

## Sine Beam

Also known as "ye olde beam model", this beam model has the antenna gains
proportional to $\sin({\rm el})^\alpha$, where $\alpha$ is set to $1.6$ by default.
This is a surprisingly good first-order approximation to the LWA beam, but it
does not include any frequency dependence.

\[
    J = \begin{pmatrix}
        \sin({\rm el})^{\alpha/2} && 0 \\
        0 && \sin({\rm el})^{\alpha/2} \\
    \end{pmatrix}
\]

Use `ttcal.jl ... --beam sine ...` to select the sine beam model (with $\alpha = 1.6$)
from the command line or use `beam = SineBeam(α)` from within a Julia script.

## LWA Memo 178 Beam

An even better approximation of the LWA beam comes from the LWA Memo Series.
[LWA Memo 178](http://www.faculty.ece.vt.edu/swe/lwa/memo/lwa0178a.pdf) by
Jayce Dowell does a parametric fit to the results of EM simulations of the
antennas.

The dipole gain is modeled as

\[
    p(\theta) = \left[1-\left(\frac{\theta}{\pi/2}\right)^\alpha\right]\cos^\beta\theta
              + \gamma\left(\frac{\theta}{\pi/2}\right)\cos^\delta\theta,
\]

where $\theta$ is the zenith angle, and $\alpha$, $\beta$, $\gamma$, and $\delta$ are
parameters that were fit for in the E- and H-planes of the dipole.

Use `ttcal.jl ... --beam memo178 ...` to select the LWA Memo 178 beam model from the
command line or use `beam = Memo178Beam()` from within a Julia script.

