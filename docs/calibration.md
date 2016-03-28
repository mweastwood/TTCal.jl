# Calibration

Calibration is necessary to remove instrumental effects on the measured visibilities
prior to imaging. The process of calibration can generally be divided into two categories:

1. Direction independent calibration
2. Direction dependent calibration

Direction independent calibration attempts to remove instrumental effects arising from
the electronics between the antenna and the correlator. For example amplifiers and cables
can add to the amplitude and phase of the voltage respectively.
Whereas direction dependent calibration attempts to remove instrumental effects arising from
the antenna beam pattern and propagation through the ionosphere.

Both types of calibration can be viewed as trying to solve for best Jones matrices $J$
in the following equation

\[
    V_{ij}^\text{measured} = J_i V_{ij}^\text{model} J_j^* + \text{noise},
\]

where $V_{ij}$ is a 2x2 matrix containing the `xx`, `xy`, `yx`, and `yy` correlations
between antenna $i$ and antenna $j$, and $J_i$ is the Jones matrix associated with antenna $i$.
In direction independent calibration the Jones matrices are a function of frequency only.
In direction dependent calibration the Jones matrices are a function of frequency and direction.

## Iterative Least Squares

TTCal uses the method of iterative least squares independently described by
[Mitchell et al. 2008](http://adsabs.harvard.edu/abs/2008ISTSP...2..707M) and
[Salvini & Wijnholds 2014](http://adsabs.harvard.edu/abs/2014A%26A...571A..97S)
to solve for the Jones matrices.
Given a set of model visibilities we seek to minimize

\[
    \sum_{ij} \left\| V_{ij}^\text{measured} - J_i V_{ij}^\text{model} J_j^* \right\|^2
\]

If $J_i$ is assumed to be constant, the least squares solution for $J_j^*$ is

\[
    J_j^* = \frac{\sum_i \left(J_i V_{ij}^\text{model}\right)^* V_{ij}^\text{measured}}
                 {\sum_i \left\|J_i V_{ij}^\text{model}\right\|^2}.
\]

This can be computed rapidly but iterating will tend to oscillate around the
solution to the original optimization problem. These oscillations need to be damped
for rapid convergence. TTCal damps these oscillations by using Runge-Kutta steps.

## Gain Calibration

With gain calibration TTCal will solve for a set of diagonal Jones matrices
(one per antenna and frequency channel). The diagonal terms of the Jones matrix
represent the complex gain (amplitude and phase) of each antenna.

See the [cookbook](cookbook.md) for an example gain calibration on an LWA dataset.

