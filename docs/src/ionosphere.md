# The Ionosphere

Note that the [TEC](https://en.wikipedia.org/wiki/Total_electron_content) of the
model ionosphere is currently hard coded to zero inside TTCal.

## Refraction

The index of refraction of a cold, collisionless, and unmagnetized plasma is

\[
    n = \sqrt{1 - \frac{\omega_p^2}{\omega^2}},
\]

where $\omega_p$ is the
[plasma frequency](https://farside.ph.utexas.edu/teaching/plasma/lectures1/node6.html)
and $\omega = 2\pi\nu$ is the angular frequency of the incident radiation.

TTCal models the ionosphere as a single uniform spherical shell. That is it does not
account for the fact that ionosphere has multiple layers or that the density can
vary within a layer. However with this simplification the refraction off the inner and
outer edges of the ionosphere can be calculated with a simple ray tracer (see the
diagram below).

| Figure 1. An illustration of the ray tracing procedure used by TTCal to model ionospheric refraction. |
|-----------------------------------------------------------|
| ![Ray tracing through the ionosphere](img/ionosphere.png) |

## Absorption

Ionospheric absorption is not currently accounted for.

