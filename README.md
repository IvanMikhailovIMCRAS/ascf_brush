# Description of the computational circuit

## Potential

Effective chemical potential of the polymer monomer units in brush:
$` u(z) = V(z) - V(H) = \beta \cdot \ln(\cos(K \cdot z)/ \cos(K \cdot H))`$,

where $`\beta = 3 / \eta^2`$, and $`K = \pi/(2N)`$.

## Volume fraction

Under good solvent condition:

$` u(z) = -\ln[1-\varphi(z)]`$,

$` \varphi(z) = \exp[-u(z)]`$,

Polymer volume fraction:

$` \varphi(z)= \left[\frac{\cos(K \cdot H)}{\cos(K \cdot z)} \right]^\beta`$

## Distribution of free ends

In general

$` g(z) = \frac{a}{N \sigma} \cdot \frac{dV(z)}{dz} \int_{0}^{\tilde{z}(z)} \frac{d \varphi[V(z)-V(z')]}{dV} dz' `$

where $` V(\tilde{z}) = V(H) - V(z) `$

Under good ($`\chi=0`$) solvent conditions:

$` g(z) = - \frac{a}{N \sigma} \cdot \frac{dU(z)}{dz} \exp{[-U(z)]} \int_{0}^{\tilde{z}(z)} \exp{[-U(z')]} dz' `$

