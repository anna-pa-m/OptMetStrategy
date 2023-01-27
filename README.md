## OptMetStrategy: Optimal metabolic strategies for microbial growth.

This is a Matlab implementation of the microbial metabolic model in [arXiv.2210.11167](https://doi.org/10.48550/arXiv.2210.11167) addressing optimal growth strategies under several environmental conditions. The model assumes a minimal microbial model suitable for a population of *E.Coli* in a carbon-limited media, whose full description is provided by two quantities, the expenditure $\epsilon(x)$ and the specific uptake $q(x)$ governed by the internal degree of freedom $x$. Under this assumption, the growth rate is given by
#### $\qquad \qquad \mu(x,s) = \frac{\phi}{w + sq(x) + \epsilon(x)}$

where $s$ is a *stress* variable modelling the environment.

### Static environments
In `main_script.m` it is possible to select each of the four environmental conditions studied in [arXiv.2210.11167](https://doi.org/10.48550/arXiv.2210.11167) and characterized by a set of *stress* distribution $p(s)$. By adjusting the variable `type_stress` one can select one of the following stress distribution:

- `uniform` 
#### $\qquad \qquad p(s) \propto \mathbb{I}\left[ s_{min} < s < s_{max}\right]$
- `exp`
#### $\qquad \qquad p(s) \propto \exp(-s/s_{0}) \mathbb{I}\left[s_{min} < s < s_{max} \right]$
- `power_law`
#### $\qquad \qquad p(s) \propto \frac{1}{s}  \mathbb{I}\left[s_{min} < s < s_{max} \right] $
- `two_states`
#### $\qquad \qquad p(s) \propto \rho \mathcal{N}\left(s_{0}, \sigma\right) + (1 - \rho) \mathcal{N}\left(s_{1}, \sigma\right)$

The script returns the optimal strategy $p(x | s)$ where $x$ is the internal degree of freedom governing the transition between the fermentation/respiration strategies (see text for more details). The growth rate 

