## OptMetStrategy: Optimal metabolic strategies for microbial growth.

This is a Matlab implementation of the microbial metabolic model in [arXiv.2210.11167](https://doi.org/10.48550/arXiv.2210.11167) addressing optimal growth strategies under several environmental conditions.

### A brief introduction to the theoretical framework

We assume a minimal microbial model suitable for a population of *E.Coli* in a carbon-limited media, whose full description is provided by two quantities, the expenditure $\epsilon(x)$ and the specific uptake $q(x)$ governed by the internal degree of freedom $x$. Under this assumption, the growth rate is given by (see [npj Sys Biol Appl 5, 1-9](https://www.nature.com/articles/s41540-019-0093-4))
### $\qquad \qquad \mu(x,s) = \frac{\phi}{w + sq(x) + \epsilon(x)}$

where $s$ is a *stress* variable modeling the environment, $\phi$ and $w$ are the fraction of proteome devoted to constitutively expressed proteins and the proteome share to be allocated to ribosome-affiliated proteins per unit of growth rate respectively. There exists an optimal value $\hat{x}\left(s\right)$ such that the growth rate is maximized. However, being $s$ a random variable distributed according to $p(s)$, cells face a certain degree of uncertainty about the environment and, therefore, the strategy to adopt. We denote as $p(x|s)$ the stochastic rule to choose $x$ for a given $s$ which at optimality is given by 

### $\qquad \qquad p(x|s) \propto p(x)e^{-\beta \mu\left(x,s\right)}$

for any given $\beta$. See [arXiv.2210.11167](https://doi.org/10.48550/arXiv.2210.11167) for a detailed derivation.

### Documentation
#### Static environments
In `main_script.m` it is possible to select each of the four environmental conditions studied in [arXiv.2210.11167](https://doi.org/10.48550/arXiv.2210.11167) and characterized by a set of *stress* distribution $p(s)$. By adjusting the variable `type_stress` one can select one of the following stress distributions:

- `uniform` 
### $\qquad \qquad p(s) \propto \mathbb{I}\left[ s_{min} < s < s_{max}\right]$
- `exp`
### $\qquad \qquad p(s) \propto \exp(-s/s_{0}) \mathbb{I}\left[s_{min} < s < s_{max} \right]$
- `power_law`
### $\qquad \qquad p(s) \propto \frac{1}{s}  \mathbb{I}\left[s_{min} < s < s_{max} \right] $
- `two_states`
### $\qquad \qquad p(s) \propto \rho \mathcal{N}\left(s_{0}, \sigma\right) + (1 - \rho) \mathcal{N}\left(s_{1}, \sigma\right)$

The script returns the optimal strategies $p(x | s)$, the marginal distribution $p(x)$ as well as the distribution of the growth rate $p\left(\mu\right)$ for several values of the hyper-parameters $\beta \in \left[0, 3\times10^{2}\right]$ (see text for more details).

