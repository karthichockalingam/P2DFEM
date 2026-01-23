# Constants in the paper [https://doi.org/10.1016/j.est.2023.107512](https://doi.org/10.1016/j.est.2023.107512 "Persistent link using digital object identifier")

Ai, Weilong, and Yuan Liu. "Improving the convergence rate of Newman’s battery model using 2nd order finite element method." _Journal of Energy Storage_ 67 (2023): 107512. 


### Constants


| Code variable | Variable | Value   | Units                          | Description               |
| ------------- | -------- | ------- | ------------------------------ | ------------------------- |
| F             | $F$      | $96485$ | $\mathrm{C \  mol^{-1}}$       | Faraday constant          |
| R             | $R$      | $8.314$ | $\mathrm{J\ mol^{-1}\ K^{-1}}$ | Universal gas constant    |
| T_ref         | $T_0$    | $298$   | $\mathrm{K}$                   | Reference temperature     |
| r0            | $r_0$    | 1e-6    | $\mathrm{m}$                   | Reference particle radius |
| I_typ         | $I_0$    |         | $\mathrm{A}$                   | Reference total current   |


### Reference scales
$i=n,e,p,sep$, where 
- $e$ - electrolyte phase
-  $n$ - negative electrode
- $p$ - positive electrode
- $sep$ - separator

| Code variable | Variable  | Formula                                  | Unit calculation                                                                                                      | Description                                                                 |
| ------------- | --------- | ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| a0            | $a_0$     | $a_0=\dfrac{1}{r_0}$                     | $\mathrm{m^{-1}}$                                                                                                     | specific surface area scale                                                 |
| L             | $L_0$     | $L_0=L_n+L_{sep}+L_p$                    | $\mathrm{m}$                                                                                                          | total cell thickness scale                                                  |
| phi_scale     | $\phi_0$  | $\phi_0=\dfrac{R \  T_0}{F}$             | $\dfrac{\mathrm{J\ mol^{-1}\ K^{-1}}\cdot\mathrm{K}}{\mathrm{C\ mol^{-1}}}=\dfrac{\mathrm{J}}{\mathrm{C}}=\mathrm{V}$ | thermal voltage scale                                                       |
| j_scale       | $j_0$ | $j_0=\dfrac{I_0}{a_0 A L_0}$             | $\dfrac{\mathrm{A}}{(\mathrm{m^{-1}})(\mathrm{m^2})(\mathrm{m})}={\mathrm{A\ m^{-2}}}$                     | volumetric current density scale                                            |
| tn<br>tp      | $t_{i0}$  | $t_{i0}=\dfrac{F A L_0 c_{i,\max}}{I_0}$ | $\dfrac{(\mathrm{C\ mol^{-1}})(\mathrm{m^2})(\mathrm{m})(\mathrm{mol\ m^{-3}})}{\mathrm{A}}=\mathrm{s}$               | characteristic solid $(i=n,p)$ diffusion time scale for electrode particles |
| te            | $t_{e0}$  | $t_{e0}=\dfrac{F A L_0 c_{e0}}{I_0}$     | $\dfrac{(\mathrm{C\ mol^{-1}})(\mathrm{m^2})(\mathrm{m})(\mathrm{mol\ m^{-3}})}{\mathrm{A}}=\mathrm{s}$               | characteristic electrolyte diffusion time scale                             |




### Parameters of the LG M50 battery cells: Negative electrode
| Code variable                | Parameter                          |  Negative electrode | Units                                    | Description                                  |
| ---------------------------- | ---------------------------------- | ------------------: | ---------------------------------------- | -------------------------------------------- |
| negative_electrode_thickness | $L_k$                              |                85.2 | $\mu\mathrm{m}$                          | Thickness                                    |
| eps_n                        | $\epsilon_{\mathrm{e}}$            |                  25 | \%                                       | Electrolyte volume fraction                  |
|                              | $\epsilon_k$                       |                  75 | \%                                       | Active material volume fraction              |
|                              | $b_k$                              |                2.91 | —                                        | Bruggeman coefficient                        |
| cnmax                        | $c_{\text{max}}$                   |               33133 | $\mathrm{mol}/\mathrm{m}^3$              | Maximum lithium concentration                |
| cn0                          | $c_0$                              |               29866 | $\mathrm{mol}/\mathrm{m}^3$              | Initial electrode lithium concentration      |
|  sig_n                       | $\sigma_k$                         |                 215 | $\mathrm{S}/\mathrm{m}$                  | Electrode conductivity                       |
| Dn                           | $D_k$                              | $3.3\times10^{-14}$ | $\mathrm{m}^2/\mathrm{s}$                | Electrode diffusivity                        |
| rn                           | $R_k$                              |                5.86 | $\mu\mathrm{m}$                          | Particle radius                              |
| kn_dim                       | $m_k$                              | $6.48\times10^{-7}$ | $\mathrm{A\ m}^{2.5}/\mathrm{mol}^{1.5}$ | Reaction rate                                |
|                              | $c_{\mathrm{e}}^{\text{typ}}$      |                     | $\mathrm{mol}/\mathrm{m}^3$              | Typical lithium concentration in electrolyte |
|                              | $D_{\mathrm{e}}^{\text{typ}}$      |                     | $\mathrm{m}^2/\mathrm{s}$                | Typical electrolyte diffusivity              |
|                              | $c_{\mathrm{e}0}$                  |                     | $\mathrm{mol}/\mathrm{m}^3$              | Initial lithium concentration in electrolyte |
|                              | $\kappa_{\mathrm{e}}^{\text{typ}}$ |                     | $\mathrm{S}/\mathrm{m}$                  | Typical electrolyte conductivity             |
|                              | $t_{+}^0$                          |                     | —                                        | Transference number                          |






 
