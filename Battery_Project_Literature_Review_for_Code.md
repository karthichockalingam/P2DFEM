# Constants in the paper [https://doi.org/10.1016/j.est.2023.107512](https://doi.org/10.1016/j.est.2023.107512 "Persistent link using digital object identifier")

Ai, Weilong, and Yuan Liu. "Improving the convergence rate of Newman’s battery model using 2nd order finite element method." _Journal of Energy Storage_ 67 (2023): 107512. 

## Units
| Unit (symbol)    | Most common unit relation                          | Name       | Description                                                             |
| ---------------- | -------------------------------------------------- | ---------- | ----------------------------------------------------------------------- |
| $\mathrm{S}$     | $\mathrm{S}=\mathrm{A\ V^{-1}}$                    | Siemens    | Electrical **conductance** (inverse of resistance).                     |
| $\mathrm{V}$     | $\mathrm{V}=\mathrm{J\ C^{-1}}=\mathrm{W\ A^{-1}}$ | Volt       | Electric **potential difference** (“voltage”).                          |
| $\mathrm{A}$     | $\mathrm{A}=\mathrm{C\ s^{-1}}$                    | Ampere     | Electric **current** (flow of charge).                                  |
| $\mathrm{m}$     | —                                                  | meter      | **Length**.                                                             |
| $\mathrm{s}$     | —                                                  | second     | **Time**.                                                               |
| $\mathrm{mol}$   | —                                                  | mole       | Amount of substance (number of entities scaled by Avogadro’s constant). |
| $\mathrm{K}$     | —                                                  | Kelvin     | **Thermodynamic temperature**.                                          |
| $\mathrm{J}$     | $\mathrm{J}=\mathrm{N\ m}$                         | Joule      | **Energy** (work/heat).                                                 |
| $\mathrm{C}$     | $\mathrm{C}=\mathrm{A\ s}$                         | Coulomb    | **Electric charge**.                                                    |
| $\mathrm{N}$     | $\mathrm{N}=\mathrm{kg\ m\ s^{-2}}$                | Newton     | **Force**.                                                              |
| $\mathrm{W}$     | $\mathrm{W}=\mathrm{J\ s^{-1}}$                    | Watt       | **Power** (rate of doing work / energy per time).                       |
| $\mu \mathrm{m}$ | $\mathrm{\mu m}=10^{-6}\ \mathrm{m}$               | micrometre | Length scale commonly used for **particle radii** in battery models.    |
|                  |                                                    |            |                                                                         |

## scaled constants

| Variable | Code variable | Value   | Units                          | Description               |
| -------- | ------------- | ------- | ------------------------------ | ------------------------- |
| $F$      | F             | $96485$ | $\mathrm{C \  mol^{-1}}$       | Faraday constant          |
| $R$      | R             | $8.314$ | $\mathrm{J\ mol^{-1}\ K^{-1}}$ | Universal gas constant    |
| $T_0$    | T_ref         | $298$   | $\mathrm{K}$                   | Reference temperature     |
| $r_0$    | r0            | 1e-6    | $\mathrm{m}$                   | Reference particle radius |
| $I_0$    | I_typ         |         | $\mathrm{A}$                   | Reference total current   |


###  Parameters of the LG M50 battery cells

| Parameter                          |  Negative electrode |            Separator |  Positive electrode | Units                                    | Description                                  |
| ---------------------------------- | ------------------: | -------------------: | ------------------: | ---------------------------------------- | -------------------------------------------- |
| $L_k$                              |                85.2 |                   12 |                75.6 | $\mu\mathrm{m}$                          | Thickness                                    |
| $\epsilon_{\mathrm{e}}$            |                  25 |                   47 |                33.5 | \%                                       | Electrolyte volume fraction                  |
| $\epsilon_k$                       |                  75 |                    - |                66.5 | \%                                       | Active material volume fraction              |
| $b_k$                              |                2.91 |                 2.57 |                2.43 | —                                        | Bruggeman coefficient                        |
| $c_{\text{max}}$                   |               33133 |                    - |              63,104 | $\mathrm{mol}/\mathrm{m}^3$              | Maximum lithium concentration                |
| $c_0$                              |               29866 |                    - |              17,038 | $\mathrm{mol}/\mathrm{m}^3$              | Initial electrode lithium concentration      |
| $\sigma_k$                         |                 215 |                    - |                0.18 | $\mathrm{S}/\mathrm{m}$                  | Electrode conductivity                       |
| $D_k$                              | $3.3\times10^{-14}$ |                    - |   $4\times10^{-15}$ | $\mathrm{m}^2/\mathrm{s}$                | Electrode diffusivity                        |
| $R_k$                              |                5.86 |                    - |                5.22 | $\mu\mathrm{m}$                          | Particle radius                              |
| $m_k$                              | $6.48\times10^{-7}$ |                    - | $3.42\times10^{-6}$ | $\mathrm{A\ m}^{2.5}/\mathrm{mol}^{1.5}$ | Reaction rate                                |
| $c_{\mathrm{e}}^{\text{typ}}$      |                     |        $1\times10^3$ |                     | $\mathrm{mol}/\mathrm{m}^3$              | Typical lithium concentration in electrolyte |
| $D_{\mathrm{e}}^{\text{typ}}$      |                     | $5.34\times10^{-10}$ |                     | $\mathrm{m}^2/\mathrm{s}$                | Typical electrolyte diffusivity              |
| $c_{\mathrm{e}0}$                  |                     |                 1000 |                     | $\mathrm{mol}/\mathrm{m}^3$              | Initial lithium concentration in electrolyte |
| $\kappa_{\mathrm{e}}^{\text{typ}}$ |                     |                  1.1 |                     | $\mathrm{S}/\mathrm{m}$                  | Typical electrolyte conductivity             |
| $t_{+}^0$                          |                     |               0.2594 |                     | —                                        | Transference number                          |



### Negative electrode code variables

| Code variable                | Parameter                          |  Negative electrode | Units                                    | Description                                  |
| ---------------------------- | ---------------------------------- | ------------------: | ---------------------------------------- | -------------------------------------------- |
| negative_electrode_thickness | $L_k$                              |                85.2 | $\mu\mathrm{m}$                          | Thickness                                    |
| eps_n                        | $\epsilon_{\mathrm{e}}$            |                  25 | \%                                       | Electrolyte volume fraction                  |
|                              | $\epsilon_k$                       |                  75 | \%                                       | Active material volume fraction              |
|                              | $b_k$                              |                2.91 | —                                        | Bruggeman coefficient                        |
| cnmax                        | $c_{\text{max}}$                   |               33133 | $\mathrm{mol}/\mathrm{m}^3$              | Maximum lithium concentration                |
| cn0                          | $c_0$                              |               29866 | $\mathrm{mol}/\mathrm{m}^3$              | Initial electrode lithium concentration      |
|                              | $\sigma_k$                         |                 215 | $\mathrm{S}/\mathrm{m}$                  | Electrode conductivity                       |
| Dn                           | $D_k$                              | $3.3\times10^{-14}$ | $\mathrm{m}^2/\mathrm{s}$                | Electrode diffusivity                        |
| rn                           | $R_k$                              |                5.86 | $\mu\mathrm{m}$                          | Particle radius                              |
| kn_dim                       | $m_k$                              | $6.48\times10^{-7}$ | $\mathrm{A\ m}^{2.5}/\mathrm{mol}^{1.5}$ | Reaction rate                                |
|                              | $c_{\mathrm{e}}^{\text{typ}}$      |                     | $\mathrm{mol}/\mathrm{m}^3$              | Typical lithium concentration in electrolyte |
|                              | $D_{\mathrm{e}}^{\text{typ}}$      |                     | $\mathrm{m}^2/\mathrm{s}$                | Typical electrolyte diffusivity              |
|                              | $c_{\mathrm{e}0}$                  |                     | $\mathrm{mol}/\mathrm{m}^3$              | Initial lithium concentration in electrolyte |
|                              | $\kappa_{\mathrm{e}}^{\text{typ}}$ |                     | $\mathrm{S}/\mathrm{m}$                  | Typical electrolyte conductivity             |
|                              | $t_{+}^0$                          |                     | —                                        | Transference number                          |


## Parameters

$i=n,e,p,sep$, where 
- $e$ - electrolyte phase
-  $n$ - negative electrode
- $p$ - positive electrode
- $sep$ - separator

| Variable          | Units                                       | Description                                |
| ----------------- | ------------------------------------------- | ------------------------------------------ |
| $\kappa$          | $\mathrm{S \  m^{-1}}$                        | conductivity of electrolyte                |
| $\phi_i$          | $\mathrm{V}$                                | potential of phase $i$                     |
| $\sigma_i$        | $\mathrm{S \  m^{-1}}$                        | conductivity of phase $i$                  |
| $A$               | $\mathrm{m^2}$                              | cell surface area                          |
| $a_i$             | $\mathrm{m^{-1}}$                           | specific surface area of phase $i$         |
| $c_i$             | $\mathrm{mol\  m^{-3}}$                       | concentration of lithium in phase $i$      |
| $D_i$             | $\mathrm{m^2 \  s^{-1}}$                       | lithium diffusion coefficient of phase $i$ |
| $D_i^{\mathrm{eff}}$ | $\mathrm{m^2 \  s^{-1}}$                   | effective diffusion coefficient of phase $i$ |
| $I$               | $\mathrm{A}$                                | electronic current                         |
| $j$               | $\mathrm{A\ m^{-3}}$                         | reaction current density per volume        |
| $L$               | $\mathrm{m}$                                | cell thickness                             |
| $L_i$             | $\mathrm{m}$                                | thickness of domain $i$                    |
| $m_i$             | $\mathrm{A\ m^{2.5}\ mol^{-1.5}}$             | reaction rate in phase $i$                 |
| $R_i$             | $\mathrm{\mu m}$                            | particle radius of phase $i$               |
| $T$               | $\mathrm{K}$                                | absolute temperature                       |
| $t_+^{0}$         | —                                           | transference number of lithium-ion         |
| $t_+$             | —                                           | transference number of lithium-ion         |
| $U_i$             | $\mathrm{V}$                                | open circuit potential of phase $i$        |


### Reference scales

| Code variable | Variable  | Formula                                  | Unit calculation                                                                                                      | Description                                                                 |
| ------------- | --------- | ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| a0            | $a_0$     | $a_0=\dfrac{1}{r_0}$                     | $\mathrm{m^{-1}}$                                                                                                     | specific surface area scale                                                 |
| L             | $L_0$     | $L_0=L_n+L_{sep}+L_p$                    | $\mathrm{m}$                                                                                                          | total cell thickness scale                                                  |
| phi_scale     | $\phi_0$  | $\phi_0=\dfrac{R \  T_0}{F}$             | $\dfrac{\mathrm{J\ mol^{-1}\ K^{-1}}\cdot\mathrm{K}}{\mathrm{C\ mol^{-1}}}=\dfrac{\mathrm{J}}{\mathrm{C}}=\mathrm{V}$ | thermal voltage scale                                                       |
| j_scale       | $j_0$ | $j_0=\dfrac{I_0}{a_0 A L_0}$             | $\dfrac{\mathrm{A}}{(\mathrm{m^{-1}})(\mathrm{m^2})(\mathrm{m})}={\color{red}\mathrm{A\ m^{-2}}}$                     | volumetric current density scale                                            |
| tn<br>tp      | $t_{i0}$  | $t_{i0}=\dfrac{F A L_0 c_{i,\max}}{I_0}$ | $\dfrac{(\mathrm{C\ mol^{-1}})(\mathrm{m^2})(\mathrm{m})(\mathrm{mol\ m^{-3}})}{\mathrm{A}}=\mathrm{s}$               | characteristic solid $(i=n,p)$ diffusion time scale for electrode particles |
| te            | $t_{e0}$  | $t_{e0}=\dfrac{F A L_0 c_{e0}}{I_0}$     | $\dfrac{(\mathrm{C\ mol^{-1}})(\mathrm{m^2})(\mathrm{m})(\mathrm{mol\ m^{-3}})}{\mathrm{A}}=\mathrm{s}$               | characteristic electrolyte diffusion time scale                             |

## How variables are scaled


| Variable                   | Formula                                                                | Unit calculation (from formula)                                       | Description                                         |
| -------------------------- | ---------------------------------------------------------------------- | --------------------------------------------------------------------- | --------------------------------------------------- |
| $\bar{c}_e$                | $`\bar{c}_e = \dfrac{c_e}{c_{e0}}`$                                      | $\dfrac{\mathrm{mol\ m^{-3}}}{\mathrm{mol\ m^{-3}}}=1$                | scaled electrolyte concentration                    |
| $\bar{x}$                  | $\bar{x} = \dfrac{x}{L_0}$                                             | $\dfrac{\mathrm{m}}{\mathrm{m}}=1$                                    | scaled spatial coordinate                           |
| $\bar{a}_i$                | $\bar{a}_i = \dfrac{a_i}{a_0}$                                         | $\dfrac{\mathrm{m^{-1}}}{\mathrm{m^{-1}}}=1$                          | scaled specific surface area of phase $i$           |
| $\bar{D}_e^{\mathrm{eff}}$ | $`\bar{D}_e^{\mathrm{eff}} = D_e^{\mathrm{eff}}\dfrac{t_{e0}}{L_0^{2}}`$ | $\left(\mathrm{m^2\ s^{-1}}\right)\dfrac{\mathrm{s}}{\mathrm{m^2}}=1$ | scaled effective electrolyte diffusivity            |
| $\bar{t}_e$                | $`\bar{t}_e = \dfrac{t}{t_{e0}}`$                                        | $\dfrac{\mathrm{s}}{\mathrm{s}}=1$                                    | scaled time                                         |
|                            |                                                                        |                                                                       |                                                     |
| $\bar{r}$                  | $\bar{r} = \dfrac{r}{r_0}$                                             | $\dfrac{\mathrm{m}}{\mathrm{m}}=1$                                    | scaled radial coordinate                            |
| $\bar{R}_i$                | $\bar{R}_i = \dfrac{R_i}{r_0}$                                         | $\dfrac{\mathrm{\mu m}}{\mathrm{\mu m}}=1$                            | scaled particle radius of phase $i$                 |
| $\bar{c}_i$                | $`\bar{c}_i = \dfrac{c_i}{c_{i,\max}}`$                                  | $\dfrac{\mathrm{mol\ m^{-3}}}{\mathrm{mol\ m^{-3}}}=1$                | scaled lithium concentration in phase $i$           |
| $\bar{D}_i^{\mathrm{eff}}$ | $`\bar{D}_i^{\mathrm{eff}} = D_i^{\mathrm{eff}}\dfrac{t_{i0}}{r_0^{2}}`$ | $\left(\mathrm{m^2\ s^{-1}}\right)\dfrac{\mathrm{s}}{\mathrm{m^2}}=1$ | scaled effective diffusion coefficient in phase $i$ |
| $\bar{t}_i$                | $`\bar{t}_i = \dfrac{t}{t_{i0}}`$                                        | $\dfrac{\mathrm{s}}{\mathrm{s}}=1$                                    | scaled time for phase-$i$ diffusion                 |


## Potentials, current, and reaction terms


| Variable                          | Formula                                                                              | Unit calculation (from formula)                                                                                                                                                                                | Description                                 |
| --------------------------------- | ------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| $\bar{T}$                         | $\bar{T}=\dfrac{T}{T_0}$                                                             | $\dfrac{\mathrm{K}}{\mathrm{K}}=1$                                                                                                                                                                             | scaled temperature                          |
| $\bar{\phi}_e$                    | $`\bar{\phi}_e=\dfrac{\phi_e}{\phi_0}`$                                                | $\dfrac{\mathrm{V}}{\mathrm{V}}=1$                                                                                                                                                                             | scaled electrolyte potential                |
| $\bar{\phi}_i$                    | $`\bar{\phi}_i=\dfrac{\phi_i}{\phi_0}`$                                                | $\dfrac{\mathrm{V}}{\mathrm{V}}=1$                                                                                                                                                                             | scaled potential of phase $i$               |
| $\bar{\kappa}^{\mathrm{eff}}$     | $`\bar{\kappa}^{\mathrm{eff}}=\kappa^{\mathrm{eff}}\dfrac{\phi_0 A}{I_0 L_0}`$         | $\left(\dfrac{\mathrm{S}}{\mathrm{m}}\right)\dfrac{\mathrm{V}\ \mathrm{m^2}}{\mathrm{A}\ \mathrm{m}}=\left(\dfrac{\mathrm{A}}{\mathrm{V\ m}}\right)\dfrac{\mathrm{V}\ \mathrm{m^2}}{\mathrm{A}\ \mathrm{m}}=1$ | scaled effective electrolyte conductivity   |
| $\bar{\kappa}^{\mathrm{eff}}_{D}$ | $`\bar{\kappa}^{\mathrm{eff}}_{D}=2\ \bar{T}\ \bar{\kappa}^{\mathrm{eff}}(1-t_+)`$     | $1\times 1\times 1=1$                                                                                                                                                                                          | scaled “diffusional” conductivity factor    |
| $\bar{\sigma}^{\mathrm{eff}}_{i}$ | $`\bar{\sigma}^{\mathrm{eff}}_{i}=\sigma^{\mathrm{eff}}_{i}\dfrac{\phi_0 A}{I_0 L_0}`$ | $\left(\dfrac{\mathrm{S}}{\mathrm{m}}\right)\dfrac{\mathrm{V}\ \mathrm{m^2}}{\mathrm{A}\ \mathrm{m}}=1$                                                                                                        | scaled effective conductivity of phase $i$  |
| $\bar{j}$                         | $\bar{j}=\dfrac{j}{j_0}$                                                             | $\dfrac{\mathrm{A\ m^{-3}}}{\mathrm{A\ m^{-3}}}=1$                                                                                                                                                             | scaled reaction current density per volume  |
| $\bar{I}$                         | $\bar{I}=\dfrac{I}{I_0}$                                                             | $\dfrac{\mathrm{A}}{\mathrm{A}}=1$                                                                                                                                                                             | scaled applied/electronic current           |
| $\bar{m}_i$                       | $`\bar{m}_i=\dfrac{m_i \ (c_{i,\max}\ c_{e0}^{0.5})}{j_0}`$                            | $\dfrac{\mathrm{A\ m^{2.5}\ mol^{-1.5}}\left(\mathrm{mol\ m^{-3}}\right)^{3/2}}{\mathrm{A\ m^{-2}}}=1$                                                                                                         | scaled reaction-rate parameter in phase $i$ |
| $\bar{U}_i$                       | $\bar{U}_i=\dfrac{U_i}{\phi_0}$                                                      | $\dfrac{\mathrm{V}}{\mathrm{V}}=1$                                                                                                                                                                             | scaled open-circuit potential of phase $i$  |
 
