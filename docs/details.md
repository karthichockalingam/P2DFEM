## Constants in the paper [https://doi.org/10.1016/j.est.2023.107512](https://doi.org/10.1016/j.est.2023.107512 "Persistent link using digital object identifier")

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
- $n$ - negative electrode
- $p$ - positive electrode
- $sep$ - separator

| Code variable | Variable  | Formula                                  | Units                                   | Description                                                                 |
| ------------- | --------- | ---------------------------------------- | --------------------------------------- | --------------------------------------------------------------------------- |
| a0            | $a_0$     | $a_0=\dfrac{1}{r_0}$                     | $\mathrm{m^{-1}}$                       | specific surface area scale                                                 |
| L             | $L_0$     | $L_0=L_n+L_{sep}+L_p$                    | $\mathrm{m}$                            | total cell thickness scale                                                  |
| phi_scale     | $\phi_0$  | $\phi_0=\dfrac{R \  T_0}{F}$             | $\mathrm{V}$                            | thermal voltage scale                                                       |
| j_scale       | $j_0$     | $j_0=\dfrac{I_0}{a_0 A L_0}$             | $\mathrm{A\ m^{-2}}$                    | volumetric current density scale                                            |
| tn<br>tp      | $t_{i0}$  | $t_{i0}=\dfrac{F A L_0 c_{i,\max}}{I_0}$ | $\mathrm{s}$                            | characteristic solid $(i=n,p)$ diffusion time scale for electrode particles |
| te            | $t_{e0}$  | $t_{e0}=\dfrac{F A L_0 c_{e0}}{I_0}$     | $\mathrm{s}$                            | characteristic electrolyte diffusion time scale                             |


### Parameters

| Variable          | Units                                       | Description                                |
| ----------------- | ------------------------------------------- | ------------------------------------------ |
| $\kappa$          | $\mathrm{S \  m^{-1}}$                      | conductivity of electrolyte                |
| $\phi_i$          | $\mathrm{V}$                                | potential of phase $i$                     |
| $\sigma_i$        | $\mathrm{S \  m^{-1}}$                      | conductivity of phase $i$                  |
| $A$               | $\mathrm{m^2}$                              | cell surface area                          |
| $a_i$             | $\mathrm{m^{-1}}$                           | specific surface area of phase $i$         |
| $c_i$             | $\mathrm{mol\  m^{-3}}$                     | concentration of lithium in phase $i$      |
| $D_i$             | $\mathrm{m^2 \  s^{-1}}$                    | lithium diffusion coefficient of phase $i$ |
| $D_i^{\mathrm{eff}}$ | $\mathrm{m^2 \  s^{-1}}$                 | effective diffusion coefficient of phase $i$ |
| $I$               | $\mathrm{A}$                                | electronic current                         |
| $j$               | $\mathrm{A\ m^{-3}}$                        | reaction current density per volume        |
| $L$               | $\mathrm{m}$                                | cell thickness                             |
| $L_i$             | $\mathrm{m}$                                | thickness of domain $i$                    |
| $m_i$             | $\mathrm{A\ m^{2.5}\ mol^{-1.5}}$           | reaction rate in phase $i$                 |
| $R_i$             | $\mathrm{\mu m}$                            | particle radius of phase $i$               |
| $T$               | $\mathrm{K}$                                | absolute temperature                       |
| $t_+^{0}$         | —                                           | transference number of lithium-ion         |
| $t_+$             | —                                           | transference number of lithium-ion         |
| $U_i$             | $\mathrm{V}$                                | open circuit potential of phase $i$        |



### Scaled parameters


| Variable                   | Formula                                                                | Unit calculation (from formula)                                       | Description                                         |
| -------------------------- | ---------------------------------------------------------------------- | --------------------------------------------------------------------- | --------------------------------------------------- |
| $\bar{c}_e$                | $`\bar{c}_e = \dfrac{c_e}{c_{e0}}`$                                      | $\dfrac{\mathrm{mol\ m^{-3}}}{\mathrm{mol\ m^{-3}}}=1$                | scaled electrolyte concentration                    |
| $\bar{x}$                  | $\bar{x} = \dfrac{x}{L_0}$                                               | $\dfrac{\mathrm{m}}{\mathrm{m}}=1$                                    | scaled spatial coordinate                           |
| $\bar{a}_i$                | $\bar{a}_i = \dfrac{a_i}{a_0}$                                           | $\dfrac{\mathrm{m^{-1}}}{\mathrm{m^{-1}}}=1$                          | scaled specific surface area of phase $i$           |
| $\bar{D}_e^{\mathrm{eff}}$ | $`\bar{D}_e^{\mathrm{eff}} = D_e^{\mathrm{eff}}\dfrac{t_{e0}}{L_0^{2}}`$ | $\left(\mathrm{m^2\ s^{-1}}\right)\dfrac{\mathrm{s}}{\mathrm{m^2}}=1$ | scaled effective electrolyte diffusivity            |
| $\bar{t}_e$                | $`\bar{t}_e = \dfrac{t}{t_{e0}}`$                                        | $\dfrac{\mathrm{s}}{\mathrm{s}}=1$                                    | scaled time                                         |
| $\bar{r}$                  | $\bar{r} = \dfrac{r}{r_0}$                                               | $\dfrac{\mathrm{m}}{\mathrm{m}}=1$                                    | scaled radial coordinate                            |
| $\bar{R}_i$                | $\bar{R}_i = \dfrac{R_i}{r_0}$                                           | $\dfrac{\mathrm{\mu m}}{\mathrm{\mu m}}=1$                            | scaled particle radius of phase $i$                 |
| $\bar{c}_i$                | $`\bar{c}_i = \dfrac{c_i}{c_{i,\max}}`$                                  | $\dfrac{\mathrm{mol\ m^{-3}}}{\mathrm{mol\ m^{-3}}}=1$                | scaled lithium concentration in phase $i$           |
| $\bar{D}_i^{\mathrm{eff}}$ | $`\bar{D}_i^{\mathrm{eff}} = D_i^{\mathrm{eff}}\dfrac{t_{i0}}{r_0^{2}}`$ | $\left(\mathrm{m^2\ s^{-1}}\right)\dfrac{\mathrm{s}}{\mathrm{m^2}}=1$ | scaled effective diffusion coefficient in phase $i$ |
| $\bar{t}_i$                | $`\bar{t}_i = \dfrac{t}{t_{i0}}`$                                        | $\dfrac{\mathrm{s}}{\mathrm{s}}=1$                                    | scaled time for phase-$i$ diffusion                 |
| $\bar{T}$                  | $\bar{T}=\dfrac{T}{T_0}$                                                 | $\dfrac{\mathrm{K}}{\mathrm{K}}=1$                                    | scaled temperature                                  |
| $\bar{\phi}_e$             | $`\bar{\phi}_e=\dfrac{\phi_e}{\phi_0}`$                                  | $\dfrac{\mathrm{V}}{\mathrm{V}}=1$                                    | scaled electrolyte potential                        |
| $\bar{\phi}_i$             | $`\bar{\phi}_i=\dfrac{\phi_i}{\phi_0}`$                                  | $\dfrac{\mathrm{V}}{\mathrm{V}}=1$                                    | scaled potential of phase $i$               |
| $\bar{\kappa}^{\mathrm{eff}}$     | $`\bar{\kappa}^{\mathrm{eff}}=\kappa^{\mathrm{eff}}\dfrac{\phi_0 A}{I_0 L_0}`$   | $\left(\dfrac{\mathrm{S}}{\mathrm{m}}\right)\dfrac{\mathrm{V}\ \mathrm{m^2}}{\mathrm{A}\ \mathrm{m}}=\left(\dfrac{\mathrm{A}}{\mathrm{V\ m}}\right)\dfrac{\mathrm{V}\ \mathrm{m^2}}{\mathrm{A}\ \mathrm{m}}=1$  | scaled effective electrolyte conductivity   |
| $\bar{\kappa}^{\mathrm{eff}}_{D}$ | $`\bar{\kappa}^{\mathrm{eff}}_{D}=2\ \bar{T}\ \bar{\kappa}^{\mathrm{eff}}(1-t_+)`$ | $1\times 1\times 1=1$                                | scaled “diffusional” conductivity factor    |
| $\bar{\sigma}^{\mathrm{eff}}_{i}$ | $`\bar{\sigma}^{\mathrm{eff}}_{i}=\sigma^{\mathrm{eff}}_{i}\dfrac{\phi_0 A}{I_0 L_0}`$ | $\left(\dfrac{\mathrm{S}}{\mathrm{m}}\right)\dfrac{\mathrm{V}\ \mathrm{m^2}}{\mathrm{A}\ \mathrm{m}}=1$                                                                                                        | scaled effective conductivity of phase $i$  |
| $\bar{j}$                  | $\bar{j}=\dfrac{j}{j_0}$                                                  | $\dfrac{\mathrm{A\ m^{-3}}}{\mathrm{A\ m^{-3}}}=1$                   | scaled reaction current density per volume  |
| $\bar{I}$                  | $\bar{I}=\dfrac{I}{I_0}$                                                  | $\dfrac{\mathrm{A}}{\mathrm{A}}=1$                                   | scaled applied/electronic current           |
| $\bar{m}_i$                | $`\bar{m}_i=\dfrac{m_i \ (c_{i,\max}\ c_{e0}^{0.5})}{j_0}`$               | $\dfrac{\mathrm{A\ m^{2.5}\ mol^{-1.5}}\left(\mathrm{mol\ m^{-3}}\right)^{3/2}}{\mathrm{A\ m^{-2}}}=1$                                                                                                                           | scaled reaction-rate parameter in phase $i$ |
| $\bar{U}_i$                | $\bar{U}_i=\dfrac{U_i}{\phi_0}$                                           | $\dfrac{\mathrm{V}}{\mathrm{V}}=1$                                    | scaled open-circuit potential of phase $i$  |


### MFEM scaling

The problem is solved on a uniform grid $\bar{x} \in[0,1]$ with spacing $\Delta \bar{x}_{\mathrm{ref}}=1 / N X$. 


In real physical space $x$, we have different normalized lengths $L N E, L S E P, L P E$ for different regions. Then the physical element sizes become:
- NE elements: $\Delta x_{N E}=L N E / N N E$
- SEP elements: $\Delta x_{S E P}=L S E P / N S E P$
- PE elements: $\Delta x_{P E}=L P E / N P E$

To address this, the physical coordinate is represented by a piecewise linear mapping $$x=F(\bar{x}),$$ 

$$
F(\bar{x})= 
\begin{cases}
\frac{L N E \cdot N X}{N N E} \bar{x}, & \bar{x} \in\left[0,\frac{N N E}{N X}\right],\\
L N E+\frac{L S E P \cdot N X}{N S E P}\left(\bar{x}-\frac{N N E}{N X}\right), & \bar{x} \in\left[\frac{N N E}{N X},\frac{N N E+N S E P}{N X}\right],\\
L N E+L S E P+\frac{L P E \cdot N X}{N P E}\left(\bar{x}-\frac{N N E+N S E P}{N X}\right), & \bar{x} \in\left[\frac{N N E+N S E P}{N X},1\right],
\end{cases}
$$

with the piecewise constant Jacobian

$$
F'(\bar{x})=\omega_i=
\begin{cases}
\frac{L N E \cdot N X}{N N E}, &\quad i = N E,\\
\frac{L S E P \cdot N X}{N S E P}, &\quad i = S E P,\\
\frac{L P E \cdot N X}{N P E}, &\quad i = P E.
\end{cases}
$$

Let $\bar{\phi}=\phi \circ F$ and $\bar{v}=v \circ F$. We then apply the change of variables to the integrals below

$$
\int \sigma \phi_x v_x d x=\int \sigma\left(\frac{1}{\omega_i} \bar\phi_{\bar{x}}\right)\left(\frac{1}{\omega_i} \bar v_{\bar{x}}\right) \omega_i d \bar{x}=\int \frac{\sigma}{\omega_i} \ \bar\phi_{\bar{x}} \bar v_{\bar{x}} d \bar{x},
$$

$$
\int f v(x) d x=\int f \bar{v}(\bar{x}) \ \omega_i d \bar{x}=\int f \omega_i \ \bar{v}(\bar{x}) d \bar{x},
$$

where $\sigma$ and $f$ are constants.

### Parameters of the LG M50 battery cells
[https://iopscience.iop.org/article/10.1149/1945-7111/ab9050](https://iopscience.iop.org/article/10.1149/1945-7111/ab9050)

Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W. Dhammika Widanage and Emma Kendrick "Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models"

#### Negative electrode
| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
| cn0                          | $c_0$                         |               29866 | $\mathrm{mol}/\mathrm{m}^3$              | initial electrode lithium concentration      |
| cnmax                        | $c_{\text{max}}$              |               33133 | $\mathrm{mol}/\mathrm{m}^3$              | maximum lithium concentration                |
|                              |                               |                     |                                          |                                              |
| rn                           | $R_k$                         | $5.86\times10^{-6}$ | $\mathrm{m}$                             | particle radius, negative electrode          |
| negative_electrode_thickness | $L_k$                         | $85.2\times10^{-6}$ | $\mathrm{m}$                             | negative electrode thickness                 |
|                              |                               |                     |                                          |                                              |
| eps_n                        | $\epsilon_{\mathrm{e}}$       |                  25 | \%                                       | electrolyte volume fraction                  |
| eps_n_s                      | $\epsilon_k$                  |                  75 | \%                                       | active material volume fraction              |
|                              |                               |                     |                                          |                                              |
| kn_dim                       | $m_k$                         | $6.48\times10^{-7}$ | $\mathrm{A\ m}^{2.5}/\mathrm{mol}^{1.5}$ | reaction rate                                |
| Dn                           | $D_k$                         | $3.3\times10^{-14}$ | $\mathrm{m}^2/\mathrm{s}$                | electrode diffusivity                        |
| An                           | $a_k=3 \epsilon_{k}/R_k$      |                     | $\mathrm{m}^{-1}$                        | specific surface area of negative electrode  |
|                              |                               |                     |                                          |                                              |
|  sig_n                       | $\sigma_k$                    |                 215 | $\mathrm{S}/\mathrm{m}$                  | negative electrode conductivity              |
|                              |                               |                     |                                          |                                              |
| Un()                         | $U_{-}()$                     |                     | $\mathrm{V}$                             | OCV curve for the negative electrode         |

#### Positive electrode

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
| cp0                          | $c_0$                         |               17038 | $\mathrm{mol}/\mathrm{m}^3$              | initial electrode lithium concentration      |
| cpmax                        | $c_{\text{max}}$              |               63104 | $\mathrm{mol}/\mathrm{m}^3$              | maximum lithium concentration                |
|                              |                               |                     |                                          |                                              |
| rp                           | $R_k$                         | $5.22\times10^{-6}$ | $\mathrm{m}$                             | particle radius, negative electrode          |
| positive_electrode_thickness | $L_k$                         | $75.6\times10^{-6}$ | $\mathrm{m}$                             | positive electrode thickness                 |
|                              |                               |                     |                                          |                                              |
| eps_p                        | $\epsilon_{\mathrm{e}}$       |                33.5 | \%                                       | electrolyte volume fraction                  |
| eps_p_s                      | $\epsilon_k$                  |                66.5 | \%                                       | active material volume fraction              |
|                              |                               |                     |                                          |                                              |
| kp_dim                       | $m_k$                         | $3.42\times10^{-6}$ | $\mathrm{A\ m}^{2.5}/\mathrm{mol}^{1.5}$ | reaction rate                                |
| Dp                           | $D_k$                         | $4.0\times10^{-15}$ | $\mathrm{m}^2/\mathrm{s}$                | electrode diffusivity                        |
| Ap                           | $a_k=3 \epsilon_{k}/R_k$      |                     | $\mathrm{m}^{-1}$                        | specific surface area of positive electrode  |
|                              |                               |                     |                                          |                                              |
|  sig_p                       | $\sigma_k$                    |                0.18 | $\mathrm{S}/\mathrm{m}$                  | negative electrode conductivity              |
|                              |                               |                     |                                          |                                              |
| Up()                         | $U_{+}()$                     |                     | $\mathrm{V}$                             | OCV curve for the positive electrode         |

#### Separator

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
| separator_thickness          | $L_s$                         | $12\times10^{-6}$   | $\mathrm{m}$                             | separator thickness                          | 
| tplus                        | $t^{+}$                       | 0.2594              | -                                        | transference number of lithium ions          |
| eps_s                        | $\epsilon_{\mathrm{s}}$       |                  45 | \%                                       | electrolyte volume fraction, separator       |

#### Electrolyte

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
|  ce0                         | $c_{\mathrm{e}0}$             |   1000              | $\mathrm{mol}/\mathrm{m}^3$              | initial electrolyte concentration            |
|  De()                        | $`D_{\mathrm{e}}()`$          |                     | $\mathrm{m}^2/\mathrm{s}$                | diffusivity of lithium ions in the electrolyte|
|  kappa()                     | $`\sigma_{\mathrm{e}}()`$     |                     | $\mathrm{S}/\mathrm{m}$                  | electronic conductivity                       |

#### Cell parameters

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
|  cell_length                 |                                |  1.58               |  $\mathrm{m}$                            | electrode length                             |
|  cell_width                  |                                |  $6.5\times10^{-2}$ |  $\mathrm{m}$                            | electrode width                              |
|  I_typ                       |                                |  5                  |  $\mathrm{A}$                            | typical current                              |
|  brugg                       |                                |  1.5                |  -                                       | Bruggeman coefﬁcients (theoretical value)    |
