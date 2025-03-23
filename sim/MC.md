Finite sample properties
================

For the fractional unobserved components model, this section examines
the finite sample properties of the proposed estimation methods relative
to popular competitors for the data-generating mechanism  
$$
y_t = d_t + x_t + c_t + u_t, \qquad \Delta_+^d x_t = \eta_t, \qquad c_t - b_1 c_{t-1} - b_{2}c_{t-2} = \epsilon_t,
$$ where $d_t$ controls for deterministic terms, $u_t$ is an additional
measurement error, $b_{1,0} = 1.6$, $b_{2,0} = -0.8$ reflect strong
cyclical patterns, and  
$$
    \mathrm{Var} \left(\begin{matrix} \eta_t \\ \epsilon_t\end{matrix} \right) = Q.
$$ The subsections below vary over $Q$, $d_t$, $u_t$. They consider
sample sizes $n \in \{100, 200, 300\}$, integration order
$d_0 \in \{0.75, 1.00, 1.75\}$, and variance ratios  
$$
\nu_0 \in \left\{1, \frac{n^{-1} \sum_{t=1}^{n}\sum_{j=0}^{t-1} \pi_j^2(-d_0)}{\sum_{j=0}^{\infty} a_j^2(\varphi_0)}  r^{-1}\right\}, \qquad r \in \left\{1, 10, 30\right\},
$$ with $\varphi_0 = (b_{1,0}, b_{2,0})'$, and
$a(L, \varphi_0) = \sum_{j=0}^{\infty} a_j(\varphi_0)L^j = (1 - b_{1,0}L - b_{2,0}L^2)^{-1}$.
While the choices for $n$ and $d_0$ cover empirically relevant sample
sizes and integration orders in macroeconomics and finance and allow for
a comparison with the $I(1)$ UC model when $d_0=1$, the choice for
$\nu_0$ is justified as follows: Trivially, setting
$\nu_0 = \sigma_{\epsilon, 0}^2/\sigma_{\eta, 0}^2 = 1$ assigns equal
variation to long- and short-run innovations. By its non-stationary
nature, the trend then dominates the overall variance of $y_t$,
i.e. $\mathrm{Var}(x_t)/\mathrm{Var}(c_t) = O(t^{2d_0 - 1})$, which
constitutes a favorable scenario for estimating $d_0$. At contrast,
letting $\nu_0$ depend on $d_0$ and $\varphi_0$ controls for the
diverging variance ratio: The numerator is the mean variance of the
trend component (under $\sigma_{\eta, 0}^2 = 1$ as in the simulations),
hence
$\mathrm{Var}(c_t) = \mathrm{Var}(a(L, \varphi_0) \epsilon_t) = \nu \sum_{j=0}^{\infty} a_j^2(\varphi_0) = r^{-1} n^{-1} \sum_{t=1}^{n}\sum_{j=0}^{t-1} \pi_j^2(-d_0)$
is proportional to the mean variance of $x_t$, and
$n^{-1} \sum_{t=1}^{n} \mathrm{Var}(x_t)/ \mathrm{Var}(c_t) = r$. This
fixes the variance ratio of $x_t$ and $c_t$ (instead of
$\sigma^2_{\epsilon, 0}/\sigma_{\eta, 0}^2$), and the lower $r$, the
weaker the relative contribution of the trend to the overall variation,
and the less favorable the scenario for estimating $d_0$ and $x_t$.

Each Monte Carlo simulation consists of $1000$ replications. For the QML
estimator, the trend is initialized with variance zero (as implied by
the type II definition of long memory), whereas the cycle is initialized
with its long-run variance. Once the prediction error variance satisfies
$\left| \frac{\mathrm{Var}_\theta(v_{t+1}(\theta)|y_1,...,y_t)-\mathrm{Var}_\theta(v_{t}(\theta)|y_1,...,y_{t-1})}{\mathrm{Var}_\theta(v_{t}(\theta)|y_1,...,y_{t-1})} \right| < 0.01,$
the optimization switches to the steady-state Kalman filter, which
assumes the prediction error variance to be constant from that point on.
Both CSS and QML estimator are initialized by first evaluating their
objective function at a large grid for the model parameters, and the
grid point referring to the lowest value of the CSS objective function
or the lowest negative likelihood is chosen as the starting point for
numerical optimization. As a benchmark, the exact local Whittle
estimator of Shimotsu and Phillips (2005) is introduced, using
$m = \lfloor n^{j} \rfloor$ Fourier frequencies,
$j \in \{.50, .60, .70\}$. Moreover, CSS and QML estimation results from
the $I(1)$ UC model (setting $d=1$) are reported, as well as QML
estimates from the approximate fractional UC model of Hartl and
Jucknewitz (2022). The latter approximates the fractional differencing
operator by an ARMA($3$, $3$) polynomial, which yields a low-dimensional
state space model so that Kalman filter and smoother remain
computationally feasible. All benchmarks are initialized analogously to
the fractional UC model, and starting values are chosen analogously by
evaluating the objective functions at the same grid (but setting $d=1$
for the integer-integrated models). Parameter estimates are compared by
the root mean squared error (RMSE), as well as by the bias. To judge
trend and cycle estimates, the coefficients of determination $R_x^2$ and
$R_c^2$ from regressing $x_t$ and $c_t$ on their estimates from the
Kalman smoother are reported.

# 1 Uncorrelated innovations

The first simulation considers the prototypical fractional UC model with
no measurement error $u_t = 0$, no deterministic terms $d_t=0$, and
diagonal covariance matrix $Q$. The resulting fractional UC model is
thus  
$$
y_t = x_t + c_t, \qquad \Delta_+^d x_t = \eta_t, \qquad c_t - b_1 c_{t-1} - b_{2}c_{t-2} = \epsilon_t, \qquad Q = \begin{pmatrix}
            1 & 0 \\
            0 & \nu
        \end{pmatrix}.
$$

For a better comparison of QML and CSS, the QML simulation assumes
$\sigma_{\eta, 0}^2 = 1$ known, so that optimization is conducted over
$(d, \nu, b_1, b_2)'$ both for CSS and QML estimator.

Table <a href="#tab:tab1">1.1</a> shows the RMSE and the bias for the
estimated integration orders. As can be expected, bias and RMSE decrease
both in $n$ and $r$, and are significantly smaller as compared to the
nonparametric Whittle estimators. The difference is particularly
striking when the signal of the trend is drowned by the cycle (i.e. $r$
small / $\nu_0$ large), which biases the Whittle estimates towards zero,
whereas the estimates for the fractional UC model are hardly affected.
Noticeably, for the fractional UC model QML is slightly superior to CSS,
and using ARMA approximations as suggested by Hartl and Jucknewitz
(2022) yields estimates that are close to the exact fractional UC models
in terms of RMSE and bias.

<details>
<summary>

Table <a href="#tab:tab1">1.1</a>: Simulation with uncorrelated
innovations: root mean squared errors (RMSE) and bias for the
integration order estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab1"></span>Table 1.1:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

RMSE

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

bias

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

-0.36

</td>
<td style="text-align:right;">

0.52

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.60

</td>
<td style="text-align:right;">

-0.56

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

-0.32

</td>
<td style="text-align:right;">

0.45

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.55

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

1.58

</td>
<td style="text-align:right;">

1.36

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.51

</td>
<td style="text-align:right;">

-0.61

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

1.13

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

1.14

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.65

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

0.33

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.50

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.83

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

0.45

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-1.57

</td>
<td style="text-align:right;">

-1.35

</td>
<td style="text-align:right;">

-0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.28

</td>
<td style="text-align:right;">

0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

1.60

</td>
<td style="text-align:right;">

1.58

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-1.29

</td>
<td style="text-align:right;">

-1.11

</td>
<td style="text-align:right;">

-0.22

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-1.12

</td>
<td style="text-align:right;">

-0.96

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

1.35

</td>
<td style="text-align:right;">

1.35

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.57

</td>
<td style="text-align:right;">

-0.55

</td>
<td style="text-align:right;">

-0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.34

</td>
<td style="text-align:right;">

-0.41

</td>
<td style="text-align:right;">

-0.22

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

-0.86

</td>
<td style="text-align:right;">

-0.82

</td>
<td style="text-align:right;">

-0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

1.20

</td>
<td style="text-align:right;">

1.21

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-1.60

</td>
<td style="text-align:right;">

-1.58

</td>
<td style="text-align:right;">

-0.95

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.53

</td>
<td style="text-align:right;">

-0.50

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-1.34

</td>
<td style="text-align:right;">

-1.34

</td>
<td style="text-align:right;">

-0.88

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.34

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-1.19

</td>
<td style="text-align:right;">

-1.20

</td>
<td style="text-align:right;">

-0.82

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

-0.57

</td>
<td style="text-align:right;">

-0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

1.55

</td>
<td style="text-align:right;">

1.63

</td>
<td style="text-align:right;">

1.27

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.41

</td>
<td style="text-align:right;">

-0.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

-0.77

</td>
<td style="text-align:right;">

-0.88

</td>
<td style="text-align:right;">

-0.48

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.31

</td>
<td style="text-align:right;">

-1.54

</td>
<td style="text-align:right;">

-1.63

</td>
<td style="text-align:right;">

-1.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

1.28

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

1.14

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.44

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-1.27

</td>
<td style="text-align:right;">

-1.41

</td>
<td style="text-align:right;">

-1.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

-0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

1.14

</td>
<td style="text-align:right;">

1.29

</td>
<td style="text-align:right;">

1.05

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-1.12

</td>
<td style="text-align:right;">

-1.28

</td>
<td style="text-align:right;">

-1.04

</td>
</tr>
</tbody>
</table>
</details>

  

Tables <a href="#tab:tab2">1.2</a> and <a href="#tab:tab3">1.3</a>
detail RMSE and bias for $\nu_0$ and the autoregressive parameters. In
addition to the fractional UC model, the tables also display the
estimation results for the $I(1)$ UC benchmark. For $b_{1,0}$ and
$b_{2,0}$, CSS and QML estimates for the fractional UC model behave
equally well, while those from the approximate fractional UC model and
the integer-integrated benchmarks come with a slightly higher RMSE.
Interestingly, the estimates for $\nu_0$ show a clear dominance of QML
over CSS, as the latter comes with a much higher RMSE and strong,
positive bias. QML - at contrast - does not appear biased. A violation
$d_0 \neq 1$ in integer-integrated models has a rather small effect on
the estimates for the cyclical autoregressive coefficients, but a
comparably strong effect on the estimate for $\nu_0$ also for the QML
estimator, thus shifting variation from trend to cycle.

<details>
<summary>

Table <a href="#tab:tab2">1.2</a>: Simulation with uncorrelated
innovations: root mean squared errors (RMSE) for the other parameter
estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab2"></span>Table 1.2:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

83463.24

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

77405.62

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

70230.62

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

28468.22

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

574.98

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

24534.57

</td>
<td style="text-align:right;">

273.26

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

80051.01

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

49436.47

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

164887.90

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

54935.16

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

295220.34

</td>
<td style="text-align:right;">

153.13

</td>
<td style="text-align:right;">

148.49

</td>
<td style="text-align:right;">

213664.23

</td>
<td style="text-align:right;">

312.86

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

39234.37

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

14237.66

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.38

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

39021.02

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

13767.54

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

111570.51

</td>
<td style="text-align:right;">

18.52

</td>
<td style="text-align:right;">

15.09

</td>
<td style="text-align:right;">

73987.46

</td>
<td style="text-align:right;">

196.75

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

15956.11

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

24703.56

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

1.35

</td>
<td style="text-align:right;">

1.44

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.44

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

7165.15

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

14749.44

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

67151.54

</td>
<td style="text-align:right;">

6.87

</td>
<td style="text-align:right;">

6.70

</td>
<td style="text-align:right;">

50228.68

</td>
<td style="text-align:right;">

260.83

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.11

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

9549.81

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

74011.20

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

32785.57

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

18972.56

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

21.97

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

35411.76

</td>
<td style="text-align:right;">

373.03

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

45783.49

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

61783.03

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

211139.75

</td>
<td style="text-align:right;">

4.62

</td>
<td style="text-align:right;">

1.59

</td>
<td style="text-align:right;">

70366.61

</td>
<td style="text-align:right;">

2.18

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

445061.91

</td>
<td style="text-align:right;">

3386.10

</td>
<td style="text-align:right;">

600.24

</td>
<td style="text-align:right;">

116294.45

</td>
<td style="text-align:right;">

13791.91

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

31647.18

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

28358.00

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

24604.27

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

2200.72

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

207156.11

</td>
<td style="text-align:right;">

71.21

</td>
<td style="text-align:right;">

60.16

</td>
<td style="text-align:right;">

30899.04

</td>
<td style="text-align:right;">

241.31

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

462.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

29009.51

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

1.14

</td>
<td style="text-align:right;">

1.22

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.41

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

413.71

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

14432.61

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

107572.01

</td>
<td style="text-align:right;">

27.63

</td>
<td style="text-align:right;">

26.51

</td>
<td style="text-align:right;">

7810.52

</td>
<td style="text-align:right;">

302.94

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1658.75

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

43705.69

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

725.27

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

5238.88

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

10.34

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

4933.27

</td>
<td style="text-align:right;">

448.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

2188.36

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

62138.43

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

160382.22

</td>
<td style="text-align:right;">

2.82

</td>
<td style="text-align:right;">

1.77

</td>
<td style="text-align:right;">

106537.07

</td>
<td style="text-align:right;">

1.88

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

470811.39

</td>
<td style="text-align:right;">

13285.17

</td>
<td style="text-align:right;">

1805.85

</td>
<td style="text-align:right;">

83163.77

</td>
<td style="text-align:right;">

27792.57

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

143.77

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

17542.51

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

653.98

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

1083.51

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

299400.82

</td>
<td style="text-align:right;">

165.09

</td>
<td style="text-align:right;">

133.76

</td>
<td style="text-align:right;">

1617.94

</td>
<td style="text-align:right;">

703.29

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

31527.73

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

18430.29

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

1.05

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.44

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

2.05

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

1.02

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

155494.38

</td>
<td style="text-align:right;">

59.27

</td>
<td style="text-align:right;">

46.85

</td>
<td style="text-align:right;">

532.65

</td>
<td style="text-align:right;">

253.75

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
</tbody>
</table>
</details>

  

<details>
<summary>

Table <a href="#tab:tab3">1.3</a>: Simulation with uncorrelated
innovations: bias for the other parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab3"></span>Table 1.3:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

8657.12

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

15815.48

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

7825.77

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

3048.26

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

37.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

1866.07

</td>
<td style="text-align:right;">

149.86

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

7283.78

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

8468.21

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

35291.94

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

8828.59

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

118564.40

</td>
<td style="text-align:right;">

25.24

</td>
<td style="text-align:right;">

-21.74

</td>
<td style="text-align:right;">

54801.47

</td>
<td style="text-align:right;">

204.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

1896.53

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

1295.86

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

2971.95

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

1183.22

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

24839.93

</td>
<td style="text-align:right;">

7.65

</td>
<td style="text-align:right;">

-2.22

</td>
<td style="text-align:right;">

8068.99

</td>
<td style="text-align:right;">

111.38

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

518.44

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

2296.62

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

495.39

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

1168.09

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

8422.09

</td>
<td style="text-align:right;">

2.47

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

5954.43

</td>
<td style="text-align:right;">

146.51

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

371.73

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

8846.97

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1660.48

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

649.94

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

4.10

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

1967.31

</td>
<td style="text-align:right;">

244.18

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

2374.59

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

7116.17

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

49749.51

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

9136.07

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

215127.38

</td>
<td style="text-align:right;">

214.94

</td>
<td style="text-align:right;">

-57.78

</td>
<td style="text-align:right;">

15728.76

</td>
<td style="text-align:right;">

1764.81

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

1068.13

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

1080.56

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

1034.93

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

139.18

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

59375.80

</td>
<td style="text-align:right;">

33.63

</td>
<td style="text-align:right;">

-5.97

</td>
<td style="text-align:right;">

995.95

</td>
<td style="text-align:right;">

208.24

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

18.56

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

1879.45

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

22.40

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

773.54

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

20082.20

</td>
<td style="text-align:right;">

17.04

</td>
<td style="text-align:right;">

-1.74

</td>
<td style="text-align:right;">

141.41

</td>
<td style="text-align:right;">

227.05

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

67.04

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

3487.79

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

24.49

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

166.30

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1.93

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

267.26

</td>
<td style="text-align:right;">

324.61

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

84.58

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

5875.34

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

29507.76

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.73

</td>
<td style="text-align:right;">

15376.12

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

237663.70

</td>
<td style="text-align:right;">

1441.07

</td>
<td style="text-align:right;">

225.80

</td>
<td style="text-align:right;">

-4223.49

</td>
<td style="text-align:right;">

5284.73

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

5.90

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

925.16

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

25.71

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

44.23

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

106265.64

</td>
<td style="text-align:right;">

71.73

</td>
<td style="text-align:right;">

-16.15

</td>
<td style="text-align:right;">

-1587.94

</td>
<td style="text-align:right;">

555.97

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

998.65

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

894.31

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

33898.88

</td>
<td style="text-align:right;">

35.92

</td>
<td style="text-align:right;">

-3.85

</td>
<td style="text-align:right;">

-532.59

</td>
<td style="text-align:right;">

230.71

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
</tbody>
</table>
</details>

  

Table <a href="#tab:tab4">1.4</a> compares the estimates for $x_t$ and
$c_t$ for the different models by regressing the respective Kalman
smoother-based estimates on the true trend and cycle and reporting the
coefficient of determination. Note that the coefficient of determination
for $x_t$ should be interpreted with caution: Under correct
specification, $\hat x_t$ and $x_t$ are cointegrated, so that regressing
the estimate on the true, simulated $x_t$ yields an interpretable $R^2$.
This no longer holds under model misspecification, e.g. when
$x_t \sim I(d_0)$, $d_0 \neq 1$, and $x_t$ is regressed on an $I(1)$
trend: The respective regression then becomes spurious, and the $R^2$ is
no longer interpretable. Consequently, a high $R^2$ for
integer-integrated models when $d_0 \neq 1$ does not necessarily
indicate a good fit. Fortunately, the simulated $c_t$ are always $I(0)$,
so that their coefficients of determination allow for a valid comparison
among the fractional and integer-integrated UC models whenever
$d_0 \neq 1$.

As can be seen, differences between the coefficients of determination
are almost negligible for CSS and QML estimator of the fractional UC
model, with the latter exhibiting slightly larger coefficients of
determination. Strikingly, for $d_0=1$ the fractional UC model shows no
loss in efficiency compared to the $I(1)$ UC model. For non-integer
$d_0$, the fractional model shows a higher $R^2$ than the
integer-integrated benchmarks, particularly when $d_0 = 1.75$. However,
integer-integrated UC models often provide a good approximation to
fractionally integrated trends as long as correlation between trend and
cycle innovations is ruled out.

<details>
<summary>

Table <a href="#tab:tab4">1.4</a>: Simulation with uncorrelated
innovations: Coefficient of determination from regressing true trend and
cycle on their respective estimates from the Kalman smoother.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab4"></span>Table 1.4:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Trend

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Cycle

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.80

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.78

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.79

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.65

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.78

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.81

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.61

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.78

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.41

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.86

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.79

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.81

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.69

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.49

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.42

</td>
</tr>
</tbody>
</table>
</details>

  

# 2 Correlated innovations

The second simulation generalizes the prototypical fractional UC model
to correlated innovations, allowing for a non-diagonal covariance matrix
$Q$. The resulting fractional UC model is thus  
$$
y_t = x_t + c_t, \qquad \Delta_+^d x_t = \eta_t, \qquad c_t - b_1 c_{t-1} - b_{2}c_{t-2} = \epsilon_t, \qquad Q = \left( \begin{matrix}
            1 & \rho \sqrt{\nu} \\
            \rho\sqrt{\nu} & \nu
        \end{matrix} \right),
$$ where $\rho = -0.8$ is set to mimic strong but not perfect
correlation. Optimization is now conducted over
$(d, \nu, \nu_2, b_1, b_2)'$ for the CSS estimator, and
$(d, \sigma_{\eta}^2, \sigma_{\eta \epsilon}, \sigma_\epsilon^2, b_1, b_2)'$
for the QML estimator of the fractional UC model, however in both cases
bias and RMSE are reported for
$\hat \rho = \widehat{\mathrm{Corr}(\eta_t, \epsilon_t)}$ for
comprehensibility.

Table <a href="#tab:tab21">2.1</a> shows the RMSE and the bias for the
estimated integration orders. As before, bias and RMSE decrease both in
$n$ and $r$, and are significantly smaller as compared to the
nonparametric Whittle estimators. Compared to the uncorrelated case, the
threshold between CSS and QML in terms of performance is somewhat
smaller, which may be due to the additional parameter that QML now has
to estimate. Nonetheless, QML is clearly favored by the Monte Carlo
study, also because its bias is comparably small.

<details>
<summary>

Table <a href="#tab:tab21">2.1</a>: Simulation with correlated
innovations: root mean squared errors (RMSE) and bias for the
integration order estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab21"></span>Table 2.1:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

RMSE

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

bias

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

-0.37

</td>
<td style="text-align:right;">

0.64

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

0.58

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.40

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.50

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.57

</td>
<td style="text-align:right;">

-0.31

</td>
<td style="text-align:right;">

0.53

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.81

</td>
<td style="text-align:right;">

-0.52

</td>
<td style="text-align:right;">

0.65

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

1.55

</td>
<td style="text-align:right;">

1.30

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

-0.73

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-1.54

</td>
<td style="text-align:right;">

-1.29

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.34

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.36

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

1.27

</td>
<td style="text-align:right;">

1.05

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

-0.68

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-1.25

</td>
<td style="text-align:right;">

-1.04

</td>
<td style="text-align:right;">

-0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

1.11

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

-0.51

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-1.08

</td>
<td style="text-align:right;">

-0.88

</td>
<td style="text-align:right;">

-0.07

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.56

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.53

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

0.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

-0.55

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.84

</td>
<td style="text-align:right;">

-0.80

</td>
<td style="text-align:right;">

-0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

1.58

</td>
<td style="text-align:right;">

1.56

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

-0.71

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-1.57

</td>
<td style="text-align:right;">

-1.55

</td>
<td style="text-align:right;">

-0.92

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.47

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

1.33

</td>
<td style="text-align:right;">

1.32

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

-0.61

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-1.32

</td>
<td style="text-align:right;">

-1.31

</td>
<td style="text-align:right;">

-0.83

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

1.18

</td>
<td style="text-align:right;">

1.18

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-1.17

</td>
<td style="text-align:right;">

-1.17

</td>
<td style="text-align:right;">

-0.75

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.47

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

-0.52

</td>
<td style="text-align:right;">

-0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.36

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.61

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.77

</td>
<td style="text-align:right;">

-0.87

</td>
<td style="text-align:right;">

-0.45

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

1.52

</td>
<td style="text-align:right;">

1.62

</td>
<td style="text-align:right;">

1.25

</td>
<td style="text-align:right;">

-0.64

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-1.51

</td>
<td style="text-align:right;">

-1.61

</td>
<td style="text-align:right;">

-1.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.40

</td>
<td style="text-align:right;">

-0.54

</td>
<td style="text-align:right;">

-0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

1.26

</td>
<td style="text-align:right;">

1.40

</td>
<td style="text-align:right;">

1.11

</td>
<td style="text-align:right;">

-0.75

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-1.25

</td>
<td style="text-align:right;">

-1.39

</td>
<td style="text-align:right;">

-1.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.31

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

1.12

</td>
<td style="text-align:right;">

1.27

</td>
<td style="text-align:right;">

1.01

</td>
<td style="text-align:right;">

-0.56

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-1.11

</td>
<td style="text-align:right;">

-1.27

</td>
<td style="text-align:right;">

-1.01

</td>
</tr>
</tbody>
</table>
</details>

  

Tables <a href="#tab:tab22">2.2</a> to <a href="#tab:tab25">2.5</a>
detail RMSE and bias for the covariance and autoregressive parameters,
both for the fractional UC model and the integer-integrated benchmark.
As before, estimates for $\nu_0$ from the CSS estimator come with a
large RMSE and bias, so they should be interpreted with caution.
Moreover, results from the QML estimator show that the covariance matrix
$Q_0$ is notoriously difficult to estimate, in particular when $n$ is
small and $\nu_0$ is large. For larger $n$, both exact and approximate
fractional UC model yield a small bias for the parameters in $Q_0$ even
when $\nu_0$ is large, and the RMSE for $\sigma_{\epsilon, 0}^2 = \nu_0$
is at least moderate compared to the level of $\nu_0$, unlike the RMSE
for $\sigma_{\eta, 0}^2 = 1$ which is high relative to
$\sigma_{\eta, 0}^2$ whenever $\nu_0$ is large. Noticeably, estimates
$\hat \rho$ from the integer-integrated UC models are massively biased
whenever $d_0 \neq 1$, and often were found to converge to corner
solutions (where $|\hat \rho | = 1$). Hence, empirical results in the UC
literature that indicate (almost) perfect correlation between long- and
short-run innovations may also be an artifact generated by the
misspecification of the integration order of the UC model. Results for
the autoregressive parameters are similar to those from the uncorrelated
model, where all fractional and integer-integrated UC models showed
stable and reliable behavior.

<details>
<summary>

Table <a href="#tab:tab22">2.2</a>: Simulation with correlated
innovations: root mean squared errors (RMSE) for the covariance
parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab22"></span>Table 2.2:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\rho$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\rho$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\rho$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\rho^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\rho^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\eta}^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\eta}^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\sigma_{\eta}^{2^{I(1)}}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\epsilon}^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\epsilon}^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\sigma_{\epsilon}^{2^{I(1)}}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

341.93

</td>
<td style="text-align:right;">

373.53

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

1.38

</td>
<td style="text-align:right;">

1.35

</td>
<td style="text-align:right;">

1.41

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

1.18

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.48

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

218.35

</td>
<td style="text-align:right;">

111.15

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

2.77

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

14.14

</td>
<td style="text-align:right;">

1.07

</td>
<td style="text-align:right;">

76.45

</td>
<td style="text-align:right;">

13.67

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

802.30

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

187.72

</td>
<td style="text-align:right;">

328.32

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

1.55

</td>
<td style="text-align:right;">

1.48

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.43

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

272.11

</td>
<td style="text-align:right;">

106.32

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

1.52

</td>
<td style="text-align:right;">

1.16

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

1.79

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

1.29

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

758.96

</td>
<td style="text-align:right;">

1000.81

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

1.11

</td>
<td style="text-align:right;">

75.60

</td>
<td style="text-align:right;">

25.20

</td>
<td style="text-align:right;">

282.06

</td>
<td style="text-align:right;">

352.43

</td>
<td style="text-align:right;">

218.11

</td>
<td style="text-align:right;">

867.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

1.47

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

2.19

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

1.75

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

63.46

</td>
<td style="text-align:right;">

55.05

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

2.34

</td>
<td style="text-align:right;">

1.20

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

2.52

</td>
<td style="text-align:right;">

1.32

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

231.43

</td>
<td style="text-align:right;">

116.80

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

23.36

</td>
<td style="text-align:right;">

6.78

</td>
<td style="text-align:right;">

53.57

</td>
<td style="text-align:right;">

60.60

</td>
<td style="text-align:right;">

25.28

</td>
<td style="text-align:right;">

470.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

1.49

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

5.97

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

5.29

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

31.63

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

1.95

</td>
<td style="text-align:right;">

2.82

</td>
<td style="text-align:right;">

3.64

</td>
<td style="text-align:right;">

1.73

</td>
<td style="text-align:right;">

2.73

</td>
<td style="text-align:right;">

3.58

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

96.00

</td>
<td style="text-align:right;">

55.01

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

8.00

</td>
<td style="text-align:right;">

2.63

</td>
<td style="text-align:right;">

16.12

</td>
<td style="text-align:right;">

19.63

</td>
<td style="text-align:right;">

9.86

</td>
<td style="text-align:right;">

681.40

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

265.29

</td>
<td style="text-align:right;">

283.25

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

1.45

</td>
<td style="text-align:right;">

1.47

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.37

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

119.60

</td>
<td style="text-align:right;">

7.85

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.35

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1.55

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

6.58

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

94.79

</td>
<td style="text-align:right;">

6.22

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

2958.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

224.34

</td>
<td style="text-align:right;">

280.91

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

1.49

</td>
<td style="text-align:right;">

1.52

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.36

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

199.40

</td>
<td style="text-align:right;">

99.99

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

1.36

</td>
<td style="text-align:right;">

1.48

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

1.58

</td>
<td style="text-align:right;">

1.53

</td>
<td style="text-align:right;">

1.52

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

5167.60

</td>
<td style="text-align:right;">

5772.34

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

1.34

</td>
<td style="text-align:right;">

258.24

</td>
<td style="text-align:right;">

99.42

</td>
<td style="text-align:right;">

713.42

</td>
<td style="text-align:right;">

1082.93

</td>
<td style="text-align:right;">

723.99

</td>
<td style="text-align:right;">

3313.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

8.07

</td>
<td style="text-align:right;">

6.03

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

493.73

</td>
<td style="text-align:right;">

580.56

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

80.56

</td>
<td style="text-align:right;">

13.27

</td>
<td style="text-align:right;">

163.06

</td>
<td style="text-align:right;">

242.61

</td>
<td style="text-align:right;">

78.28

</td>
<td style="text-align:right;">

1735.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

1.23

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

1.20

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

1.26

</td>
<td style="text-align:right;">

0.65

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

319.78

</td>
<td style="text-align:right;">

193.53

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

19.14

</td>
<td style="text-align:right;">

4.95

</td>
<td style="text-align:right;">

72.49

</td>
<td style="text-align:right;">

69.22

</td>
<td style="text-align:right;">

28.07

</td>
<td style="text-align:right;">

2491.60

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

171.32

</td>
<td style="text-align:right;">

187.49

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

1.49

</td>
<td style="text-align:right;">

1.50

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

78.29

</td>
<td style="text-align:right;">

2.39

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1.12

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

4.92

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

88.99

</td>
<td style="text-align:right;">

4.78

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

4439.58

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

229.49

</td>
<td style="text-align:right;">

212.40

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

1.38

</td>
<td style="text-align:right;">

1.37

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.35

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

136.68

</td>
<td style="text-align:right;">

55.76

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

1.78

</td>
<td style="text-align:right;">

2.08

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

1.99

</td>
<td style="text-align:right;">

1.80

</td>
<td style="text-align:right;">

1.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

15173.55

</td>
<td style="text-align:right;">

15978.76

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

1.28

</td>
<td style="text-align:right;">

567.33

</td>
<td style="text-align:right;">

53.39

</td>
<td style="text-align:right;">

1120.05

</td>
<td style="text-align:right;">

2185.04

</td>
<td style="text-align:right;">

1517.30

</td>
<td style="text-align:right;">

5638.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

1.56

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

1.23

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

135.42

</td>
<td style="text-align:right;">

2.64

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

1077.27

</td>
<td style="text-align:right;">

1599.49

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

148.44

</td>
<td style="text-align:right;">

9.82

</td>
<td style="text-align:right;">

389.28

</td>
<td style="text-align:right;">

342.63

</td>
<td style="text-align:right;">

163.53

</td>
<td style="text-align:right;">

2298.89

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

498.43

</td>
<td style="text-align:right;">

533.28

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

38.24

</td>
<td style="text-align:right;">

4.88

</td>
<td style="text-align:right;">

144.07

</td>
<td style="text-align:right;">

127.66

</td>
<td style="text-align:right;">

56.43

</td>
<td style="text-align:right;">

2877.87

</td>
</tr>
</tbody>
</table>
</details>
<details>
<summary>

Table <a href="#tab:tab23">2.3</a>: Simulation with correlated
innovations: bias for the covariance parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab23"></span>Table 2.3:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\rho$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\rho$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\rho$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\rho^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\rho^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\eta}^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\eta}^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\sigma_{\eta}^{2^{I(1)}}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\epsilon}^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\sigma_{\epsilon}^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\sigma_{\epsilon}^{2^{I(1)}}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

127.52

</td>
<td style="text-align:right;">

165.94

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

1.15

</td>
<td style="text-align:right;">

1.13

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.73

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

55.08

</td>
<td style="text-align:right;">

16.60

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

2.15

</td>
<td style="text-align:right;">

-0.32

</td>
<td style="text-align:right;">

31.36

</td>
<td style="text-align:right;">

1.91

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

278.62

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

37.98

</td>
<td style="text-align:right;">

116.68

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

1.39

</td>
<td style="text-align:right;">

1.28

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

117.51

</td>
<td style="text-align:right;">

24.89

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

-584.37

</td>
<td style="text-align:right;">

-974.58

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

27.82

</td>
<td style="text-align:right;">

8.38

</td>
<td style="text-align:right;">

153.69

</td>
<td style="text-align:right;">

138.06

</td>
<td style="text-align:right;">

44.82

</td>
<td style="text-align:right;">

245.53

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

4.41

</td>
<td style="text-align:right;">

3.65

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

-7.35

</td>
<td style="text-align:right;">

-96.24

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

7.61

</td>
<td style="text-align:right;">

1.71

</td>
<td style="text-align:right;">

31.09

</td>
<td style="text-align:right;">

28.04

</td>
<td style="text-align:right;">

2.15

</td>
<td style="text-align:right;">

186.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

1.24

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.66

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

1.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.42

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

-31.43

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

2.70

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

8.04

</td>
<td style="text-align:right;">

8.29

</td>
<td style="text-align:right;">

-1.24

</td>
<td style="text-align:right;">

236.92

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

75.55

</td>
<td style="text-align:right;">

104.31

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

1.38

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.84

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

17.81

</td>
<td style="text-align:right;">

1.18

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

-0.50

</td>
<td style="text-align:right;">

38.16

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

1451.91

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

54.57

</td>
<td style="text-align:right;">

99.41

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

1.36

</td>
<td style="text-align:right;">

1.44

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.82

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

83.36

</td>
<td style="text-align:right;">

22.38

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

-5112.83

</td>
<td style="text-align:right;">

-5764.63

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

55.85

</td>
<td style="text-align:right;">

13.88

</td>
<td style="text-align:right;">

381.10

</td>
<td style="text-align:right;">

401.26

</td>
<td style="text-align:right;">

95.16

</td>
<td style="text-align:right;">

153.83

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

-105.97

</td>
<td style="text-align:right;">

-580.27

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

16.32

</td>
<td style="text-align:right;">

2.78

</td>
<td style="text-align:right;">

115.88

</td>
<td style="text-align:right;">

97.74

</td>
<td style="text-align:right;">

12.20

</td>
<td style="text-align:right;">

787.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

-2.98

</td>
<td style="text-align:right;">

-193.39

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

5.59

</td>
<td style="text-align:right;">

1.03

</td>
<td style="text-align:right;">

45.94

</td>
<td style="text-align:right;">

36.74

</td>
<td style="text-align:right;">

2.41

</td>
<td style="text-align:right;">

1083.28

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

33.51

</td>
<td style="text-align:right;">

57.61

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

1.40

</td>
<td style="text-align:right;">

1.44

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.84

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

7.76

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

37.56

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

-0.37

</td>
<td style="text-align:right;">

2554.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

60.40

</td>
<td style="text-align:right;">

77.62

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

1.26

</td>
<td style="text-align:right;">

1.29

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.86

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

59.78

</td>
<td style="text-align:right;">

17.27

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

-15108.59

</td>
<td style="text-align:right;">

-15972.50

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

88.55

</td>
<td style="text-align:right;">

14.01

</td>
<td style="text-align:right;">

687.92

</td>
<td style="text-align:right;">

724.83

</td>
<td style="text-align:right;">

176.29

</td>
<td style="text-align:right;">

-501.90

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

22.40

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

-856.50

</td>
<td style="text-align:right;">

-1598.85

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

28.36

</td>
<td style="text-align:right;">

2.62

</td>
<td style="text-align:right;">

291.23

</td>
<td style="text-align:right;">

172.16

</td>
<td style="text-align:right;">

22.48

</td>
<td style="text-align:right;">

1205.65

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

-55.19

</td>
<td style="text-align:right;">

-533.08

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

8.73

</td>
<td style="text-align:right;">

1.08

</td>
<td style="text-align:right;">

97.87

</td>
<td style="text-align:right;">

74.22

</td>
<td style="text-align:right;">

4.73

</td>
<td style="text-align:right;">

1417.96

</td>
</tr>
</tbody>
</table>
</details>
<details>
<summary>

Table <a href="#tab:tab26">2.4</a>: Simulation with correlated
innovations: root mean squared errors (RMSE) for the autoregressive
parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab26"></span>Table 2.4:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.41

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.17

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.38

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.31

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
</tbody>
</table>
</details>
<details>
<summary>

Table <a href="#tab:tab25">2.5</a>: Simulation with correlated
innovations: bias for the autoregressive parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab25"></span>Table 2.5:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.54

</td>
<td style="text-align:right;">

-0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.55

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.53

</td>
<td style="text-align:right;">

-0.22

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
</tbody>
</table>
</details>

Table <a href="#tab:tab24">2.6</a> compares the coefficient of
determination for the estimates of $x_t$ and $c_t$ by regressing the
respective Kalman smoother-based estimates on the true trend and cycle.
As before, the QML estimates for the fractional UC model dominate,
however the gap to CSS gets somewhat bigger. Moreover, once correlation
is allowed, integer-integrated UC models perform poorer in approximating
the fractional model  
whenever $d_0 \neq 1$, as manifested by a lower $R^2$ compared to the
uncorrelated simulations.

<details>
<summary>

Table <a href="#tab:tab24">2.6</a>: Simulation with correlated
innovations: Coefficient of determination from regressing true trend and
cycle on their respective estimates from the Kalman smoother.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab24"></span>Table 2.6:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Trend

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Cycle

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.64

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.71

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.35

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.59

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.54

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.71

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.86

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.72

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.60

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.42

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.82

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.58

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.76

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.88

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.75

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.44

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.65

</td>
</tr>
</tbody>
</table>
</details>

  

# 3 Deterministic trends

The third simulation adds a linear, deterministic trend to the
fractional UC model $$
y_t =  t + x_t + c_t, \qquad \Delta_+^d x_t = \eta_t, \qquad c_t - b_1 c_{t-1} - b_{2}c_{t-2} = \epsilon_t, \qquad  Q = \begin{pmatrix}
            1 & 0 \\
            0 & \nu
        \end{pmatrix},
$$ Note that since $\sigma_\eta^2 = 1$, the variation generated by the
deterministic trend is proportional to the variation generated by the
stochastic trend innovations.

Table <a href="#tab:tab31">3.1</a> shows the RMSE and the bias for the
estimated integration orders. Noticeably, including a linear
deterministic trend generates a negative bias for $\hat d$ when $r=1$.
In this setting, the overall variation generated by the cycle is large,
whereas the stochastic trend attributes little to the overall dynamics
and is relatively smooth. The fractional UC model then attributes some
of the variation of the stochastic trend to the deterministic term, so
that the estimate for the trend component can take up some of the
variation generated by the cycle. This reduces the memory estimate
$\hat d$ and explains the results for $r=1$. This bias vanishes as $n$
and $r$ increase. All other conclusions are similar to those of table
<a href="#tab:tab1">1.1</a>.

<details>
<summary>

Table <a href="#tab:tab31">3.1</a>: Simulation with deterministic trend:
root mean squared errors (RMSE) and bias for the integration order
estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab31"></span>Table 3.1:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

RMSE

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

bias

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.70

</td>
<td style="text-align:right;">

-0.43

</td>
<td style="text-align:right;">

0.52

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.81

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

0.33

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.51

</td>
<td style="text-align:right;">

-0.37

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.67

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

0.45

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.95

</td>
<td style="text-align:right;">

-0.67

</td>
<td style="text-align:right;">

0.45

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

1.04

</td>
<td style="text-align:right;">

1.07

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

1.73

</td>
<td style="text-align:right;">

1.52

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

-1.00

</td>
<td style="text-align:right;">

-0.96

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-1.73

</td>
<td style="text-align:right;">

-1.50

</td>
<td style="text-align:right;">

-0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.35

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

1.61

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

-0.61

</td>
<td style="text-align:right;">

-0.44

</td>
<td style="text-align:right;">

-0.31

</td>
<td style="text-align:right;">

-1.59

</td>
<td style="text-align:right;">

-1.28

</td>
<td style="text-align:right;">

-0.22

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.39

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

1.13

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

-0.39

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-1.39

</td>
<td style="text-align:right;">

-1.10

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.69

</td>
<td style="text-align:right;">

-0.63

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.68

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

-0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.32

</td>
<td style="text-align:right;">

-0.40

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.68

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.96

</td>
<td style="text-align:right;">

-0.92

</td>
<td style="text-align:right;">

-0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

1.01

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

1.74

</td>
<td style="text-align:right;">

1.73

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

-0.95

</td>
<td style="text-align:right;">

-0.84

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-1.74

</td>
<td style="text-align:right;">

-1.73

</td>
<td style="text-align:right;">

-0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.32

</td>
<td style="text-align:right;">

-0.25

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.63

</td>
<td style="text-align:right;">

-0.56

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

1.63

</td>
<td style="text-align:right;">

1.59

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

-0.31

</td>
<td style="text-align:right;">

-0.46

</td>
<td style="text-align:right;">

-1.62

</td>
<td style="text-align:right;">

-1.58

</td>
<td style="text-align:right;">

-0.90

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.40

</td>
<td style="text-align:right;">

-0.35

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

1.46

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-1.43

</td>
<td style="text-align:right;">

-1.40

</td>
<td style="text-align:right;">

-0.84

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

-0.67

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.47

</td>
<td style="text-align:right;">

-0.64

</td>
<td style="text-align:right;">

-0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.39

</td>
<td style="text-align:right;">

-0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.68

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

-0.90

</td>
<td style="text-align:right;">

-0.97

</td>
<td style="text-align:right;">

-0.49

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

1.73

</td>
<td style="text-align:right;">

1.75

</td>
<td style="text-align:right;">

1.32

</td>
<td style="text-align:right;">

-0.92

</td>
<td style="text-align:right;">

-0.70

</td>
<td style="text-align:right;">

-0.64

</td>
<td style="text-align:right;">

-1.72

</td>
<td style="text-align:right;">

-1.74

</td>
<td style="text-align:right;">

-1.31

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.32

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.50

</td>
<td style="text-align:right;">

-0.66

</td>
<td style="text-align:right;">

-0.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

1.54

</td>
<td style="text-align:right;">

1.65

</td>
<td style="text-align:right;">

1.20

</td>
<td style="text-align:right;">

-0.36

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-1.52

</td>
<td style="text-align:right;">

-1.64

</td>
<td style="text-align:right;">

-1.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.31

</td>
<td style="text-align:right;">

-0.45

</td>
<td style="text-align:right;">

-0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

1.34

</td>
<td style="text-align:right;">

1.51

</td>
<td style="text-align:right;">

1.11

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-1.32

</td>
<td style="text-align:right;">

-1.50

</td>
<td style="text-align:right;">

-1.10

</td>
</tr>
</tbody>
</table>
</details>

  

Tables <a href="#tab:tab32">3.2</a> and <a href="#tab:tab33">3.3</a>
detail RMSE and bias for $\nu_0$ and the autoregressive parameters. They
show a similar performance as in tables <a href="#tab:tab2">1.2</a> and
<a href="#tab:tab3">1.3</a>.

<details>
<summary>

Table <a href="#tab:tab32">3.2</a>: Simulation with deterministic trend:
root mean squared errors (RMSE) for the other parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab32"></span>Table 3.2:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

37084.68

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

233015.34

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

4328.72

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

134455.30

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1451.85

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

31278.73

</td>
<td style="text-align:right;">

23.79

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

3354.90

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

204318.78

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

14800.41

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

2.13

</td>
<td style="text-align:right;">

161834.02

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

142022.84

</td>
<td style="text-align:right;">

148.49

</td>
<td style="text-align:right;">

683.74

</td>
<td style="text-align:right;">

623808.94

</td>
<td style="text-align:right;">

150.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

1742.43

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

99553.79

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

1.32

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.35

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

528.81

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

93115.62

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

124092.42

</td>
<td style="text-align:right;">

17.51

</td>
<td style="text-align:right;">

90.21

</td>
<td style="text-align:right;">

428174.88

</td>
<td style="text-align:right;">

22.63

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

73.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

126831.51

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

1.41

</td>
<td style="text-align:right;">

1.34

</td>
<td style="text-align:right;">

1.45

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.44

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

1579.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

71171.36

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

67250.57

</td>
<td style="text-align:right;">

6.82

</td>
<td style="text-align:right;">

28.41

</td>
<td style="text-align:right;">

313352.90

</td>
<td style="text-align:right;">

10.71

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

28202.98

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

224140.37

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

2045.36

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

80328.72

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

46.83

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

33055.67

</td>
<td style="text-align:right;">

142.21

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.22

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

133.08

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

232532.84

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

106984.49

</td>
<td style="text-align:right;">

7.50

</td>
<td style="text-align:right;">

4.91

</td>
<td style="text-align:right;">

372289.05

</td>
<td style="text-align:right;">

3.10

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

184100.12

</td>
<td style="text-align:right;">

10156.18

</td>
<td style="text-align:right;">

3893.20

</td>
<td style="text-align:right;">

375651.84

</td>
<td style="text-align:right;">

1371.57

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

2.34

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

80203.22

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

2497.53

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

73147.33

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

227176.14

</td>
<td style="text-align:right;">

69.01

</td>
<td style="text-align:right;">

527.52

</td>
<td style="text-align:right;">

265816.66

</td>
<td style="text-align:right;">

108.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

140368.21

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

1.16

</td>
<td style="text-align:right;">

1.10

</td>
<td style="text-align:right;">

1.26

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

246.52

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

50468.26

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

174428.75

</td>
<td style="text-align:right;">

27.37

</td>
<td style="text-align:right;">

20.46

</td>
<td style="text-align:right;">

183430.42

</td>
<td style="text-align:right;">

64.35

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

8926.24

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

176352.05

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1319.48

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

68638.24

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

20.72

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

20484.15

</td>
<td style="text-align:right;">

266.76

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.36

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

128.19

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

194688.71

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

115778.54

</td>
<td style="text-align:right;">

6.61

</td>
<td style="text-align:right;">

1.45

</td>
<td style="text-align:right;">

400679.56

</td>
<td style="text-align:right;">

1.79

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

203200.12

</td>
<td style="text-align:right;">

12999.00

</td>
<td style="text-align:right;">

1440.41

</td>
<td style="text-align:right;">

347241.48

</td>
<td style="text-align:right;">

20249.31

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

31600.30

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

89628.01

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

3760.12

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

50423.17

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

372978.71

</td>
<td style="text-align:right;">

155.25

</td>
<td style="text-align:right;">

133.46

</td>
<td style="text-align:right;">

241742.53

</td>
<td style="text-align:right;">

302.73

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

255.15

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

122773.32

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

1.04

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.42

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

1.23

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

15384.83

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

267947.05

</td>
<td style="text-align:right;">

58.05

</td>
<td style="text-align:right;">

46.37

</td>
<td style="text-align:right;">

81016.26

</td>
<td style="text-align:right;">

147.91

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
</tbody>
</table>
</details>

  

<details>
<summary>

Table <a href="#tab:tab33">3.3</a>: Simulation with deterministic trend:
bias for the other parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab33"></span>Table 3.3:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1731.64

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

121856.20

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

605.84

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.23

</td>
<td style="text-align:right;">

51326.28

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

134.75

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.23

</td>
<td style="text-align:right;">

9902.13

</td>
<td style="text-align:right;">

9.65

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

202.29

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

95551.01

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

3157.28

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

-2.04

</td>
<td style="text-align:right;">

76496.16

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

60318.22

</td>
<td style="text-align:right;">

-7.96

</td>
<td style="text-align:right;">

-449.94

</td>
<td style="text-align:right;">

473621.50

</td>
<td style="text-align:right;">

-5.15

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

57.24

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

21594.18

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.94

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

89.69

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

28391.57

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

47647.15

</td>
<td style="text-align:right;">

5.04

</td>
<td style="text-align:right;">

-86.16

</td>
<td style="text-align:right;">

251078.43

</td>
<td style="text-align:right;">

10.62

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

2.88

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

30446.28

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-1.11

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

98.42

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

16320.71

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

23952.46

</td>
<td style="text-align:right;">

2.07

</td>
<td style="text-align:right;">

-28.17

</td>
<td style="text-align:right;">

147953.18

</td>
<td style="text-align:right;">

7.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

929.29

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

81725.31

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

125.33

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

13275.54

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

7.05

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

7200.27

</td>
<td style="text-align:right;">

50.87

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

23.54

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

82998.33

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

15778.35

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

-4.86

</td>
<td style="text-align:right;">

203112.60

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

59291.76

</td>
<td style="text-align:right;">

781.86

</td>
<td style="text-align:right;">

-2551.27

</td>
<td style="text-align:right;">

234074.55

</td>
<td style="text-align:right;">

101.73

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

11895.43

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

100.51

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

10154.53

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

85091.05

</td>
<td style="text-align:right;">

26.94

</td>
<td style="text-align:right;">

-500.98

</td>
<td style="text-align:right;">

95361.60

</td>
<td style="text-align:right;">

64.93

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

34215.64

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.88

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

9.86

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

5365.54

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

56056.97

</td>
<td style="text-align:right;">

15.70

</td>
<td style="text-align:right;">

-1.85

</td>
<td style="text-align:right;">

43503.05

</td>
<td style="text-align:right;">

43.29

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

288.90

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

50052.72

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

64.08

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

6353.83

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

3.23

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

2470.59

</td>
<td style="text-align:right;">

132.11

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.17

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

13.21

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

60526.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

16714.28

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

197376.29

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

48580.05

</td>
<td style="text-align:right;">

1208.40

</td>
<td style="text-align:right;">

133.01

</td>
<td style="text-align:right;">

200131.85

</td>
<td style="text-align:right;">

1355.43

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

999.37

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

13393.82

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

199.56

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

5322.68

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

156361.73

</td>
<td style="text-align:right;">

61.68

</td>
<td style="text-align:right;">

-13.79

</td>
<td style="text-align:right;">

73858.02

</td>
<td style="text-align:right;">

188.87

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

8.13

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

25592.07

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

824.94

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

93907.84

</td>
<td style="text-align:right;">

33.43

</td>
<td style="text-align:right;">

-2.10

</td>
<td style="text-align:right;">

10526.19

</td>
<td style="text-align:right;">

118.62

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
</tbody>
</table>
</details>

  

Table <a href="#tab:tab34">3.4</a> compares the estimates for $x_t$ and
$c_t$ for the different models by regressing the respective Kalman
smoother-based estimates on the true trend and cycle and reporting the
coefficient of determination. Results are again similar to table
<a href="#tab:tab4">1.4</a>.

<details>
<summary>

Table <a href="#tab:tab34">3.4</a>: Simulation with deterministic trend:
Coefficient of determination from regressing true trend and cycle on
their respective estimates from the Kalman smoother.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab34"></span>Table 3.4:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Trend

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Cycle

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.79

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.90

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.79

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.53

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.65

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.72

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.93

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.81

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.61

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.55

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.78

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.74

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.86

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.81

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.69

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.54

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.49

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.82

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.74

</td>
</tr>
</tbody>
</table>
</details>

  

# 4 Deterministic trends with breaks

The fourth simulation modifies the DGP of the previous simulation by
imposing a trend break at $\lfloor T/2 \rfloor$

$$
y_t = d_t + x_t + c_t, \qquad \Delta_+^d x_t = \eta_t, \qquad c_t - b_1 c_{t-1} - b_{2}c_{t-2} = \epsilon_t, \qquad Q = \begin{pmatrix}
            1 & 0 \\
            0 & \nu
        \end{pmatrix},
$$ where $$
    d_t = \begin{cases}
         t & \text{ if } t \leq \lfloor T/2 \rfloor \\
         t - (1 / 2) (t -  \lfloor T/2 \rfloor) & \text{else.}
    \end{cases}
$$ This break is not anticipated by the fractional UC model, i.e. the
model of section <a href="#sec:trend">3</a> is fit to the DGP with a
trend break at $\lfloor T/2 \rfloor$. The impact of the model
misspecification is then studied in what follows.

Table <a href="#tab:tab41">4.1</a> shows the RMSE and the bias for the
estimated integration orders. Adding an (unanticipated) trend break to
the DGP slightly increases the RMSE for $d_0=0.75$, but not for the
higher $d_0$ (where the trend is differenced out during the estimation).
However, it leaves the relative performance of the estimators under
study unaffected: Again, the parametric models outperform the
benchmarks, QML outperforms CSS, and ARMA approximations yield
surprisingly good estimates. Noticeably, the unanticipated trend break
does not induce an additional bias to the estimate of $d_0$, as a
comparison to section <a href="#sec:trend">3</a> shows.

<details>
<summary>

Table <a href="#tab:tab41">4.1</a>: Simulation with a deterministic
trend with unanticipated break: root mean squared errors (RMSE) and bias
for the integration order estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab41"></span>Table 4.1:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

RMSE

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

bias

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.57

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.50

</td>
<td style="text-align:right;">

-0.25

</td>
<td style="text-align:right;">

0.36

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.51

</td>
<td style="text-align:right;">

-0.36

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.51

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.80

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

0.46

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

1.04

</td>
<td style="text-align:right;">

1.06

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

1.73

</td>
<td style="text-align:right;">

1.52

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

-0.99

</td>
<td style="text-align:right;">

-0.96

</td>
<td style="text-align:right;">

-0.50

</td>
<td style="text-align:right;">

-1.73

</td>
<td style="text-align:right;">

-1.50

</td>
<td style="text-align:right;">

-0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

1.60

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

-0.44

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-1.59

</td>
<td style="text-align:right;">

-1.28

</td>
<td style="text-align:right;">

-0.22

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.32

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

1.12

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

-0.40

</td>
<td style="text-align:right;">

-0.25

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-1.38

</td>
<td style="text-align:right;">

-1.09

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.32

</td>
<td style="text-align:right;">

-0.39

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.70

</td>
<td style="text-align:right;">

-0.66

</td>
<td style="text-align:right;">

-0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

1.01

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

1.74

</td>
<td style="text-align:right;">

1.73

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

-0.95

</td>
<td style="text-align:right;">

-0.83

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

-1.74

</td>
<td style="text-align:right;">

-1.73

</td>
<td style="text-align:right;">

-0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.33

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

1.62

</td>
<td style="text-align:right;">

1.58

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

-0.40

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

-0.46

</td>
<td style="text-align:right;">

-1.61

</td>
<td style="text-align:right;">

-1.57

</td>
<td style="text-align:right;">

-0.90

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

1.45

</td>
<td style="text-align:right;">

1.41

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-1.43

</td>
<td style="text-align:right;">

-1.40

</td>
<td style="text-align:right;">

-0.84

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.28

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.39

</td>
<td style="text-align:right;">

-0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

-0.70

</td>
<td style="text-align:right;">

-0.33

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

1.73

</td>
<td style="text-align:right;">

1.75

</td>
<td style="text-align:right;">

1.32

</td>
<td style="text-align:right;">

-0.93

</td>
<td style="text-align:right;">

-0.72

</td>
<td style="text-align:right;">

-0.28

</td>
<td style="text-align:right;">

-1.72

</td>
<td style="text-align:right;">

-1.74

</td>
<td style="text-align:right;">

-1.31

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.31

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

1.53

</td>
<td style="text-align:right;">

1.65

</td>
<td style="text-align:right;">

1.20

</td>
<td style="text-align:right;">

-0.36

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

-1.52

</td>
<td style="text-align:right;">

-1.64

</td>
<td style="text-align:right;">

-1.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

1.34

</td>
<td style="text-align:right;">

1.50

</td>
<td style="text-align:right;">

1.11

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

-1.32

</td>
<td style="text-align:right;">

-1.49

</td>
<td style="text-align:right;">

-1.10

</td>
</tr>
</tbody>
</table>
</details>

  

Tables <a href="#tab:tab42">4.2</a> and <a href="#tab:tab43">4.3</a>
contain RMSE and bias for $\nu_0$ and the autoregressive parameters.
Again, for $d_0 = 0.75$ adding an (unanticipated) break increases the
RMSE for the estimates of $\nu_0$ and the autoregressive parameters. The
latter in addition become biased whenever $\nu_0$ is very small,
i.e. when the contribution of the cycle to the overall variation is
small. All further conclusions are similar to those drawn from tables
<a href="#tab:tab32">3.2</a> and <a href="#tab:tab33">3.3</a>, i.e. the
simulations with a linear deterministic trend without break.

<details>
<summary>

Table <a href="#tab:tab42">4.2</a>: Simulation with a deterministic
trend with unanticipated break: root mean squared errors (RMSE) for the
other parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab42"></span>Table 4.2:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

19320.16

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

18574.75

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

3687.32

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

49485.01

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

2105.51

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

55537.62

</td>
<td style="text-align:right;">

25.61

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

1474.61

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

5450.62

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

17010.63

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

2.04

</td>
<td style="text-align:right;">

102721.87

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

124337.13

</td>
<td style="text-align:right;">

148.82

</td>
<td style="text-align:right;">

692.03

</td>
<td style="text-align:right;">

605695.54

</td>
<td style="text-align:right;">

150.15

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

13.85

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

88.79

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

1.34

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.37

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

2281.87

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

20452.73

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

118796.30

</td>
<td style="text-align:right;">

17.43

</td>
<td style="text-align:right;">

90.41

</td>
<td style="text-align:right;">

427151.90

</td>
<td style="text-align:right;">

27.25

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

115.40

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

1.03

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

1.63

</td>
<td style="text-align:right;">

1.40

</td>
<td style="text-align:right;">

1.18

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.44

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

317.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

29548.88

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

73260.22

</td>
<td style="text-align:right;">

6.88

</td>
<td style="text-align:right;">

28.51

</td>
<td style="text-align:right;">

290252.34

</td>
<td style="text-align:right;">

11.30

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

4851.26

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

2965.11

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

27380.10

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

93.50

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

38163.48

</td>
<td style="text-align:right;">

140.25

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

4318.44

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

203147.71

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

4.83

</td>
<td style="text-align:right;">

59290.38

</td>
<td style="text-align:right;">

2.85

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

170800.87

</td>
<td style="text-align:right;">

18752.62

</td>
<td style="text-align:right;">

3957.99

</td>
<td style="text-align:right;">

376724.15

</td>
<td style="text-align:right;">

3611.40

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

886.58

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

1.14

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.41

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

2569.48

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

13635.70

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

221429.35

</td>
<td style="text-align:right;">

69.13

</td>
<td style="text-align:right;">

528.97

</td>
<td style="text-align:right;">

270567.75

</td>
<td style="text-align:right;">

109.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

1632.23

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

1.22

</td>
<td style="text-align:right;">

1.25

</td>
<td style="text-align:right;">

1.61

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

1.15

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

5012.19

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

1101.81

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

198546.21

</td>
<td style="text-align:right;">

29.19

</td>
<td style="text-align:right;">

177.87

</td>
<td style="text-align:right;">

192039.06

</td>
<td style="text-align:right;">

89.48

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

12199.36

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

7069.35

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

33.96

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

36621.55

</td>
<td style="text-align:right;">

270.92

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.37

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

14568.38

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

383452.74

</td>
<td style="text-align:right;">

2.18

</td>
<td style="text-align:right;">

7.83

</td>
<td style="text-align:right;">

26506.80

</td>
<td style="text-align:right;">

2.49

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

197836.84

</td>
<td style="text-align:right;">

22431.18

</td>
<td style="text-align:right;">

15576.56

</td>
<td style="text-align:right;">

358071.29

</td>
<td style="text-align:right;">

11704.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

607.74

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.41

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

5788.93

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

372176.03

</td>
<td style="text-align:right;">

151.69

</td>
<td style="text-align:right;">

1515.95

</td>
<td style="text-align:right;">

238112.25

</td>
<td style="text-align:right;">

343.88

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

3436.22

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

1.45

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

1.11

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

1.57

</td>
<td style="text-align:right;">

1.32

</td>
<td style="text-align:right;">

1.01

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.45

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

5370.22

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

256502.25

</td>
<td style="text-align:right;">

58.15

</td>
<td style="text-align:right;">

504.05

</td>
<td style="text-align:right;">

104436.51

</td>
<td style="text-align:right;">

147.89

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
</tr>
</tbody>
</table>
</details>

  

<details>
<summary>

Table <a href="#tab:tab43">4.3</a>: Simulation with a deterministic
trend with unanticipated break: bias for the other parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab43"></span>Table 4.3:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

788.57

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

1044.03

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

286.85

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

7608.82

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

180.09

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.23

</td>
<td style="text-align:right;">

13146.25

</td>
<td style="text-align:right;">

10.77

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

82.77

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

366.68

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.07

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

3749.95

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-2.02

</td>
<td style="text-align:right;">

32486.43

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

51245.54

</td>
<td style="text-align:right;">

-7.51

</td>
<td style="text-align:right;">

-455.93

</td>
<td style="text-align:right;">

451006.39

</td>
<td style="text-align:right;">

-4.37

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

2.78

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-1.18

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

139.67

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

2400.94

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

48561.11

</td>
<td style="text-align:right;">

5.20

</td>
<td style="text-align:right;">

-86.69

</td>
<td style="text-align:right;">

249847.17

</td>
<td style="text-align:right;">

11.19

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

9.57

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.52

</td>
<td style="text-align:right;">

-0.41

</td>
<td style="text-align:right;">

-1.48

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

13.74

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

2332.44

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.51

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

24610.87

</td>
<td style="text-align:right;">

1.96

</td>
<td style="text-align:right;">

-28.35

</td>
<td style="text-align:right;">

132914.63

</td>
<td style="text-align:right;">

7.33

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

334.32

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.47

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

515.96

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

1215.82

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

13.13

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

5905.70

</td>
<td style="text-align:right;">

50.95

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

319.73

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.44

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

62977.59

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-4.83

</td>
<td style="text-align:right;">

6000.67

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

51628.99

</td>
<td style="text-align:right;">

1304.45

</td>
<td style="text-align:right;">

-2562.97

</td>
<td style="text-align:right;">

233449.89

</td>
<td style="text-align:right;">

222.18

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

48.27

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.25

</td>
<td style="text-align:right;">

-1.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

363.21

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

460.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

81697.18

</td>
<td style="text-align:right;">

27.52

</td>
<td style="text-align:right;">

-503.55

</td>
<td style="text-align:right;">

97078.06

</td>
<td style="text-align:right;">

66.54

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

164.29

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.69

</td>
<td style="text-align:right;">

-0.76

</td>
<td style="text-align:right;">

-1.53

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

333.29

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

58.59

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

66989.99

</td>
<td style="text-align:right;">

16.16

</td>
<td style="text-align:right;">

-177.45

</td>
<td style="text-align:right;">

44573.55

</td>
<td style="text-align:right;">

48.95

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1418.85

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

771.49

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

4.46

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

3147.25

</td>
<td style="text-align:right;">

137.72

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

2038.29

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.54

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

162533.15

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

-7.81

</td>
<td style="text-align:right;">

1166.06

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

46287.79

</td>
<td style="text-align:right;">

1868.21

</td>
<td style="text-align:right;">

-15174.28

</td>
<td style="text-align:right;">

203946.38

</td>
<td style="text-align:right;">

1135.36

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

34.40

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.88

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

756.26

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

157495.64

</td>
<td style="text-align:right;">

60.79

</td>
<td style="text-align:right;">

-1466.74

</td>
<td style="text-align:right;">

69834.96

</td>
<td style="text-align:right;">

195.40

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

550.47

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

-0.85

</td>
<td style="text-align:right;">

-1.51

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

223.23

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

88772.05

</td>
<td style="text-align:right;">

33.50

</td>
<td style="text-align:right;">

-501.54

</td>
<td style="text-align:right;">

14093.70

</td>
<td style="text-align:right;">

118.73

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.05

</td>
</tr>
</tbody>
</table>
</details>

  

Table <a href="#tab:tab44">4.4</a> compares the estimates for $x_t$ and
$c_t$ for the different models by regressing the respective Kalman
smoother-based estimates on the true trend and cycle and reporting the
coefficient of determination. Again, allowing for an unanticipated trend
break does not show deviations from the relative ordering in table
<a href="#tab:tab34">3.4</a>; however for low $d_0$ the coefficients of
determination for the trend are somewhat smaller

<details>
<summary>

Table <a href="#tab:tab44">4.4</a>: Simulation with a deterministic
trend with unanticipated break: Coefficient of determination from
regressing true trend and cycle on their respective estimates from the
Kalman smoother.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab44"></span>Table 4.4:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Trend

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Cycle

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.77

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.88

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.82

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.77

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.37

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.52

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.60

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.71

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.80

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.57

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.54

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.37

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.77

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.86

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.93

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.67

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.54

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.48

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.82

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
</tbody>
</table>
</details>

  

# 5 Outliers

The fifth simulation considers a data-generating mechanism with
unanticipated outlier

$$
y_t = u_t + x_t + c_t, \qquad \Delta_+^d x_t = \eta_t, \qquad c_t - b_1 c_{t-1} - b_{2}c_{t-2} = \epsilon_t, \qquad Q = \begin{pmatrix}
            1 & 0 \\
            0 & \nu
        \end{pmatrix},
$$ where $$
    u_t = \begin{cases}
        - 10 \sqrt{\nu} & \text{ if } t = t^*, \\
        0 & \text{else,}
    \end{cases}
$$ and $t^*$ is drawn from
$\{\lfloor 0.9 T \rfloor, \lfloor 0.9 T \rfloor + 1, ..., T-1, T\}$ to
mimic a situation with an outlier at the end of the sample, a situation
where detecting an outlier is particularly difficult. The magnitude of
the outlier is $10$ times the standard deviation of the short-run
innovations. The fUC model does not anticipate the outlier, i.e. the
model of section <a href="#sec:baseline">1</a> is fit to the above DGP.

Table <a href="#tab:tab51">5.1</a> details RMSE and bias for the
estimated integration orders. Naturally, introducing the outlier
increases the RMSE as compared to subsection
<a href="#sec:baseline">1</a>, particularly for large $\nu_0$. Moreover,
$\hat d$ is downward-biased, and the higher the $\nu_0$ (and thus the
magnitude of the bias), the larger the bias. For $\nu_0$ fixed, the bias
decreases as $n$ increases. The relative performance of the different
estimators is in line with the results of table
<a href="#tab:tab1">1.1</a>, the simulation without outlier.

<details>
<summary>

Table <a href="#tab:tab51">5.1</a>: Simulation with unanticipated
outlier: root mean squared errors (RMSE) and bias for the integration
order estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab51"></span>Table 5.1:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

RMSE

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

bias

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$d$<sub>QML</sub>

</th>
<th style="text-align:right;">

$d$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$d^{.5}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.6}$<sub>EW</sub>

</th>
<th style="text-align:right;">

$d^{.7}$<sub>EW</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

-0.38

</td>
<td style="text-align:right;">

0.32

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.65

</td>
<td style="text-align:right;">

-0.44

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.51

</td>
<td style="text-align:right;">

-0.41

</td>
<td style="text-align:right;">

-0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

-0.58

</td>
<td style="text-align:right;">

-0.34

</td>
<td style="text-align:right;">

0.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

-0.83

</td>
<td style="text-align:right;">

-0.60

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

1.07

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

1.58

</td>
<td style="text-align:right;">

1.38

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

-1.02

</td>
<td style="text-align:right;">

-0.37

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-1.57

</td>
<td style="text-align:right;">

-1.37

</td>
<td style="text-align:right;">

-0.58

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

0.09

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

1.14

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

-0.83

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-1.29

</td>
<td style="text-align:right;">

-1.13

</td>
<td style="text-align:right;">

-0.56

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.04

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

0.05

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

1.15

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

-0.73

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-1.12

</td>
<td style="text-align:right;">

-0.98

</td>
<td style="text-align:right;">

-0.52

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

-0.60

</td>
<td style="text-align:right;">

-0.56

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.57

</td>
<td style="text-align:right;">

-0.55

</td>
<td style="text-align:right;">

-0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.23

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.35

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

-0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.55

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

-0.86

</td>
<td style="text-align:right;">

-0.82

</td>
<td style="text-align:right;">

-0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

1.07

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

1.60

</td>
<td style="text-align:right;">

1.58

</td>
<td style="text-align:right;">

1.02

</td>
<td style="text-align:right;">

-1.05

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.37

</td>
<td style="text-align:right;">

-1.60

</td>
<td style="text-align:right;">

-1.58

</td>
<td style="text-align:right;">

-1.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

0.00

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.53

</td>
<td style="text-align:right;">

-0.50

</td>
<td style="text-align:right;">

-0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

1.35

</td>
<td style="text-align:right;">

1.35

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

-0.88

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-1.34

</td>
<td style="text-align:right;">

-1.34

</td>
<td style="text-align:right;">

-0.93

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.34

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

1.21

</td>
<td style="text-align:right;">

1.22

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

-0.78

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-1.19

</td>
<td style="text-align:right;">

-1.20

</td>
<td style="text-align:right;">

-0.86

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

-0.57

</td>
<td style="text-align:right;">

-0.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.41

</td>
<td style="text-align:right;">

-0.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.62

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

-0.51

</td>
<td style="text-align:right;">

-0.61

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.77

</td>
<td style="text-align:right;">

-0.88

</td>
<td style="text-align:right;">

-0.49

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

1.07

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

1.55

</td>
<td style="text-align:right;">

1.63

</td>
<td style="text-align:right;">

1.28

</td>
<td style="text-align:right;">

-1.05

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.28

</td>
<td style="text-align:right;">

-1.54

</td>
<td style="text-align:right;">

-1.63

</td>
<td style="text-align:right;">

-1.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

-0.06

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

-0.44

</td>
<td style="text-align:right;">

-0.59

</td>
<td style="text-align:right;">

-0.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

1.28

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

1.15

</td>
<td style="text-align:right;">

-0.89

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-1.27

</td>
<td style="text-align:right;">

-1.41

</td>
<td style="text-align:right;">

-1.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

-0.42

</td>
<td style="text-align:right;">

-0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

1.14

</td>
<td style="text-align:right;">

1.29

</td>
<td style="text-align:right;">

1.06

</td>
<td style="text-align:right;">

-0.75

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-1.12

</td>
<td style="text-align:right;">

-1.28

</td>
<td style="text-align:right;">

-1.05

</td>
</tr>
</tbody>
</table>
</details>

  

Tables <a href="#tab:tab52">5.2</a> and <a href="#tab:tab53">5.3</a>
contain RMSE and bias for $\nu_0$ and the autoregressive parameters.
While the outlier increases the RMSE significantly, the relative
performance of all estimators is similar to the results in tables
<a href="#tab:tab2">1.2</a> and <a href="#tab:tab3">1.3</a>. Naturally,
the outlier upward-biases the estimate for $\nu_0$ for the QML
estimator, and this bias decreases in $n$.

<details>
<summary>

Table <a href="#tab:tab52">5.2</a>: Simulation with unanticipated
outlier: root mean squared errors (RMSE) for the other parameter
estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab52"></span>Table 5.2:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

45494.30

</td>
<td style="text-align:right;">

2.09

</td>
<td style="text-align:right;">

2.78

</td>
<td style="text-align:right;">

67748.39

</td>
<td style="text-align:right;">

3.07

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.53

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

6411.45

</td>
<td style="text-align:right;">

2.58

</td>
<td style="text-align:right;">

2.74

</td>
<td style="text-align:right;">

33573.41

</td>
<td style="text-align:right;">

2.94

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.50

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

6.93

</td>
<td style="text-align:right;">

3.72

</td>
<td style="text-align:right;">

3.32

</td>
<td style="text-align:right;">

2792.80

</td>
<td style="text-align:right;">

284.81

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

31816.67

</td>
<td style="text-align:right;">

1.25

</td>
<td style="text-align:right;">

1.84

</td>
<td style="text-align:right;">

71078.50

</td>
<td style="text-align:right;">

2.04

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.55

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

44506.61

</td>
<td style="text-align:right;">

36.10

</td>
<td style="text-align:right;">

11.18

</td>
<td style="text-align:right;">

56384.96

</td>
<td style="text-align:right;">

11.27

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.50

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

73356.80

</td>
<td style="text-align:right;">

3402.89

</td>
<td style="text-align:right;">

3078.36

</td>
<td style="text-align:right;">

107142.74

</td>
<td style="text-align:right;">

3387.54

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.53

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

31639.92

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

16250.14

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

1.01

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.47

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

4559.32

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

42570.77

</td>
<td style="text-align:right;">

1.09

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.50

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

66493.91

</td>
<td style="text-align:right;">

325.70

</td>
<td style="text-align:right;">

318.58

</td>
<td style="text-align:right;">

44160.56

</td>
<td style="text-align:right;">

410.82

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.56

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

44707.25

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

33436.54

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

1.46

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.45

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

542.65

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

3969.49

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.43

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.42

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

14949.99

</td>
<td style="text-align:right;">

106.86

</td>
<td style="text-align:right;">

106.39

</td>
<td style="text-align:right;">

23775.39

</td>
<td style="text-align:right;">

359.46

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.49

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

55841.25

</td>
<td style="text-align:right;">

1.04

</td>
<td style="text-align:right;">

1.65

</td>
<td style="text-align:right;">

40248.60

</td>
<td style="text-align:right;">

1.79

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.33

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

172.85

</td>
<td style="text-align:right;">

1.38

</td>
<td style="text-align:right;">

1.39

</td>
<td style="text-align:right;">

8291.25

</td>
<td style="text-align:right;">

1.56

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

29.96

</td>
<td style="text-align:right;">

2.13

</td>
<td style="text-align:right;">

1.97

</td>
<td style="text-align:right;">

1456.45

</td>
<td style="text-align:right;">

384.78

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.52

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

71375.22

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

1.49

</td>
<td style="text-align:right;">

24456.38

</td>
<td style="text-align:right;">

1.67

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.33

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

68709.39

</td>
<td style="text-align:right;">

28.04

</td>
<td style="text-align:right;">

12.81

</td>
<td style="text-align:right;">

68162.34

</td>
<td style="text-align:right;">

13.22

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.31

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

78752.93

</td>
<td style="text-align:right;">

10764.37

</td>
<td style="text-align:right;">

10332.75

</td>
<td style="text-align:right;">

54054.61

</td>
<td style="text-align:right;">

42329.53

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.33

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.37

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

31669.12

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

7913.40

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.54

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

79.42

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

4251.68

</td>
<td style="text-align:right;">

1.15

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.27

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

48379.51

</td>
<td style="text-align:right;">

1078.09

</td>
<td style="text-align:right;">

1005.70

</td>
<td style="text-align:right;">

40302.72

</td>
<td style="text-align:right;">

1213.19

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.39

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

513.07

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

21636.66

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

1.15

</td>
<td style="text-align:right;">

1.27

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.42

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

39259.69

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

5686.00

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

31561.98

</td>
<td style="text-align:right;">

373.78

</td>
<td style="text-align:right;">

340.56

</td>
<td style="text-align:right;">

195.06

</td>
<td style="text-align:right;">

476.64

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

55802.81

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

16649.37

</td>
<td style="text-align:right;">

1.24

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

1.27

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.18

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

21.57

</td>
<td style="text-align:right;">

1.43

</td>
<td style="text-align:right;">

1.49

</td>
<td style="text-align:right;">

10518.35

</td>
<td style="text-align:right;">

452.16

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.12

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

46082.68

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

1.25

</td>
<td style="text-align:right;">

32412.05

</td>
<td style="text-align:right;">

1.41

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

45106.42

</td>
<td style="text-align:right;">

20.54

</td>
<td style="text-align:right;">

13.52

</td>
<td style="text-align:right;">

73840.84

</td>
<td style="text-align:right;">

13.64

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.23

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

81302.97

</td>
<td style="text-align:right;">

19798.31

</td>
<td style="text-align:right;">

20094.30

</td>
<td style="text-align:right;">

81894.39

</td>
<td style="text-align:right;">

60952.46

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

27.26

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

28647.29

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.52

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

1.04

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

1.40

</td>
<td style="text-align:right;">

1.15

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.19

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

34121.62

</td>
<td style="text-align:right;">

3093.92

</td>
<td style="text-align:right;">

1946.71

</td>
<td style="text-align:right;">

2985.87

</td>
<td style="text-align:right;">

3011.49

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.17

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.30

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

31609.95

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

30859.61

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

1.09

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.47

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

14684.31

</td>
<td style="text-align:right;">

720.51

</td>
<td style="text-align:right;">

653.19

</td>
<td style="text-align:right;">

535.75

</td>
<td style="text-align:right;">

866.76

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

0.31

</td>
</tr>
</tbody>
</table>
</details>

  

<details>
<summary>

Table <a href="#tab:tab53">5.3</a>: Simulation with unanticipated
outlier: bias for the other parameter estimates.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab53"></span>Table 5.3:

</caption>
<thead>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>QML</sub>

</th>
<th style="text-align:right;">

$\nu$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$\nu^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_1$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_1^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$b_2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$b_2^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

2439.83

</td>
<td style="text-align:right;">

1.61

</td>
<td style="text-align:right;">

2.59

</td>
<td style="text-align:right;">

9804.01

</td>
<td style="text-align:right;">

2.94

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.52

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.49

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

283.16

</td>
<td style="text-align:right;">

2.32

</td>
<td style="text-align:right;">

2.58

</td>
<td style="text-align:right;">

2198.46

</td>
<td style="text-align:right;">

2.81

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.36

</td>
<td style="text-align:right;">

-0.46

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.46

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

3.32

</td>
<td style="text-align:right;">

2.91

</td>
<td style="text-align:right;">

167.80

</td>
<td style="text-align:right;">

169.32

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.41

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

1213.98

</td>
<td style="text-align:right;">

1.03

</td>
<td style="text-align:right;">

1.66

</td>
<td style="text-align:right;">

9128.29

</td>
<td style="text-align:right;">

1.94

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.51

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

3411.81

</td>
<td style="text-align:right;">

23.89

</td>
<td style="text-align:right;">

10.77

</td>
<td style="text-align:right;">

6065.93

</td>
<td style="text-align:right;">

10.88

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.39

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.42

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.47

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

8801.92

</td>
<td style="text-align:right;">

3136.03

</td>
<td style="text-align:right;">

2976.89

</td>
<td style="text-align:right;">

14955.45

</td>
<td style="text-align:right;">

3246.89

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.49

</td>
<td style="text-align:right;">

-0.48

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.50

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

1045.41

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

1594.75

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.49

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.14

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

148.67

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

3154.55

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

-0.32

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.34

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.43

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

5451.58

</td>
<td style="text-align:right;">

313.76

</td>
<td style="text-align:right;">

306.81

</td>
<td style="text-align:right;">

2519.40

</td>
<td style="text-align:right;">

392.80

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.47

</td>
<td style="text-align:right;">

-0.47

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

-0.53

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

2042.67

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

3152.09

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.24

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

30.74

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

240.20

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.15

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.20

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

779.04

</td>
<td style="text-align:right;">

103.07

</td>
<td style="text-align:right;">

102.69

</td>
<td style="text-align:right;">

1055.89

</td>
<td style="text-align:right;">

289.97

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.46

</td>
<td style="text-align:right;">

-0.46

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.45

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

3647.44

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

1.47

</td>
<td style="text-align:right;">

3052.25

</td>
<td style="text-align:right;">

1.72

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.31

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

7.75

</td>
<td style="text-align:right;">

1.28

</td>
<td style="text-align:right;">

1.29

</td>
<td style="text-align:right;">

371.93

</td>
<td style="text-align:right;">

1.49

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.26

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

2.37

</td>
<td style="text-align:right;">

2.01

</td>
<td style="text-align:right;">

1.83

</td>
<td style="text-align:right;">

68.00

</td>
<td style="text-align:right;">

263.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

-0.33

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.17

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.45

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

5620.90

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

1.31

</td>
<td style="text-align:right;">

2067.42

</td>
<td style="text-align:right;">

1.61

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.27

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.31

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

5135.01

</td>
<td style="text-align:right;">

19.71

</td>
<td style="text-align:right;">

12.48

</td>
<td style="text-align:right;">

5867.70

</td>
<td style="text-align:right;">

12.73

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.26

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.30

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

1859.20

</td>
<td style="text-align:right;">

10078.36

</td>
<td style="text-align:right;">

10072.01

</td>
<td style="text-align:right;">

991.00

</td>
<td style="text-align:right;">

13715.65

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.31

</td>
<td style="text-align:right;">

-0.29

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

1060.46

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

361.56

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.07

</td>
<td style="text-align:right;">

-0.30

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

3.15

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

150.14

</td>
<td style="text-align:right;">

1.08

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.24

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

2447.35

</td>
<td style="text-align:right;">

1055.31

</td>
<td style="text-align:right;">

982.46

</td>
<td style="text-align:right;">

1708.28

</td>
<td style="text-align:right;">

1175.11

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.37

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

18.75

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

1067.47

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

-0.06

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.08

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

1790.92

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

214.86

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.11

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.12

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.17

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

1099.96

</td>
<td style="text-align:right;">

365.87

</td>
<td style="text-align:right;">

332.69

</td>
<td style="text-align:right;">

-194.93

</td>
<td style="text-align:right;">

464.21

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.30

</td>
<td style="text-align:right;">

-0.04

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.30

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.38

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

3625.41

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

1040.74

</td>
<td style="text-align:right;">

1.19

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.09

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

-0.52

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

-0.52

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

1.57

</td>
<td style="text-align:right;">

1.37

</td>
<td style="text-align:right;">

1.42

</td>
<td style="text-align:right;">

612.09

</td>
<td style="text-align:right;">

335.52

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

-0.14

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

-0.19

</td>
<td style="text-align:right;">

-0.08

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.32

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.01

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

2651.06

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

1.04

</td>
<td style="text-align:right;">

2057.19

</td>
<td style="text-align:right;">

1.36

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

2334.29

</td>
<td style="text-align:right;">

16.34

</td>
<td style="text-align:right;">

13.18

</td>
<td style="text-align:right;">

7086.22

</td>
<td style="text-align:right;">

13.37

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.20

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.21

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

-7395.95

</td>
<td style="text-align:right;">

19409.84

</td>
<td style="text-align:right;">

19682.32

</td>
<td style="text-align:right;">

-4724.74

</td>
<td style="text-align:right;">

27027.24

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.21

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.25

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

1.20

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

1164.24

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.05

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

0.06

</td>
<td style="text-align:right;">

-0.34

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

-0.62

</td>
<td style="text-align:right;">

1.10

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.13

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.14

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.17

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

547.45

</td>
<td style="text-align:right;">

2219.61

</td>
<td style="text-align:right;">

1913.25

</td>
<td style="text-align:right;">

-1524.33

</td>
<td style="text-align:right;">

2562.83

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.23

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.23

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

-0.29

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

1002.06

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

1724.37

</td>
<td style="text-align:right;">

0.03

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.47

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

-0.05

</td>
<td style="text-align:right;">

0.00

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

-0.16

</td>
<td style="text-align:right;">

-0.15

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

-0.10

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.19

</td>
<td style="text-align:right;">

-0.12

</td>
<td style="text-align:right;">

0.25

</td>
<td style="text-align:right;">

0.01

</td>
<td style="text-align:right;">

-0.09

</td>
<td style="text-align:right;">

-0.07

</td>
<td style="text-align:right;">

-0.01

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

-0.02

</td>
<td style="text-align:right;">

0.10

</td>
<td style="text-align:right;">

0.08

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

1798.93

</td>
<td style="text-align:right;">

709.36

</td>
<td style="text-align:right;">

641.92

</td>
<td style="text-align:right;">

-535.74

</td>
<td style="text-align:right;">

839.20

</td>
<td style="text-align:right;">

0.02

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.22

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

-0.03

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.04

</td>
<td style="text-align:right;">

-0.30

</td>
</tr>
</tbody>
</table>
</details>

  

Table <a href="#tab:tab54">5.4</a> again compares the estimates for
$x_t$ and $c_t$ for the different models by means of the coefficient of
determination. The unanticipated outlier reduces the coefficient of
determination both for trend and cycle, which is natural as the outlier
is not captured by the model. All further conclusions from table
<a href="#tab:tab54">5.4</a> are similar to those from table
<a href="#tab:tab4">1.4</a>.

<details>
<summary>

Table <a href="#tab:tab54">5.4</a>: Simulation with unanticipated
outlier: Coefficient of determination from regressing true trend and
cycle on their respective estimates from the Kalman smoother.

</summary>
<table class=" lightable-paper lightable-hover" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>

<span id="tab:tab54"></span>Table 5.4:

</caption>
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="4">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Trend

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="5">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Cycle

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">

$n$

</th>
<th style="text-align:left;">

$r$

</th>
<th style="text-align:left;">

$d$<sub>0</sub>

</th>
<th style="text-align:left;">

$\nu$<sub>0</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>CSS</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>QML</sub>

</th>
<th style="text-align:right;">

$R^2$<sub>ARMA</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>CSS</sub>

</th>
<th style="text-align:right;">

${R^2}^{I(1)}$<sub>QML</sub>

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.55

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.76

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.89

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.65

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.50

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.80

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

3.82

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.71

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1053.88

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.41

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.37

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.10

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.07

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.18

</td>
<td style="text-align:right;">

0.16

</td>
<td style="text-align:right;">

0.11

</td>
<td style="text-align:right;">

0.11

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.39

</td>
<td style="text-align:right;">

0.38

</td>
<td style="text-align:right;">

0.40

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

105.39

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.70

</td>
<td style="text-align:right;">

0.56

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.02

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.36

</td>
<td style="text-align:right;">

0.35

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.28

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.13

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.59

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.60

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

35.13

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.73

</td>
</tr>
<tr>
<td style="text-align:left;">

200

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.61

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.83

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.73

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.90

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.93

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.93

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.60

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.77

</td>
<td style="text-align:right;">

0.84

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

7.59

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.77

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

5871.84

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.98

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.44

</td>
<td style="text-align:right;">

0.26

</td>
<td style="text-align:right;">

0.03

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.09

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.21

</td>
<td style="text-align:right;">

0.22

</td>
<td style="text-align:right;">

0.20

</td>
<td style="text-align:right;">

0.13

</td>
<td style="text-align:right;">

0.13

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.76

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.58

</td>
<td style="text-align:right;">

0.57

</td>
<td style="text-align:right;">

0.60

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

587.18

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.60

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.03

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.46

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.48

</td>
<td style="text-align:right;">

0.31

</td>
<td style="text-align:right;">

0.38

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.25

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.74

</td>
<td style="text-align:right;">

0.75

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

195.73

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.78

</td>
</tr>
<tr>
<td style="text-align:left;">

300

</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.64

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.63

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.86

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.72

</td>
<td style="text-align:right;">

0.75

</td>
<td style="text-align:right;">

0.76

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.93

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:right;">

0.83

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.82

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.95

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.66

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.65

</td>
<td style="text-align:right;">

0.69

</td>
<td style="text-align:right;">

0.85

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.85

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

11.37

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.78

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

16098.90

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

1.00

</td>
<td style="text-align:right;">

0.99

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.53

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.51

</td>
<td style="text-align:right;">

0.24

</td>
<td style="text-align:right;">

0.02

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

10

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.11

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.27

</td>
<td style="text-align:right;">

0.29

</td>
<td style="text-align:right;">

0.28

</td>
<td style="text-align:right;">

0.15

</td>
<td style="text-align:right;">

0.16

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

1.14

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.67

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.68

</td>
<td style="text-align:right;">

0.68

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

1609.89

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.97

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.60

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

30

</td>
<td style="text-align:left;">

0.75

</td>
<td style="text-align:left;">

0.04

</td>
<td style="text-align:right;">

0.87

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.88

</td>
<td style="text-align:right;">

0.84

</td>
<td style="text-align:right;">

0.86

</td>
<td style="text-align:right;">

0.54

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.56

</td>
<td style="text-align:right;">

0.40

</td>
<td style="text-align:right;">

0.48

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.00

</td>
<td style="text-align:left;">

0.38

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.90

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.81

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.79

</td>
<td style="text-align:right;">

0.80

</td>
<td style="text-align:right;">

0.80

</td>
</tr>
<tr>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">

1.75

</td>
<td style="text-align:left;">

536.63

</td>
<td style="text-align:right;">

0.93

</td>
<td style="text-align:right;">

0.95

</td>
<td style="text-align:right;">

0.96

</td>
<td style="text-align:right;">

0.91

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.92

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.94

</td>
<td style="text-align:right;">

0.89

</td>
<td style="text-align:right;">

0.78

</td>
</tr>
</tbody>
</table>
</details>

  

# 6 References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-HarWei2018a" class="csl-entry">

Hartl, Tobias, and Roland Jucknewitz. 2022. “Approximate State Space
Modelling of Unobserved Fractional Components.” *Econometric Reviews* 41
(1): 75–98.

</div>

<div id="ref-ShiPhi2005" class="csl-entry">

Shimotsu, Katsumi, and Peter C. B. Phillips. 2005. “Exact Local Whittle
Estimation of Fractional Integration.” *The Annals of Statistics* 33
(4): 1890–1933.

</div>

</div>
