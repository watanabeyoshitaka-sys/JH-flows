# JH-flows

MATLAB/INTLAB code for the numerical verification (computer-assisted existence proof)
of Jeffery–Hamel flows, accompanying the paper

> Y. Watanabe, *A numerical verification approach for Jeffery–Hamel flows*
> (submitted to *Journal of Fluid Mechanics*).

The code implements the verification algorithm FN-Int (Nakao, Plum & Watanabe 2019)
for the boundary value problem

    g'' + 4g + g^2 = b0  (-alpha < theta < alpha),   g(±alpha) = 0,   (1/2)∫ g dtheta = R,

using piecewise linear finite elements and interval arithmetic.
All verified quantities are computed with rigorous rounding-error control.

## Requirements

- MATLAB R2020a or later (tested with MATLAB R2024a, Update 3)
- INTLAB Version 14 (tested). INTLAB is **not** included in this repository;
  obtain it from its official distribution and run `startintlab` before using this code.

## Files

| File | Role |
| --- | --- |
| `main.m` | Main script: sets parameters, computes an approximate solution, runs the verification |
| `FN_Int.m` | Verification algorithm (Section 3.5.4 of the paper) |
| `Newton_Raphson.m` | Newton–Raphson iteration for the approximate solution (floating point) |
| `initial_guess.m` | Initial guess: an analytic function (`gn = 1`) or a stored solution (`gn = 2`) |
| `matrices_D_L.m`, `matrix_G2.m`, `matrix_G2_intval.m` | Finite element matrices (2.14), (2.18) |
| `local_L2_norm.m` | Section 3.4.3 |
| `L2_norm_of_f.m`, `inf_norm.m` | Sections 3.5.1 and 3.5.2 |
| `generate_local_node.m` | Element and node connectivity |
| `plot_approximate_solution.m` | Plots of u_h and of the velocity profile |
| `initial_values/` | Stored approximate solutions used for the figures of the paper |

## Usage

1. Start MATLAB and INTLAB (`startintlab`).
2. Edit the parameters at the top of `main.m`:
   `a` (alpha), `R`, `eq` (transformation (2.4)/(2.5)/(2.6) as 1/2/3), `N` and `gn`.
3. Run `main`.

`main.m` prints f(0), the radius `beta` (= ‖V*‖_X) and reports
`Verification has been completed!` on success; otherwise it stops with an error.

## Reproducing the figures

`main.m` reads its initial guess from `app.mat` (`gn = 2`) and saves the computed
solution back to `app.mat`. To reproduce a figure, copy the corresponding file
from `initial_values/` to `app.mat`, set the parameters in `main.m` as in the table,
and run `main`:

```matlab
copyfile("initial_values/fig02-R03-a01pi04-eq3-I.mat", "app.mat")
main
```
`app.mat` initially contains the solution of Figure 2, so running `main`
right after cloning reproduces Figure 2.

Since `app.mat` is overwritten after each run, it can be used for a parameter
search (continuation): change `a` or `R` slightly and run `main` again.

| Figure | Type | R | alpha | eq | N | File in `initial_values/` |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | I | 3 | π/4 | 3 | 2000 | `fig02-R03-a01pi04-eq3-I.mat` |
| 3 | II2 | 3 | π/4 | 3 | 2000 | `fig03-R03-a01pi04-eq3-II2.mat` |
| 4 | IV1 | 3 | π/4 | 3 | 5000 | `fig04-R03-a01pi04-eq3-IV1.mat` |
| 5 | III1 | −3 | π/4 | 3 | 2000 | `fig05-RM03-a01pi04-eq3-III1.mat` |
| 6 | II2 | −3 | π/4 | 3 | 10000 | `fig06-RM03-a01pi04-eq3-II2.mat` |
| 7 | IV1 | −3 | π/4 | 3 | 5000 | `fig07-RM03-a01pi04-eq3-IV1.mat` |
| 10 | III3 | 1 | 3π/8 | 2 | 40000 | `fig10-R01-a03dv08Pi-eq2-III3.mat` |
| 13 | II2 | 1 | π/2 | 3 | 5000 | `fig13-R01-a01pi02-eq3-II2.mat` |
| 20 | V2 | 1 | 3π/4 | 2 | 40000 | `fig20-R01-a03dv04Pi-eq2-V2.mat` |
| 21 | II3 | 15 | 5π/8 | 3 | 5000 | `fig21-R15-a05dv08Pi-eq3-II3.mat` |

The velocity profiles in Figures 11, 12, 14, 15 and 17–19 are obtained in the same way
from the files `figNN-*.mat`. File names encode the parameters:
`R03` = 3, `RM03` = −3, `a01pi04` = π/4, `a03dv08Pi` = 3π/8, `eq3` = transformation (2.6).

## License

The code in this repository is released under the MIT License (see `LICENSE`).
INTLAB is distributed separately under its own license.

## Citation

If you use this code, please cite the paper above and this repository (see `CITATION.cff`).
