reset
TOL = 1e-11
set logscale y
plot TOL w l lc 8 notitle, "../output/CG_convergence.out" u 1:2 w l t "Conjugate Gradient", "../output/SH_convergence.out" u 1:2 w l t "SHAKE", "../output/BSH_convergence.out" u 1:2 w l t "BSHAKE"
