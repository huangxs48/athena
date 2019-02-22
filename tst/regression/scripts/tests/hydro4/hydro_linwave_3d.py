# Regression test based on Newtonian hydro linear wave convergence problem

# Confirm fourth-order convergence rate for semidiscrete integration with RK4 + PPM (char.
# and primitive) + Laplacian flux correction terms. 3D uniform square grid, no SMR,
# fourth-order accurate approximation to analytic solution used in initialization and in
# error calculations.

# Modules
import scripts.utils.athena as athena
from math import log
import numpy as np

# List of time/integrator and time/xorder combinations to test:
solvers = [('rk4', '4c'), ('ssprk5_4', '4')]
# Matching above list of solver configurations, provide bounds on error metrics:
# for each tested resolution (excluding lowest Nx1=16) and wave_flag.
# Upper bound on RMS-L1 errors:
error_tols = [((5.6e-9, 3.6e-10), (4.05e-9, 2.65e-10)),
              ((5.6e-9, 3.65e-10), (4.05e-9, 2.65e-10))
              ]
# for each wave_flag, lower bound on convergence rate at 64x32x32 asymptotic convergence
# regime. Linear hydro waves stop converging around RMS-L1 error 1e-11 to 1e-12
rate_tols = [(3.95, 3.94),  (3.95, 3.94)]
# this metric is redundant with above error_tols, but it is simpler...

# time/correct_err= time/correct_ic=false (slightly higher errors, but faster convergence)
# error_tols = [((5.7e-9, 3.7e-10), (4.1e-9, 2.7e-10)),
#               ((5.7e-9, 3.7e-10), (4.1e-9, 2.7e-10))
#               ]
# rate_tols = [(3.95, 3.94),  (3.95, 3.94)]

resolution_range = [16, 32, 64]
num_nx1 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
nrows_per_solver = 2*num_nx1 + 2


# Prepare Athena++
def prepare(**kwargs):
    athena.configure(
        nghost=4,  # required for fourth-order configurations
        prob='linear_wave',
        coord='cartesian',
        flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for (torder, xorder) in solvers:
        # L-going sound wave
        for i in resolution_range:
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'time/correct_ic=true', 'time/correct_err=true',
                         'problem/wave_flag=0', 'problem/vflow=0.0',
                         'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2),
                         'mesh/nx3=' + repr(i/2),
                         'meshblock/nx1=' + repr(i), 'meshblock/nx2=' + repr(i/2),
                         'meshblock/nx3=' + repr(i/2),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave3d', arguments)
        # L-going entropy wave
        for i in resolution_range:
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'time/correct_ic=true', 'time/correct_err=true',
                         'problem/wave_flag=3', 'problem/vflow=1.0',
                         'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2),
                         'mesh/nx3=' + repr(i/2),
                         'meshblock/nx1=' + repr(i), 'meshblock/nx2=' + repr(i/2),
                         'meshblock/nx3=' + repr(i/2),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave3d', arguments)
        # L/R-going sound wave, no SMR
        for w in (0, 4):
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'time/correct_ic=true', 'time/correct_err=true',
                         'problem/wave_flag=' + repr(w),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave3d', arguments)


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = np.loadtxt(filename)

    for ((torder, xorder), err_tol, rate_tol) in zip(solvers, error_tols, rate_tols):
        # effectively list.pop() range of rows for this solver configuration
        solver_results = np.array(data[0:nrows_per_solver])
        data = np.delete(data, np.s_[0:nrows_per_solver], 0)

        # Compute error convergence rates with Richardson extrapolation for each wave flag
        # --------------------
        print('{} + {}'.format(torder.upper(), xorder))
        # L-going sound wave
        print("Sound wave error convergence:")
        print("nx1   |   rate   |   RMS-L1")
        rms_errs = solver_results[0:num_nx1, 4]
        nx1_range = solver_results[0:num_nx1, 0]
        for i in range(1, num_nx1):
            rate = log(rms_errs[i-1]/rms_errs[i])/log(nx1_range[i]/nx1_range[i-1])
            print(int(nx1_range[i]), rate, rms_errs[i])
            # old rate calculation from hydro/hydro_linwave.py:
            # print(rms_errs[i]/rms_errs[i-1])
            if (nx1_range[i] == 128 and rate < rate_tol[0]):
                print("L-going sound wave converging at rate {} slower than {}".format(
                    rate, rate_tol[0]))
                return False
            if (rms_errs[i] > err_tol[0][i-1]):
                print("L-going sound wave error {} is larger than tolerance {}".format(
                    rms_errs[i], err_tol[0][i-1]))
                return False

        # L-going entropy wave
        print("Entropy wave error convergence:")
        print("nx1   |   rate   |   RMS-L1")
        rms_errs = solver_results[num_nx1:2*num_nx1, 4]
        nx1_range = solver_results[num_nx1:2*num_nx1, 0]
        for i in range(1, num_nx1):
            rate = log(rms_errs[i-1]/rms_errs[i])/log(nx1_range[i]/nx1_range[i-1])
            print(int(nx1_range[i]), rate, rms_errs[i])
            # old rate calculation from hydro/hydro_linwave.py:
            # print(rms_errs[i]/rms_errs[i-1])
            if (nx1_range[i] == 128 and rate < rate_tol[1]):
                print("L-going entropy wave converging at rate {} slower than {}".format(
                    rate, rate_tol[1]))
                return False
            if (rms_errs[i] > err_tol[1][i-1]):
                print("L-going entropy wave error {} is larger than tolerance {}".format(
                    rms_errs[i], err_tol[1][i-1]))
                return False

        # Check that errors are identical for sound waves in each direction at default
        # 64x32x32 resolution
        if (not np.allclose(solver_results[-2, 4], solver_results[-1, 4],
                            atol=5e-16, rtol=1e-5)):
            print(("L/R-going sound wave errors, {} and {}".format(solver_results[-2, 4],
                                                                   solver_results[-1, 4]),
                   ", have a difference that is not close to round-off"))
            return False

    return True
