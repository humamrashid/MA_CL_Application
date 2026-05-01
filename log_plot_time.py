#!/usr/bin/python3

from matplotlib import pyplot as plt

xs = [ 2**k for k in range(0, 11) ]
n_elapsed = [4.38e-07, 7.7e-07, 2.579e-06, 1.3577e-05, 0.000104062, \
        0.000920927, 0.00504276, 0.0282846, 0.222668, 1.88236, 19.4296]
s_elapsed = [6.71e-07, 1.3809e-05, 0.000119518, 0.000928153, \
        0.00506581, 0.0244157, 0.158543, 1.06294, 7.43906, 52.1426, 365.632]

fig, ax1 = plt.subplots(1)
strassen_lines, = ax1.loglog(xs, s_elapsed, color='green', lw=2, basex=2, basey=2)
naive_lines, = ax1.loglog(xs, n_elapsed, color='blue', lw=2, basex=2, basey=2)

fig.suptitle('Matrix Multiplication Plot (Time)')
ax1.set_title('Computation Time')
ax1.set_xlabel('Matrix size (n x n)')
ax1.set_ylabel('Elapsed time (in seconds)')
strassen_lines.set_label('Strassen Algorithm')
naive_lines.set_label('Naive Algorithm')
ax1.legend()
plt.show()

# EOF.
