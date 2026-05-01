#!/usr/bin/python3

from matplotlib import pyplot as plt

xs = [ 2**k for k in range(0, 11) ]
n_muls = [1, 8, 64, 512, 4096, 32768, 262144, \
        2097152, 16777216, 134217728, 1073741824]
n_adds = [0, 4, 48, 448, 3840, 31744, 258048, \
        2080768, 16711680, 133955584, 1072693248]
n_total = [1, 12, 112, 960, 7936, 64512, 520192, \
        4177920, 33488896, 268173312, 2146435072]
s_muls = [1, 7, 49, 343, 2401, 16807, \
        117649, 823543, 5764801, 40353607, 282475249]
s_adds = [0, 18, 198, 1674, 12870, 94698, \
        681318, 4842954, 34195590, 240548778, 1688560038]

s_total = [1, 25, 247, 2017, 15271, 111505, \
        798967, 5666497, 39960391, 280902385, 1971035287]

strassen_lines = {}
naive_lines = {}
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
#fig, ax3 = plt.subplots(1)
strassen_lines[0], = ax1.loglog(xs, s_muls, color='green', lw=2, basex=2, \
        basey=2)
naive_lines[0], = ax1.loglog(xs, n_muls, color='blue', lw=2, basex=2, basey=2)
strassen_lines[1], = ax2.loglog(xs, s_adds, color='green', lw=2, basex=2, \
        basey=2)
naive_lines[1], = ax2.loglog(xs, n_adds, color='blue', lw=2, basex=2, basey=2)
#strassen_lines[0], = ax3.loglog(xs, s_total, color='green', lw=2, basex=2, \
        #basey=2)
#naive_lines[0], = ax3.loglog(xs, n_total, color='blue', lw=2, basex=2, \
        #basey=2)

fig.suptitle('Matrix Multiplication Plot (Operations)')
ax1.set_title('Multiplication Operations')
ax2.set_title('Addition Operations')
ax2.set_xlabel('Matrix size (n x n)')
ax1.set_ylabel('Number of operations')
ax2.set_ylabel('Number of operations')
#ax3.set_title('Total Operations')
#ax3.set_xlabel('Matrix size (n x n)')
#ax3.set_ylabel('Number of operations')

for i in range (0, 2):
#for i in range (0, 1):
    strassen_lines[i].set_label('Strassen Algorithm')
    naive_lines[i].set_label('Naive Algorithm')
ax1.legend()
ax2.legend()
#ax3.legend()
plt.show()

# EOF.
