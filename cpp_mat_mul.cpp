/* Matrix multiplication in C++ (targeting C++11 or later).
 *
 * Uses either naive (brute-force) method or Strassen's algorithm depending on
 * user selection.
 *
 * Metrics collected include:
 *  1. Number of multiplications
 *  2. Number of additions
 *  3. Elapsed time.
 */

#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
using namespace std;

// Constants for "large" n and present output file.
#define LARGE_N 10
#define OUT_FILE "mat_mul_out.txt"

// Operation counters.
static unsigned long naive_mul = 0,
                     naive_add = 0,
                     strassen_mul = 0,
                     strassen_add = 0;

/* Print out operations metrics.
 */
static void display_metrics(unsigned long mul, unsigned long add)
{
    cout << "\nNumber of multiplications: " << mul << "\nNumber of additions: "
        << add << "\nTotal number of operations: " << (mul + add) << endl;
    return;
}

/* Print out matrix to the output file OUT_FILE.
 */
static void display_matrix(vector<vector<int>> &matrix, int n, ofstream &ofs)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ofs.width(10);
            ofs << matrix[i][j];
        }
        ofs << endl;
    }
    return;
}

/* Print out matrix to the screen.
 */
static void display_matrix(vector<vector<int>> &matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout.width(10);
            cout << matrix[i][j];
        }
        cout << endl;
    }
    return;
}

/* Add one matrix to another and put result in third matrix.
 */
static void add(
        vector<vector<int>> &A,
        vector<vector<int>> &B,
        vector<vector<int>> &C,
        int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            strassen_add += 1;
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return;
}

/* Subtract one matrix from another and put result in third matrix.
 */
static void subtract(
        vector<vector<int>> &A,
        vector<vector<int>> &B,
        vector<vector<int>> &C,
        int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            strassen_add += 1;
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return;
}

/* Use Strassen's algorithm to compute C = A x B.
 */
static void strassen(
        vector<vector<int>> &A,
        vector<vector<int>> &B,
        vector<vector<int>> &C,
        int size)
{
    /* Multiplication setup for 2 x 2 matrix:
     *
     * Matrix A:
     * [a b]   
     * [c d] 
     *
     * Matrix B:
     * [e g]
     * [f h]
     *
     * p1 = a * (g - h)
     * p2 = (a + b) * h
     * p3 = (c + d) * e
     * p4 = d * (f - e)
     * p5 = (a + d) * (e + h)
     * p6 = (b - d) * (f + h)
     * p7 = (a - c) * (e + g)
     *
     * r = p5 + p4 - p2 + p6
     * s = p1 + p2
     * t = p3 + p4
     * u = p5 + p1 - p3 - p7
     *
     * Matrix C = A x B:
     * [r s]
     * [t u]
     */
    if (size == 1)
    {
        // Base case.
        C[0][0] = A[0][0] * B[0][0];
        strassen_mul += 1;
        return;
    } else
    {
        int new_size = size / 2;
        vector<int> zeros(new_size);

        // Vectors for submatrices and partial results, with elements initialized
        // to zero-valued vectors.
        vector<vector<int>>
            a(new_size, zeros), b(new_size, zeros),
            c(new_size, zeros), d(new_size, zeros),
            e(new_size, zeros), f(new_size, zeros),
            g(new_size, zeros), h(new_size, zeros),
            r(new_size, zeros), s(new_size, zeros),
            t(new_size, zeros), u(new_size, zeros),
            p1(new_size, zeros), p2(new_size, zeros),
            p3(new_size, zeros), p4(new_size, zeros),
            p5(new_size, zeros), p6(new_size, zeros),
            p7(new_size, zeros), A_result(new_size, zeros),
            B_result(new_size, zeros);

        /* Make 8 new submatrices for matrices A and B.
         *
         * A:
         * [a b]
         * [c d]
         *
         * B:
         * [e g]
         * [f h]
         */
        for (int i = 0; i < new_size; i++)
        {
            for (int j = 0; j < new_size; j++)
            {
                a[i][j] = A[i][j];
                b[i][j] = A[i][j + new_size];
                c[i][j] = A[i + new_size][j];
                d[i][j] = A[i + new_size][j + new_size];

                e[i][j] = B[i][j];
                f[i][j] = B[i + new_size][j];
                g[i][j] = B[i][j + new_size];
                h[i][j] = B[i + new_size][j + new_size];
            }
        }

        // Recursively compute p1 through p7:
        
        // p1 = a * (g - h)
        // g - h
        subtract(g, h, B_result, new_size);
        // a * (g - h)
        strassen(a, B_result, p1, new_size);

        // p2 = (a + b) * h
        // a + b
        add(a, b, A_result, new_size);
        // (a + b) * h
        strassen(A_result, h, p2, new_size);

        // p3 = (c + d) * e
        // c + d
        add(c, d, A_result, new_size);
        // (c + d) * e
        strassen(A_result, e, p3, new_size);

        // p4 = d * (f - e)
        // f - e
        subtract(f, e, B_result, new_size);
        // d * (f - e)
        strassen(d, B_result, p4, new_size);

        // p5 = (a + d) * (e + h)
        // a + d
        add(a, d, A_result, new_size);
        // e + h
        add(e, h, B_result, new_size);
        // (a + d) * (e + h)
        strassen(A_result, B_result, p5, new_size);

        // p6 = (b - d) * (f + h)
        // b - d
        subtract(b, d, A_result, new_size);
        // f + h
        add(f, h, B_result, new_size);                
        // (b - d) * (f + h)
        strassen(A_result, B_result, p6, new_size);

        // p7 = (a - c) * (e + g)
        // a - c
        subtract(a, c, A_result, new_size);
        // e + g
        add(e, g, B_result, new_size);               
        // (a - c) * (e + g)
        strassen(A_result, B_result, p7, new_size);

        // Compute final result matrix C.

        // r = p5 + p4 - p2 + p6
        // p5 + p4
        add(p5, p4, A_result, new_size);
        // p5 + p4 + p6
        add(A_result, p6, B_result, new_size);
        // p5 + p4 - p2 + p6
        subtract(B_result, p2, r, new_size);

        // s = p1 + p2
        add(p1, p2, s, new_size);

        // t = p3 + p4
        add(p3, p4, t, new_size);

        // u = p5 + p1 - p3 - p7
        // p5 + p1
        add(p5, p1, A_result, new_size);
        // p5 + p1 - p3
        subtract(A_result, p3, B_result, new_size);
        // p5 + p1 - p3 - p7
        subtract(B_result, p7, u, new_size);

        /* Set submatrices of C.
         *
         * C:
         * [r s]
         * [t u]
         */
        for (int i = 0; i < new_size; i++)
        {
            for (int j = 0; j < new_size; j++)
            {
                C[i][j] = r[i][j];
                C[i][j + new_size] = s[i][j];
                C[i + new_size][j] = t[i][j];
                C[i + new_size][j + new_size] = u[i][j];
            }
        }
    }

    return;
}

// Get next power of 2.
static int nextpowerof2(int k)
{
    return pow(2, int(ceil(log2(k))));
}

/* Compute C = A x B after making sure matrices for any n >= 1 (not only powers
 * of 2) are taken care of and write out the results.
 */
static void compute_mul(
        vector<vector<int>> &A,
        vector<vector<int>> &B,
        int n,
        bool is_naive)
{
    // Get next higher power of 2.
    int p = nextpowerof2(n);
    vector<int> zeros(p);
    vector<vector<int>> Ap(p, zeros), Bp(p, zeros), Cp(p, zeros);

    // Fill in original values of A and B.
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Ap[i][j] = A[i][j];
            Bp[i][j] = B[i][j];
        }
    }

    cout << "\nStarting computation...";
    auto start = chrono::steady_clock::now();
    if (is_naive)
    {
        // Multiply matrices A and B using brute-force method.
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++)
            {
                Cp[i][j] = 0;
                for (int k = 0; k < p; k++)
                {
                    Cp[i][j] += Ap[i][k] * Bp[k][j];
                    naive_mul += 1;
                }
                naive_add += (p - 1);
            }
        }
    } else
    {
        // Use Strassen's algorithm.
        strassen(Ap, Bp, Cp, p);
    }
    auto end = chrono::steady_clock::now();
    cout << "done." << endl;

    cout << "Elapsed time (matrix multiplication): "
        << chrono::duration_cast<chrono::nanoseconds>(end - start).count()
        / 1000000000.0 << "s" << endl;

    // Prepare and print out result matrix.
    ofstream ofs;
    vector<int> temp(n);
    vector<vector<int>> C(n, temp);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = Cp[i][j];
    if (n < LARGE_N)
    {
        cout << "\nResult:\n\nMatrix A:\n" << endl;
        display_matrix(A, n);
        cout << "\nMatrix B:\n" << endl;
        display_matrix(B, n);
        cout << "\nMatrix C = A x B:\n" << endl;
    } else
    {
        string out_f = OUT_FILE;
        ofs.open(out_f);
        if (!ofs.is_open())
        {
            throw runtime_error("Unable to open file: " + out_f);
            exit(1);
        }
        cout << "Writing output to file: '" << out_f << "'" << endl;
        ofs << "Matrix C = A x B:\n" << endl;
    }
    if (n < LARGE_N)
        display_matrix(C, n);
    else
        display_matrix(C, n, ofs);

    return;
}

static void help_exit(string s)
{
    cerr << "Usage: " << s << " -<n|s> [file]" << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    int n;
    ifstream ifs;
    bool is_naive, from_file;

    if (argc < 2)
        help_exit(argv[0]);
    string user_opt = argv[1];
    if (user_opt.compare("-n") == 0)
    {
        cout << "Method: Naive" << endl;
        is_naive = true;
    } else if (user_opt.compare("-s") == 0)
    {
        cout << "Method: Strassen algorithm" << endl;
        is_naive = false;
    } else
        help_exit(argv[0]);

    if (argc != 3)
    {
        from_file = false;
        cout << "Enter the size of square matrices (n):\n";
        cin >> n;
    } else
    {
        string filename = argv[2];
        cout << "Reading input from " << filename << " file...";
        ifs.open(filename);
        if (!ifs.is_open())
        {
            throw runtime_error("Unable to open file: " + filename);
            exit(1);
        }
        from_file = true;
        ifs >> n;
    }
    if (from_file)
        cout << "done.\n";
    cout << "Matrix size: " << n << " x " << n << endl;

    vector<vector<int>> A, B;

    // Input for matrix A.
    if (!from_file)
        cout << "Enter the first matrix: " << endl;
    for (int i = 0; i < n; i++)
    {
        vector<int> temp;
        for (int j = 0; j < n; j++)
        {
            int i;
            if (!from_file)
                cin >> i;
            else
                ifs >> i;
            temp.push_back(i);
        }
        A.push_back(temp);
    }

    // Input for matrix B.
    if (!from_file)
        cout << "\nEnter the second matrix: " << endl;
    for (int i = 0; i < n; i++)
    {
        vector<int> temp;
        for (int j = 0; j < n; j++)
        {
            int i;
            if (!from_file)
                cin >> i;
            else
                ifs >> i;
            temp.push_back(i);
        }
        B.push_back(temp);
    }

    // Compute the product.
    compute_mul(A, B, n, is_naive);

    // Print operation metrics.
    if (is_naive)
        display_metrics(naive_mul, naive_add);
    else
        display_metrics(strassen_mul, strassen_add);

    return 0;
}

// EOF.
