/* Generates a grid of random numbers. */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

void print_rand_grid(int k, int m, int n, int r)
{
    int i;
    for (i = 0; i < k; i++)
    {
        printf("%d ", (rand() % (n - m + 1)) + m);
        if ((i + 1) % r == 0) printf("\n");
    }
    if ((k % r) != 0) printf("\n");
    return;
}

void msg_exit(char *s)
{
    fprintf(stderr, "%s", s);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    int k, m, n, r;

    if (argc != 5)
    {
        fprintf(stderr, "Usage: %s <size> <l_bound> <u_bound> <row_len>\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    k = atoi(argv[1]);
    if (k <= 0)
        msg_exit("Size must be > 0\n");

    m = atoi(argv[2]);
    if (m < 0)
        msg_exit("Lower bound must be >= 0\n");

    n = atoi(argv[3]);
    if (n < m)
        msg_exit("Upper bound must be >= lower bound\n");

    r = atoi(argv[4]);
    if (r <= 0 || r > k)
        msg_exit("Row length must be > 0 and <= size\n");

    srand(time(NULL));
    print_rand_grid(k, m, n, r);

    return EXIT_SUCCESS;
}

/* EOF. */
