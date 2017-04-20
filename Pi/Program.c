#include <stdio.h>

double thread(int, int, double);
double nthDigit(int, double);

int main() {
    double n = 2000000000;
    double ns = 1.0 / n;
    printf("%lf\n", ns * thread(0, n, ns));
}

// a thread to consolidate pi n digits from n to m
double thread(int n, int m, double ns) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for(int i = n; i < m; i++) {
        sum += nthDigit(i, ns);
    }
    return sum;
}

// pi function
double nthDigit(int n, double ns) {
    double x = (n + 0.5) * ns;
    return 4.0 / (1.0 + x * x);
}