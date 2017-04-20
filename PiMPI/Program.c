#include <stdio.h>
#include <mpi.h>
#define stag 2001
#define etag 2002
#define nstag 2003
#define rtag 2004

double thread(int, int, double);
double nthDigit(int, double);

int main() {
    MPI_Status status;
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // Print off a hello world message
    if(world_rank == 0) {
        double n = 20000000000;
        double ns = 1.0 / n;
        int load = (n / world_size);
        double pi = 0;
        for(int i = 1; i < world_size; i++) {
            int start = i * load;
            int end = i + load;
            if(end > n) {
                end = n;
            }
            // Send the Parameters to the cluster
            MPI_Send(&start, 1, MPI_INT, i, stag, MPI_COMM_WORLD);
            MPI_Send(&end, 1, MPI_INT, i, etag, MPI_COMM_WORLD);
            MPI_Send(&ns, 1, MPI_INT, i, nstag, MPI_COMM_WORLD);
        }
        for(i = 1; i < world_size; i++) {
            double part;
            // Get the result from the cluster
            MPI_Recv(&part, 1, MPI_DOUBLE, MPI_ANY_SOURCE, rtag, MPI_COMM_WORLD, &status);
            pi += part;
        }
        pi *= ns * pi;
        printf("%lf\n", ns * pi);
    } else {
        int n;
        int m;
        // Compute the value
        MPI_Recv(&n, 1, MPI_DOUBLE, MPI_ANY_SOURCE, stag, MPI_COMM_WORLD, &status);
        MPI_Recv(&m, 1, MPI_DOUBLE, MPI_ANY_SOURCE, etag, MPI_COMM_WORLD, &status);
        MPI_Recv(&ns, 1, MPI_DOUBLE, MPI_ANY_SOURCE, nstag, MPI_COMM_WORLD, &status);
        double res = thread(n, m, ns);
        MPI_Send(&res, 1, MPI_INT, i, rtag, MPI_COMM_WORLD);
    }
    // Finalize the MPI environment.
    MPI_Finalize();
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