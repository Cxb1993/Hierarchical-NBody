#include <cstdio>
#include "MLFMM.h"
#include "BHNode.h"

clock_t timer;
void tic() { timer = clock(); }
double toc() { return (double)(clock() - timer) / CLOCKS_PER_SEC + 0.001; }

void PrintFMMHeader() {
    printf("%3s %4s %8s %10s %10s %10s %10s %10s %10s\n", 
        "l", "p", "N", "FLOP", "t_direct", "t_FMM", "FLOPS", "Abs. Err", "Rel. Err");
}

void RunFMM(int levels, int degree, int N) 
{
    Potential coulomb(degree);
    MLFMM tree(levels, coulomb);

    std::vector<Point*> sources(N);
    std::vector<Point*> targets(N);
    std::vector<double> exact(N, 0.0);
    std::vector<double> approx(N, 0.0);

    double directTime = 0;
    double approxTime = 0;

    for (int index = 0; index < N; index++) {
        double x = randf();
        double y = randf();
        sources[index] = new Point(Complex(x, y), index);
        tree.AddSource(sources[index]);
        targets[index] = new Point(Complex(x, y), index);
        tree.AddTarget(targets[index]);
    }
    
    tic();
    tree.DirectSolve();
    directTime = toc();
    
    for (int index = 0; index < N; index++)
        exact[index] = targets[index]->potential;

#ifdef PIECHART

    tree.flops = 0;
    long s1,s2,s3,s4,s5; 
    tree.MultipoleExpansion(); s1 = tree.flops; tree.flops = 0;
    tree.MultipoleToMultipoleTranslation(); s2 = tree.flops; tree.flops = 0;
    tree.MultipoleToLocalTranslation(); s3 = tree.flops; tree.flops = 0; 
    tree.LocalToLocalTranslation(); s4 = tree.flops; tree.flops = 0;
    tree.LocalExpansion(); s5 = tree.flops; tree.flops = 0;
    long sum = s1+s2+s3+s4+s5;
    printf("%ld %ld %ld %ld %ld %ld\n", s1, s2, s3, s4, s5, sum);

#endif

    tic();
    tree.Solve();
    approxTime = toc();

    for (int index = 0; index < N; index++)
        approx[index] = targets[index]->potential;

    double absError = AvgAbsError(exact, approx);
    double relError = AvgRelError(exact, approx);
    double fps = (double)tree.flops / approxTime;

    printf("%3d %4d %8d %10ld %10.3f %10.3f %10.2e %10.2e %10.2e\n", levels, degree, N, tree.flops, directTime, approxTime, fps, absError, relError);
    fflush(stdout);

    // clean up
    for (int index = 0; index < N; index++){
        delete sources[index];
        delete targets[index];
    }
}

void TestFMMPerformance() {
    PrintFMMHeader();
    double tolerance = 1.0e-2;
    for (int i = 0; i < 51; i++) {
        double exponent = (0.06 * i) + 2;
        int N = round(pow(10, exponent));
        int levels = round(log((double)N) / log(4.0)); 
        if (levels <= 3) levels = 3;
        int p = round(log(1.0/tolerance) / log(2));
        RunFMM(levels, p, N);
    }
}

void TestFMMScaling() {
    PrintFMMHeader();
    for (int i = 0; i < 31; i++) {
        double exponent = (0.07 * i) + 2;
        int N = round(pow(10, exponent));
        int levels = 7;
        int p = 10;
        RunFMM(levels, p, N);
    }
}

void TestFMMLevels() {
    PrintFMMHeader();
    for (int level = 4; level < 13; level++) {
        int N = pow(4, level);
        int p = 2;
        RunFMM(level, p, N);
    }
}

void TestDelicious() {
    RunFMM(5, 6, 1024);
    RunFMM(6, 6, 4096);
    RunFMM(5, 12, 1024);
    RunFMM(6, 12, 4096);
    RunFMM(5, 6, 1024);
    RunFMM(5, 6, 4096);
}

long BHNode::flops = 0; 
void RunBH(int maxDepth, int N, double theta) {

    BHNode tree(Complex(0.5, 0.5), Complex(0.5, 0.5), 0, maxDepth);

    std::vector<Point*> sources(N);
    std::vector<Point*> targets(N);

    std::vector<double> exact(N, 0.0);
    std::vector<double> approx(N, 0.0);

    double directTime = 0;
    double approxTime = 0;

    for (int i = 0; i < N; i++) {
        double x = randf();
        double y = randf();
        sources[i] = new Point(Complex(x, y), i);
        tree.AddSource(sources[i]);
        targets[i] = new Point(Complex(x, y), i);
    }

    tic();
    for (int i = 0; i < N; i++){
        exact[i] = tree.ComputePotentialDirect(sources, targets[i]);
    }
    directTime = toc();

    BHNode::flops = 0;

    tic();
    tree.ComputeChargeDistribution();
    for (int i = 0; i < N; i++){
        approx[i] = tree.ComputePotential(targets[i], theta);
    }
    approxTime = toc();

    double error = AvgAbsError(exact, approx);
    double FLOPS = (double)BHNode::flops / (approxTime + 0.001);
    double speedup = (double)N * N / BHNode::flops;
    printf("%8d %5d %5.3f %10ld %10.3f %10.3f %10.3f %10.2e %10.2e\n",
        N, maxDepth, theta, BHNode::flops, speedup, directTime, approxTime, FLOPS, error);
    fflush(stdout);

    // clean up
    for (int index = 0; index < N; index++){
        delete sources[index];
        delete targets[index];
    }
}

void PrintBHHeader() {
    printf("%8s %5s %5s %10s %10s %10s %10s %10s %10s\n",
        "N", "depth", "theta", "FLOP", "speedup", "t_direct", "t_approx", "FLOPS", "error");
    fflush(stdout);
}

void TestBHN() {
    PrintBHHeader();
    for (int i = 0; i < 51; i++) {
        double exponent = (0.06 * i) + 2;
        int N = round(pow(10, exponent));
        int levels = round(log((double)N) / log(4.0) + 0.5); 
        if (levels <= 3) levels = 3;
        RunBH(levels, N, 4.0);
    }
}

void TestBHTheta() {
    PrintBHHeader();
    for (int i = 0; i < 51; i++) {
        double exponent = (0.06 * i) + 2;
        int N = round(pow(10, exponent));
        int levels = round(log((double)N) / log(4.0) + 0.5); 
        if (levels <= 3) levels = 3;
        RunBH(levels, N, 4.0);
    }
}

int main(int argc, char** argv) 
{
    return 0;
}
