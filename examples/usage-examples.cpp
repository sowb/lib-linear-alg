#include <iostream>
#include "../include/matrix.h"
using namespace std;

int main()
{
    // Generate 5x4 random matrix
    Matrix A(5, 4, unif(-1, 1)); //rnorm - for the normal distribution
    A.print("matrix A:");
    transpose(A).print("Transpose A:");
    // Create matrix 4x4 matrix ;
    Matrix B = {{2, 0, 8, 10, 7},
                {3, 41, 9, 4, 1},
                {1, 5, 5, 86, 4}};

    B.print("matrix B:");

    // Transpose matrix
    Matrix C = B.T() * B;
    C.print("B.T * B:");

    //Get the Diagonal of the matrix
    Matrix D = diag(B);
    D.print("Diagonal of D:");

    // Create a diagonal matrix
    Matrix Dmat = diag_matrix(D);
    Dmat.print("\nDiagonal Matrix:");

    //QR factorization
    Matrix E = {{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}};
    E.print("E:");
    print("\nQR factorization of E:\n");
    Matrix Q, R;
    QR_factorization(E, &Q, &R);
    Q.print("Q:");
    R.print("R:");
    (Q * R).print("Check QR = B");

    // eigen decomposition
    print("\nEigen decomposition of E");
    Matrix EigVal, EigVec;
    eigen_decomposition(E, &EigVal, &EigVec); // use the basic QR
    EigVal.print("EigVal:");
    EigVec.print("EigVec");

    // Eigen decomposition using QR with shift
    QR_algorithm_w_shift(E, &EigVal, &EigVec);
    EigVal.print("EigVal-QR shift:");
    EigVec.print("EigVec-QR shift");

    //Scaling a matrix : (x-x_bar)/std(x)
    Matrix F = scale(B, "col");
    F.print("B scaled:");

    // Singular Value Decomposition
    print("\nSingular Value Decomposition");
    Matrix U, S, V;
    SVD(F, &U, &S, &V);
    U.print("U:");
    S.print("S:");
    V.print("V:");

    //Principal Component Analysis
    print("\nPrincipal Components Analysis");
    Matrix Comp, Z, Var;

    //scale B
    PCA(F, &Comp, &Z, &Var, 2);
    Z.print("Singular Values :");
    Var.print("Explained Variance :");
    Comp.print("Components");

    // Save the components in a file
    //Comp.to_csv("components");

    return 0;
}
