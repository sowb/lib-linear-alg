#include <iostream>
#include "../include/matrix.h"
#include <string>
using namespace std;

int main()
{
    //Load the iris data
    Matrix Iris("iris.data", ",", false);

    // Check dim
    print("\nThe dataset dimension");
    Iris.shape();

    // Print the data
    Iris.head(5);

    // Scale the data, perform PCA
    Iris = scale(Iris, "col");

    Matrix Comp, Z, ExpVar;
    PCA(Iris, &Comp, &Z, &ExpVar, 2);

    //Print the results
    print("\nPrincipal components");
    Comp.head(5);
    Z.print("Singular Values");
    ExpVar.print("Explained variance");

    // Save the results
    Comp.to_csv("iris-pc.csv");
    Z.to_csv("iris-sing-vals.csv");
    ExpVar.to_csv("iris-exp-var.csv");

    return 0;
}
