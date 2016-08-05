#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <array>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;







/// 1. Specify the constant parameters
///
///
///
///
///
/// todo: Figure out how to use Jocobi module in Eigen package.
/// todo: Try the eigenvalues on the actual Hamiltonian.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///   Section 1 Setup
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///         Section 1.1 Constants
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const double pi = M_PI;                                         //3.1415926535897932384626433832795;
const double theta_12 = 0.579640;                               //rad, sin(theta_12)^2 = 0.3, theta_12 = 0.579640 rad = 33.2109 degrees
const double theta_13 = 0.152245;                               //rad, sin(theta_13)^2 = 0.023, theta_13 = 0.152245
const double delta_m_sq_21 = 7.5 * pow( 10.0 , -5 );            // our value 7.5 * pow(10.0, -5)eV^2
const double delta_m_sq_31 = 2.47 * pow( 10.0 , -3 );           // our value 2.47 * pow(10.0, -3) eV^2
const double epsilon = sqrt(delta_m_sq_21/delta_m_sq_31);
const double s12 = sin(theta_12);
const double c12 = cos(theta_12);
const double s13 = sin(theta_13);
const double c13 = cos(theta_13);
const double default_theta_23=40;         //degrees
const double default_s23 = sin(default_theta_23);
const double default_c23 = cos(default_theta_23);
const double default_eta = 0.0; // parameter for diagonal NSI
const double default_varepsilon = 0.0;  /// This is the modulous of the off-diaonal NSI.
const double default_omega = 0;  /// This is the phase of the off-diagonal NSI.
const double default_delta = 0;  /// This is the CP-violation phase.
const int default_Lkm = 735; //Baseline distance in km.
const double default_alpha_LL = -4; //Lower Limit of alpha (-4)
const double default_alpha_UL = 4.2; //Upper Limit of alpha. 0 is about 7GeV. (1)
const int default_hierarchy = 1; //1 for normal hierarchy, -1 for inverted hierarchy
const int default_neutrino = 1;  //1 for (neutrinos), -1 for (antineutrino)
double density;
const int default_steps = 1000; //how many steps do I want to divide this into.










////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///         Section 1.2 Functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double readDensity(int Lkm);

/// todo: check these function setup
double energy_to_alpha(double energy , double density );

double energy_to_a(double energy, double density);

double alpha_to_energy(double alpha , double density );

double alpha_to_a(double alpha);

Matrix3cd matrixR12(double theta, double delta);
Matrix3cd matrixR13(double theta, double delta);
Matrix3cd matrixR23(double theta, double delta);
Matrix3cd matrixQ1(double delta);
Matrix3cd matrixQ2(double delta);
Matrix3cd matrixQ3(double delta);
Matrix3cd matrixMult( Matrix3cd A, Matrix3cd B );
Matrix3cd matrixMult( Matrix3cd A, Matrix3cd B, Matrix3cd C );
Matrix3cd matrixMult( Matrix3cd A, Matrix3cd B, Matrix3cd C, Matrix3cd D );
Matrix3cd matrixMult( Matrix3cd A, Matrix3cd B, Matrix3cd C, Matrix3cd D, Matrix3cd E );
Matrix3cd matrixMult( Matrix3cd A, Matrix3cd B, Matrix3cd C, Matrix3cd D, Matrix3cd E, Matrix3cd F);
Matrix3cd matrixMultByReal(double a, Matrix3cd B);

Matrix3cd
initializeH(double a, double theta_23, double delta, double eta, double varepsilon, double omega, double hierarchy,
            double neutrino);

Matrix3cd initializeU(double theta_23, double delta);

Matrix3cd diagonalizeH(Matrix3cd matrix, Matrix3cd u);

Matrix3cd diagonalizeU(Matrix3cd matrix, Matrix3cd u);

Matrix3cd prob_exact(double energy, int lkm, Matrix3cd matrix, Matrix3cd u);



int main()
{

    /// Set up toggles for each graph
    bool PlotFigure1a = false;
    bool PlotNumerical = true;







    if (PlotNumerical){
        double eta = default_eta;
        double theta_23 = 40*pi/180;
        double delta = default_delta;
        double varepsilon = 0.1;
        double omega = pi/8;
        double hierarchy = default_hierarchy;
        double neutrino = default_neutrino;
        Matrix3cd H, U, probability;

        int steps = default_steps;
        int Lkm = 1300;
        double density = readDensity(Lkm);
        double energy = 10;  /// Energy = 2 GeV
        double a = energy_to_a(energy, density);
        double min_energy = 0.1; /// Gev
        double max_energy = 50; /// Gev



        H = initializeH(a, theta_23, delta, eta, varepsilon, omega, hierarchy, neutrino);

        cout << H << "\n\n";

        U = initializeU(theta_23, delta);

        H = diagonalizeH(H, U);
        U = diagonalizeU(H, U);

        probability = prob_exact(energy, Lkm, H, U);



    }




    return 0;
}







////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///     Basic matrix setup and manipulations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix3cd matrixR12(double theta, double delta) {
    Matrix3cd result;
    result(0,0) = {cos(theta),0};
    result(0,1) = {sin(theta)*cos(-delta),sin(theta)*sin(-delta)};
    result(0,2) = {0,0};
    result(1,0) = {-sin(theta)*cos(delta),-sin(theta)*sin(delta)};
    result(1,1) = {cos(theta),0};
    result(1,2) = {0,0};
    result(2,0) = {0,0};
    result(2,1) = {0,0};
    result(2,2) = {1,0};
    return result;
}
Matrix3cd matrixR13(double theta, double delta) {
    Matrix3cd result;
    result(0,0) = {cos(theta),0};
    result(0,1) = {0,0};
    result(0,2) = {sin(theta)*cos(-delta),sin(theta)*sin(-delta)};
    result(1,0) = {0,0};
    result(1,1) = {1,0};
    result(1,2) = {0,0};
    result(2,0) = {-sin(theta)*cos(delta),-sin(theta)*sin(delta)};
    result(2,1) = {0,0};
    result(2,2) = {cos(theta),0};
    return result;
}
Matrix3cd matrixR23(double theta, double delta) {
    Matrix3cd result;
    result(0,0) = {1,0};
    result(0,1) = {0,0};
    result(0,2) = {0,0};
    result(1,0) = {0,0};
    result(1,1) = {cos(theta),0};
    result(1,2) = {sin(theta)*cos(-delta),sin(theta)*sin(-delta)};
    result(2,0) = {0,0};
    result(2,1) = {-sin(theta)*cos(delta),-sin(theta)*sin(delta)};
    result(2,2) = {cos(theta),0};
    return result;
}
Matrix3cd matrixQ1(double delta)
{
    Matrix3cd result;
    result = Vector3cd({cos(delta),sin(delta)},1,1).asDiagonal();
    return result;
}
Matrix3cd matrixQ2(double delta)
{
    Matrix3cd result;
    result = Vector3cd(1,{cos(delta),sin(delta)},1).asDiagonal();
    return result;
}
Matrix3cd matrixQ3(double delta)
{
    Matrix3cd result;
    result = Vector3cd(1, 1, {cos(delta),sin(delta)}).asDiagonal();
    return result;
}
Matrix3cd matrixMult( Matrix3cd A, Matrix3cd B )
{
    Matrix3cd C;
    C = A * B;
    return C;
}
Matrix3cd matrixMult(Matrix3cd A, Matrix3cd B, Matrix3cd C)
{
    Matrix3cd D;
    D = matrixMult( A, matrixMult( B, C) );
    return D;
}
Matrix3cd matrixMult(Matrix3cd A, Matrix3cd B, Matrix3cd C, Matrix3cd D )
{
    Matrix3cd E;
    E = matrixMult( A, matrixMult( B, matrixMult( C, D) ) );
    return E;
}
Matrix3cd matrixMult(Matrix3cd A, Matrix3cd B, Matrix3cd C, Matrix3cd D, Matrix3cd E )
{
    Matrix3cd F;
    F = matrixMult( A, matrixMult( B, matrixMult( C, matrixMult( D, E) ) ) );
    return F;
}
Matrix3cd matrixMult( Matrix3cd A, Matrix3cd B, Matrix3cd C, Matrix3cd D, Matrix3cd E, Matrix3cd F)
{
    Matrix3cd G;
    G = matrixMult( A, matrixMult( B, matrixMult( C, matrixMult( D, matrixMult( E, F)))));
    return G;
}
Matrix3cd matrixMultByReal (double x, Matrix3cd A )
{
    return x * A;
}







////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double alpha_to_energy ( double alpha , double density )
{
    double a = delta_m_sq_31 * pow(epsilon, -alpha);
    double result = a/(7.63e-05 * density);
    return result;
}
///
/// \param energy Energy of the neutrino beam
/// \param density Average density along this particular baseline
/// \return

double energy_to_alpha ( double energy , double density )
{
    double result = 7.63e-05 * density * energy / delta_m_sq_31;
    result = - log(result) / log(epsilon);
    return result;
}

double alpha_to_a(double alpha)
{
    double a = delta_m_sq_31 * pow(epsilon, -alpha);
    return a;
}

double energy_to_a(double energy, double density)
{
    double a = energy * density * 7.63e-05;
    return a;
}

Matrix3cd
initializeH(double a, double theta_23, double delta, double eta, double varepsilon, double omega, double hierarchy,
            double neutrino) {
    Matrix3cd H_mass, H_NSI, U, Q, H, temp;

    H_mass = Vector3cd( 0, delta_m_sq_21, (hierarchy*delta_m_sq_31) ).asDiagonal();
    H_NSI(0,0) = {1,0};
    H_NSI(0,1) = {0,0};
    H_NSI(0,2) = {0,0};
    H_NSI(1,0) = {0,0};
    H_NSI(1,1) = {eta,0};
    H_NSI(1,2) = {varepsilon*cos(omega),varepsilon*sin(omega)};
    H_NSI(2,0) = {0,0};
    H_NSI(2,1) = {varepsilon*cos(-omega),varepsilon*sin(-omega)};
    H_NSI(2,2) = {eta,0};


    U = matrixMult( matrixR23( theta_23, 0) , matrixR13( theta_13, delta), matrixR12( theta_12, 0) );
    Q = matrixQ3(delta);
    H = matrixMult( Q.adjoint(),  U.adjoint(), H_NSI, U, Q);
    H = matrixMultByReal(neutrino * a, H);
    temp = H_mass + H;
    H = temp;
    return H;
}













Matrix3cd prob_exact(double energy, int lkm, Matrix3cd matrix, Matrix3cd u) {
    return Eigen::Matrix3cd();
}

Matrix3cd diagonalizeU(Matrix3cd matrix, Matrix3cd u) {
    return Eigen::Matrix3cd();
}

Matrix3cd diagonalizeH(Matrix3cd matrix, Matrix3cd u) {
    return Eigen::Matrix3cd();
}

Matrix3cd initializeU(double theta_23, double delta) {
    return Eigen::Matrix3cd();
}




double readDensity(int Lkm) {
    ifstream densityFile("/Users/ykao/ClionProjects/MatrixDiagonalization/average.dat");
    int readLkm;
    double readDensity;
    while (densityFile >> readLkm >> readDensity) {
        if (readLkm == Lkm) {
            density = readDensity; //  Average density in grams per cm^3.
        }
    }
    return density;
}





