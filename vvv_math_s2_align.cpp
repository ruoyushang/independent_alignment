
#include <iostream>
#include <math.h>
#include <complex>
#include "string"
#include "fstream"
//#include "/Users/rshang/pSCT/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
//#include "/Users/rshang/pSCT/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;
//using Eigen::MatrixXd;

const int N_panels = 16;

MatrixXcd fillMatrix(std::string name)
{
    MatrixXcd matrix(6,6);
    std::ifstream file;
    file.open ("S2_response_matrix.txt");
    if (!file.is_open()) {
        std::cout << "file not found!!" << std::endl;
        return matrix;
    }

    std::string word;
    while (file >> word)
    {
        if (word.find(name) != std::string::npos) {
            std::cout<< word << '\n';
            for (int i=0;i<6;i++) {
                for (int j=0;j<6;j++) {
                    file >> word;
                    matrix(i,j) = atof(word.c_str());
                }
            }
            break;
        }
    }
    return matrix;
}
MatrixXcd Hermitian(double phase, MatrixXcd E)
{
    const std::complex<double> If(0.0, 1.0);
    MatrixXcd H(6,6);
    MatrixXcd I(6,6);
    I << 1.0, 0., 0., 0., 0., 0.,
         0., 1.0, 0., 0., 0., 0.,
         0., 0., 1.0, 0., 0., 0.,
         0., 0., 0., 1.0, 0., 0.,
         0., 0., 0., 0., 1.0, 0.,
         0., 0., 0., 0., 0., 1.0;
    H = E.transpose()*E - cos(phase)*(E.transpose()+E) - If*sin(phase)*(E.transpose()-E) + I*(cos(phase)+If*sin(phase))*(cos(phase)-If*sin(phase));
    return H;
}
MatrixXcd AvgHermitian(double phase, std::vector<MatrixXcd> E)
{
    std::vector<MatrixXcd> H;
    for (int i=0;i<N_panels;i++) {
        H.push_back(MatrixXcd(6,6));
        H.at(i) = Hermitian(phase, E.at(i));
    }
    MatrixXcd Havg(6,6);
    Havg << 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0.;
    for (int i=0;i<N_panels;i++) {
        Havg += 1./double(N_panels)*H.at(i);
    }
    return Havg;
}
MatrixXcd BuildQmatrix(std::vector<MatrixXcd> E, int panel)
{
    MatrixXcd Q(N_panels*6,N_panels*6);
    MatrixXcd I(6,6);
    I << 1.0, 0., 0., 0., 0., 0.,
         0., 1.0, 0., 0., 0., 0.,
         0., 0., 1.0, 0., 0., 0.,
         0., 0., 0., 1.0, 0., 0.,
         0., 0., 0., 0., 1.0, 0.,
         0., 0., 0., 0., 0., 1.0;
    std::cout << "Begin calculating Q matrix." << std::endl;
    for (int i=0;i<N_panels;i++) {
        int j = i;
        if (panel>=0) j = panel;
        if (i!=N_panels-1) {
            Q.block(i*6,i*6,6,6) = I+E.at(j).transpose()*E.at(j);
            Q.block(i*6,(i+1)*6,6,6) = E.at(j).transpose();
            Q.block((i+1)*6,i*6,6,6) = E.at(j);
        }
        else {
            Q.block((N_panels-1)*6,(N_panels-1)*6,6,6) = I+E.at(j).transpose()*E.at(j);
            Q.block((N_panels-1)*6,0,6,6) = E.at(j).transpose();
            Q.block(0,(N_panels-1)*6,6,6) = E.at(j);
        }
    }
    std::cout << "Done calculating Q matrix." << std::endl;
    return Q;
}
VectorXd BigDeltaVector() 
{
    VectorXd Delta(6);
    Delta(0) = -121.988;
    Delta(1) = 175.078;
    Delta(2) = 145.356;
    Delta(3) = 349.083;
    Delta(4) = -38.090;
    Delta(5) = -442.441;
    //Delta(0) = 121.988;
    //Delta(1) = -175.078;
    //Delta(2) = -145.356;
    //Delta(3) = -349.083;
    //Delta(4) = 38.090;
    //Delta(5) = 442.441;
    return Delta;
}
VectorXd ActuatorVector(int panel) 
{
    VectorXd act(6);
    if (panel==0) {
    act(0)  = 437.788;
    act(1)  = 433.22;
    act(2)  = 434.316;
    act(3)  = 434.731;
    act(4)  = 438.934;
    act(5)  = 437.814;
    }
    if (panel==1) {
    act(0)  = 440.363;
    act(1)  = 439.662;
    act(2)  = 438.245;
    act(3)  = 436.43;
    act(4) = 439.481;
    act(5) = 437.166;
    }
    if (panel==2) {
    act(0) = 446.181;
    act(1) = 445.907;
    act(2) = 440.75;
    act(3) = 442.111;
    act(4) = 443.495;
    act(5) = 441.512;
    }
    if (panel==3) {
    act(0) = 450.694;
    act(1) = 435.025;
    act(2) = 441.132;
    act(3) = 439.058;
    act(4) = 438.54;
    act(5) = 445.473;
    }
    if (panel==4) {
    act(0) = 435.939;
    act(1) = 429.865;
    act(2) = 435.411;
    act(3) = 437.265;
    act(4) = 438.337;
    act(5) = 436.669;
    }
    if (panel==5) {
    act(0) = 437.431;
    act(1) = 437.416;
    act(2) = 433.252;
    act(3) = 434.347;
    act(4) = 431.887;
    act(5) = 432.346;
    }
    if (panel==6) {
    act(0) = 440.923;
    act(1) = 445.227;
    act(2) = 443.197;
    act(3) = 436.226;
    act(4) = 435.398;
    act(5) = 444.572;
    }
    if (panel==7) {
    act(0) = 446.292;
    act(1) = 444.795;
    act(2) = 441.109;
    act(3) = 440.7;
    act(4) = 440.289;
    act(5) = 450.505;
    }
    return act;
}
//VectorXd DeltaSigmaVector() 
//{
//    // Jan08
//    VectorXd DeltaSigma(48);
//    DeltaSigma(0)  =   2.22308;
//    DeltaSigma(1)  =  0.383084;
//    DeltaSigma(2)  = -0.501165;
//    DeltaSigma(3)  =   1.02764;
//    DeltaSigma(4)  =  0.390244;
//    DeltaSigma(5)  =  -4.64975;
//    DeltaSigma(6)  =  0.275801;
//    DeltaSigma(7)  =   10.3989;
//    DeltaSigma(8)  =  -3.80159;
//    DeltaSigma(9)  =   8.27924;
//    DeltaSigma(10) =   1.46952;
//    DeltaSigma(11) =  -7.62038;
//    DeltaSigma(12) =  -2.84122;
//    DeltaSigma(13) =   9.35395;
//    DeltaSigma(14) = -0.359173;
//    DeltaSigma(15) =   3.86903;
//    DeltaSigma(16) =  -3.32244;
//    DeltaSigma(17) =  -5.44716;
//    DeltaSigma(18) =   3.08892;
//    DeltaSigma(19) =  -4.73686;
//    DeltaSigma(20) = -0.457365;
//    DeltaSigma(21) =   3.54856;
//    DeltaSigma(22) =  -3.61201;
//    DeltaSigma(23) =   5.31459;
//    DeltaSigma(24) =   4.63734;
//    DeltaSigma(25) =  -12.2238;
//    DeltaSigma(26) =   5.22802;
//    DeltaSigma(27) =  -6.40906;
//    DeltaSigma(28) =     7.141;
//    DeltaSigma(29) =  -3.99606;
//    DeltaSigma(30) =   7.16758;
//    DeltaSigma(31) =   3.56684;
//    DeltaSigma(32) =   10.9275;
//    DeltaSigma(33) =    8.7883;
//    DeltaSigma(34) =   26.1361;
//    DeltaSigma(35) =  -33.4946;
//    DeltaSigma(36) =  -10.6336;
//    DeltaSigma(37) =   46.4568;
//    DeltaSigma(38) =   41.8835;
//    DeltaSigma(39) =   75.3998;
//    DeltaSigma(40) =   17.6161;
//    DeltaSigma(41) =  -104.142;
//    DeltaSigma(42) =  -2.56872;
//    DeltaSigma(43) =   20.0496;
//    DeltaSigma(44) =  -47.6552;
//    DeltaSigma(45) =  -1.24372;
//    DeltaSigma(46) =  -39.7564;
//    DeltaSigma(47) =   43.2116;
//           
//    return -1.*DeltaSigma;
//}          
//VectorXd DeltaSigmaVector() 
//{
//    // Jan09
//    VectorXd DeltaSigma(48);
//    DeltaSigma(0)  =      0.;       
//    DeltaSigma(1)  =      0.;
//    DeltaSigma(2)  =      0.;       
//    DeltaSigma(3)  =      0.;
//    DeltaSigma(4)  =      0.;       
//    DeltaSigma(5)  =      0.;
//    DeltaSigma(6)  =   0.439497;
//    DeltaSigma(7)  =    4.05985;
//    DeltaSigma(8)  =  -0.104039;
//    DeltaSigma(9)  =   -2.80682;
//    DeltaSigma(10) =   -2.36359;
//    DeltaSigma(11) =   0.587471;
//    DeltaSigma(12) =   1.63181;
//    DeltaSigma(13) =   10.2016;
//    DeltaSigma(14) =   3.31363;
//    DeltaSigma(15) =   19.9294;
//    DeltaSigma(16) =   6.84334;
//    DeltaSigma(17) =  -6.44434;
//    DeltaSigma(18) =  -3.32298;
//    DeltaSigma(19) =   -4.2736;
//    DeltaSigma(20) =  -6.09279;
//    DeltaSigma(21) =  -34.3185;
//    DeltaSigma(22) =  -3.80748;
//    DeltaSigma(23) =   13.4792;
//    DeltaSigma(24) =  -6.89052;
//    DeltaSigma(25) =   10.5784;
//    DeltaSigma(26) =   4.30211;
//    DeltaSigma(27) =   10.9598;
//    DeltaSigma(28) =   1.34734;
//    DeltaSigma(29) =  -30.4773;
//    DeltaSigma(30) =   1.80516;
//    DeltaSigma(31) =  -8.01165;
//    DeltaSigma(32) =  -8.53547;
//    DeltaSigma(33) =  -16.6432;
//    DeltaSigma(34) =  -6.18572;
//    DeltaSigma(35) =    4.0446;
//    DeltaSigma(36) =  -6.53139;
//    DeltaSigma(37) =  -6.98123;
//    DeltaSigma(38) =  -15.4135;
//    DeltaSigma(39) =  -25.0957;
//    DeltaSigma(40) =  -7.46383;
//    DeltaSigma(41) =   13.6368;
//    DeltaSigma(42) =   10.5425;
//    DeltaSigma(43) =   19.0884;
//    DeltaSigma(44) =   13.2412;
//    DeltaSigma(45) =   27.9662;
//    DeltaSigma(46) =   8.95851;
//    DeltaSigma(47) =  -3.96427;
//
//    return -1.*DeltaSigma;
//}          
VectorXd DeltaSigmaVector() 
{
    // Jan21
    VectorXd DeltaSigma(N_panels*6);
    for (int i=0;i<N_panels*6;i++)
    {
        DeltaSigma(i)  =      0.;       
    }

    DeltaSigma(90)  =      19.0324;       
    DeltaSigma(91)  =      18.2996;
    DeltaSigma(92)  =      23.4634;       
    DeltaSigma(93)  =      31.0742;
    DeltaSigma(94)  =      18.9542;       
    DeltaSigma(95)  =      -6.35156;


    return -1.*DeltaSigma;
}          
VectorXd BuildNoiseVector(VectorXd dS, std::vector<MatrixXcd> Et) 
{
    VectorXd noise(N_panels*6);
    noise = VectorXd::Zero(N_panels*6);
    VectorXd dSigma(N_panels*6);
    dSigma = VectorXd::Zero(N_panels*6);
    //dSigma.segment(42,6) = BigDeltaVector();
    dSigma = DeltaSigmaVector();
    for (int i=0;i<N_panels;i++) {
        if (i<N_panels-1) {
            noise.segment(i*6,6) = Et.at(i).real()*dS.segment(i*6,6)+dS.segment((i+1)*6,6)+dSigma.segment(i*6,6);
        }
        else {
            noise.segment(i*6,6) = Et.at(i).real()*dS.segment(i*6,6)+dS.segment((0)*6,6)+dSigma.segment(i*6,6);
        }
    }
    return noise;
}
VectorXd BuilddSmatrixIterative(MatrixXcd Q, VectorXd dS, int panel, int dim) 
{
    ComplexEigenSolver<MatrixXcd> eigensolver(Q);
    VectorXd dS_new(dim);
    dS_new = VectorXd::Zero(dim);
    VectorXd dS_final(dim);
    dS_final = VectorXd::Zero(dim);
    MatrixXd e_vectors(dim,dim);
    VectorXd lambda(dim);
    e_vectors = eigensolver.eigenvectors().real();
    lambda = eigensolver.eigenvalues().real();
    for (int i=0;i<dim;i++) {
        double projection = dS.transpose().dot(e_vectors.col(i));
        if (abs(lambda(i))<1.0) dS_new += e_vectors.col(i)*projection;
    }
    dS_final = -1.*(Q*dS_new).real();
    return dS_final;
}
VectorXd BuilddSmatrix(MatrixXcd Q, MatrixXcd dV, int panel, int dim) 
{
    ComplexEigenSolver<MatrixXcd> eigensolver(Q);
    VectorXd dS(dim);
    dS = VectorXd::Zero(dim);
    VectorXd dV_real(dim);
    VectorXd lambda(dim);
    MatrixXd e_vectors(dim,dim);
    lambda = eigensolver.eigenvalues().real();
    e_vectors = eigensolver.eigenvectors().real();
    dV_real = dV.real();
    for (int i=0;i<dim;i++) {
        if (i==0) continue;
        if (i==1) continue;
        double projection = -1.*dV_real.transpose().dot(e_vectors.col(i));
        //if (abs(projection/lambda(i))>200.) projection = 0.;
        //if (i<4) projection = 0.;
        std::cout << i << ", lambda " << lambda(i) << ", dV^{T}*e " << -projection << ", dV^{T}*e/lambda " << -projection/lambda(i) << std::endl;
        dS += e_vectors.col(i)*projection/lambda(i);
    }
    return dS;
}
MatrixXcd BuilddVmatrix(std::vector<MatrixXcd> E, int panel)
{
    MatrixXcd W(N_panels*6,N_panels*6);
    MatrixXcd W_reduced((N_panels-1)*6,N_panels*6);
    MatrixXcd dV((N_panels-1)*6,1);
    //MatrixXcd dSigma(48,1);
    VectorXd dSigma(N_panels*6);
    MatrixXcd I(6,6);
    I << 1.0, 0., 0., 0., 0., 0.,
         0., 1.0, 0., 0., 0., 0.,
         0., 0., 1.0, 0., 0., 0.,
         0., 0., 0., 1.0, 0., 0.,
         0., 0., 0., 0., 1.0, 0.,
         0., 0., 0., 0., 0., 1.0;
    std::cout << "Begin calculating dV matrix." << std::endl;
    for (int i=0;i<N_panels;i++) {
        int j = i;
        if (panel>=0) j = panel;
        if (i!=(N_panels-1)) {
            W.block(i*6,i*6,6,6) = E.at(j).transpose();
            W.block((i+1)*6,i*6,6,6) = I;
        }
        else {
            W.block((N_panels-1)*6,(N_panels-1)*6,6,6) = E.at(j).transpose();
            W.block(0,(N_panels-1)*6,6,6) = I;
        }
    }
    W_reduced = W.block(6,0,(N_panels-1)*6,N_panels*6);
    //dSigma = MatrixXcd::Zero(48,1);
    //dSigma.block(42,0,6,1) = BigDeltaVector();
    dSigma = VectorXd::Zero(N_panels*6);
    dSigma = DeltaSigmaVector();
    std::cout << "Done calculating W and dSigma matrix." << std::endl;
    dV = W_reduced*dSigma;
    std::cout << "Done calculating dV matrix." << std::endl;
    return dV;
}
MatrixXcd BuildUmatrix(MatrixXcd H1a, MatrixXcd H1b, MatrixXcd H2a, MatrixXcd H2b, MatrixXcd H3a, MatrixXcd H3b)
{
    MatrixXcd U(6,6);

    ComplexEigenSolver<MatrixXcd> eigensolver1a(H1a);
    ComplexEigenSolver<MatrixXcd> eigensolver2a(H2a);
    ComplexEigenSolver<MatrixXcd> eigensolver3a(H3a);
    ComplexEigenSolver<MatrixXcd> eigensolver1b(H1b);
    ComplexEigenSolver<MatrixXcd> eigensolver2b(H2b);
    ComplexEigenSolver<MatrixXcd> eigensolver3b(H3b);
    
    U.col(0) = eigensolver1a.eigenvectors().col(0);
    U.col(1) = eigensolver1b.eigenvectors().col(1);
    U.col(2) = eigensolver2a.eigenvectors().col(0);
    U.col(3) = eigensolver2b.eigenvectors().col(0);
    U.col(4) = eigensolver3a.eigenvectors().col(1);
    U.col(5) = eigensolver3b.eigenvectors().col(1);

    return U;
}
int main()
{
    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    std::vector<MatrixXcd> MR;
    MR.push_back(MatrixXcd(6,6));
    MR.at(0) = fillMatrix("M0R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(1) = fillMatrix("M1R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(2) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(3) = fillMatrix("M3R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(4) = fillMatrix("M4R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(5) = fillMatrix("M5R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(6) = fillMatrix("M6R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(7) = fillMatrix("M7R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(8) = fillMatrix("M8R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(9) = fillMatrix("M9R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(10) = fillMatrix("M10R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(11) = fillMatrix("M11R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(12) = fillMatrix("M12R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(13) = fillMatrix("M13R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(14) = fillMatrix("M14R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(15) = fillMatrix("M15R");

    std::vector<MatrixXcd> ML;
    ML.push_back(MatrixXcd(6,6));
    ML.at(0) = fillMatrix("M0L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(1) = fillMatrix("M1L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(2) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(3) = fillMatrix("M3L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(4) = fillMatrix("M4L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(5) = fillMatrix("M5L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(6) = fillMatrix("M6L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(7) = fillMatrix("M7L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(8) = fillMatrix("M8L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(9) = fillMatrix("M9L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(10) = fillMatrix("M10L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(11) = fillMatrix("M11L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(12) = fillMatrix("M12L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(13) = fillMatrix("M13L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(14) = fillMatrix("M14L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(15) = fillMatrix("M15L");

    //std::vector<MatrixXcd> MR;
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(0) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(1) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(2) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(3) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(4) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(5) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(6) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(7) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(8) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(9) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(10) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(11) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(12) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(13) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(14) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(15) = fillMatrix("M3R");

    //std::vector<MatrixXcd> ML;
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(0) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(1) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(2) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(3) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(4) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(5) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(6) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(7) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(8) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(9) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(10) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(11) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(12) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(13) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(14) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(15) = fillMatrix("M3L");


    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    std::cout << "Here we compute Em(i) = MR(i)*ML(i)^{-1}" << std::endl;
    std::vector<MatrixXcd> Em;
    VectorXd Em_eigenvalue_real(6);
    VectorXd Em_eigenvalue_imag(6);
    for (int i=0;i<N_panels;i++) {
        Em.push_back(MatrixXcd(6,6));
        Em.at(i) = MR.at(i)*ML.at(i).inverse();
        //std::cout << "Em(" << i << "):" << std::endl;
        //std::cout << Em.at(i) << std::endl;
        ComplexEigenSolver<MatrixXcd> eigensolver_Em(Em.at(i));
        Em_eigenvalue_real = eigensolver_Em.eigenvalues().real();
        Em_eigenvalue_imag = eigensolver_Em.eigenvalues().imag();
        std::cout << "The eigenvalues of Em("<<i<<") are:" << std::endl;
        for (int j=0;j<6;j++)
        {
            std::cout << pow(Em_eigenvalue_real(j)*Em_eigenvalue_real(j)+Em_eigenvalue_imag(j)*Em_eigenvalue_imag(j),0.5) << ", " << atan2(Em_eigenvalue_imag(j), Em_eigenvalue_real(j))/M_PI << "pi" << std::endl;
        }
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    double Phase_Choice_0 = M_PI;
    double Phase_Choice_1 = 6.*M_PI/8.;
    double Phase_Choice_2 = 7.*M_PI/8.;

    //std::cout << "Here we compute H_avg = 1/8 * sum (H(i))" << std::endl;
    //MatrixXcd Havg1(6,6);
    //MatrixXcd Havg2(6,6);
    //MatrixXcd Havg3(6,6);
    //Havg1 = AvgHermitian(M_PI, Em);
    //Havg2 = AvgHermitian(Phase_Choice, Em);
    //Havg3 = AvgHermitian(-Phase_Choice, Em);
    //ComplexEigenSolver<MatrixXcd> eigensolver1(Havg1);
    //std::cout << "The eigenvalues of avg H(pi) are:\n" << eigensolver1.eigenvalues() << std::endl;
    //ComplexEigenSolver<MatrixXcd> eigensolver2(Havg2);
    //std::cout << "The eigenvalues of avg H(3pi/4) are:\n" << eigensolver2.eigenvalues() << std::endl;
    //ComplexEigenSolver<MatrixXcd> eigensolver3(Havg3);
    //std::cout << "The eigenvalues of avg H(-3pi/4) are:\n" << eigensolver3.eigenvalues() << std::endl;

    std::vector<MatrixXcd> H1a;
    std::vector<MatrixXcd> H2a;
    std::vector<MatrixXcd> H3a;
    std::vector<MatrixXcd> H1b;
    std::vector<MatrixXcd> H2b;
    std::vector<MatrixXcd> H3b;
    for (int i=0;i<N_panels;i++) {
        H1a.push_back(MatrixXcd(6,6));
        H2a.push_back(MatrixXcd(6,6));
        H3a.push_back(MatrixXcd(6,6));
        H1b.push_back(MatrixXcd(6,6));
        H2b.push_back(MatrixXcd(6,6));
        H3b.push_back(MatrixXcd(6,6));
        H1a.at(i) = Hermitian(Phase_Choice_0, Em.at(i));
        H2a.at(i) = Hermitian(Phase_Choice_1, Em.at(i));
        H3a.at(i) = Hermitian(Phase_Choice_2, Em.at(i));
        H1b.at(i) = Hermitian(-Phase_Choice_0, Em.at(i));
        H2b.at(i) = Hermitian(-Phase_Choice_1, Em.at(i));
        H3b.at(i) = Hermitian(-Phase_Choice_2, Em.at(i));
        //std::cout << "Eigenvalues of panel " << i << " H matrix:" << std::endl;
        //std::cout << "Panel " << i << " H(pi) matrix:" << std::endl;
        //std::cout << H1.at(i) << std::endl;
        ComplexEigenSolver<MatrixXcd> eigensolver_1a(H1a.at(i));
        ComplexEigenSolver<MatrixXcd> eigensolver_1b(H1b.at(i));
        //std::cout << "The eigenvalues of H(pi) are:\n" << eigensolver_1.eigenvalues() << std::endl;
        //std::cout << "Panel " << i << " H(3pi/4) matrix:" << std::endl;
        //std::cout << H2.at(i) << std::endl;
        ComplexEigenSolver<MatrixXcd> eigensolver_2a(H2a.at(i));
        ComplexEigenSolver<MatrixXcd> eigensolver_2b(H2b.at(i));
        //std::cout << "The eigenvalues of H(3pi/4) are:\n" << eigensolver_2.eigenvalues() << std::endl;
        //std::cout << "Panel " << i << " H(-3pi/4) matrix:" << std::endl;
        //std::cout << H3.at(i) << std::endl;
        ComplexEigenSolver<MatrixXcd> eigensolver_3a(H3a.at(i));
        ComplexEigenSolver<MatrixXcd> eigensolver_3b(H3b.at(i));
        //std::cout << "The eigenvalues of H(-3pi/4) are:\n" << eigensolver_3.eigenvalues() << std::endl;
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    //std::cout << "Here we compute U_avg." << std::endl;
    //MatrixXcd Uavg(6,6);
    //Uavg = BuildUmatrix(Havg1, Havg2, Havg3);
    //std::cout << "U_avg:" << std::endl;
    //std::cout << Uavg << std::endl;
    std::cout << "Here we compute U(i) for each panel." << std::endl;
    std::vector<MatrixXcd> U;
    for (int i=0;i<N_panels;i++) {
        U.push_back(MatrixXcd(6,6));
        U.at(i) = BuildUmatrix(H1a.at(i),H1b.at(i),H2a.at(i),H2b.at(i),H3a.at(i),H3b.at(i));
        //std::cout << "U(" << i << "):" << std::endl;
        //std::cout << U.at(i) << std::endl;
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    const std::complex<double> If(0.0, 1.0);
    MatrixXcd D(6,6);
    D << cos(Phase_Choice_0) + sin(Phase_Choice_0) * If, 0., 0., 0., 0., 0.,
         0., cos(-Phase_Choice_0) + sin(-Phase_Choice_0) * If, 0., 0., 0., 0.,
         0., 0., cos(Phase_Choice_1) + sin(Phase_Choice_1) * If, 0., 0., 0.,
         0., 0., 0., cos(-Phase_Choice_1) + sin(-Phase_Choice_1) * If, 0., 0.,
         0., 0., 0., 0., cos(Phase_Choice_2) + sin(Phase_Choice_2) * If, 0.,
         0., 0., 0., 0., 0., cos(-Phase_Choice_2) + sin(-Phase_Choice_2) * If;
    //std::cout << "The D matrix is:" << std::endl;
    //std::cout << D << std::endl;

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    MatrixXcd Et_avg(6,6);
    Et_avg = MatrixXcd::Zero(6,6);
    //Et_avg = Uavg*D*Uavg.inverse();
    //std::cout << "Here we compute Et_avg = U_avg*D*U_avg^{-1}:" << std::endl;
    //std::cout << Et_avg << std::endl;
    //ComplexEigenSolver<MatrixXcd> eigensolver_Et(Et_avg);
    //std::cout << "The eigenvalues of Et_avg are:\n" << eigensolver_Et.eigenvalues() << std::endl;

    std::cout << "Here we also compute Et(i) = U(i)*D*U(i)^{-1} for each panel." << std::endl;
    std::vector<MatrixXcd> Et;
    for (int i=0;i<N_panels;i++) {
        Et.push_back(MatrixXcd(6,6));
        Et.at(i) = U.at(i)*D*U.at(i).inverse();
        //Et.at(i) << 0.86444, 0.51951, 2.9642, 0.51623, -3.6621, 1.0762,
        //            2.4213, -2.789, -3.7638, -1.4315, 0.22454, -4.1223,
        //            1.6891, -0.46836, -1.4982, -0.55488, -0.10290, -1.0299,
        //            8.4037, -3.5034, -4.3467, -3.5891, -3.6306, -7.6545,
        //            2.3642, 2.6671e-2, 2.0832, 0.19756, -3.8864, 0.31863,
        //            -4.607, 3.0437, 6.7292, 2.8721, -1.8445, 6.0698;
        std::cout << "Et(" << i << "):" << std::endl;
        std::cout << Et.at(i) << std::endl;
        Et_avg += 1./double(N_panels)*Et.at(i);
    }
    VectorXd Et_eigenvalue_real(6);
    VectorXd Et_eigenvalue_imag(6);
    for (int i=0;i<N_panels;i++) {
        ComplexEigenSolver<MatrixXcd> eigensolver_Et(Et.at(i));
        Et_eigenvalue_real = eigensolver_Et.eigenvalues().real();
        Et_eigenvalue_imag = eigensolver_Et.eigenvalues().imag();
        std::cout << "The eigenvalues of Et("<<i<<") are:" << std::endl;
        for (int j=0;j<6;j++)
        {
            std::cout << pow(Et_eigenvalue_real(j)*Et_eigenvalue_real(j)+Et_eigenvalue_imag(j)*Et_eigenvalue_imag(j),0.5) << ", " << atan2(Et_eigenvalue_imag(j), Et_eigenvalue_real(j))/M_PI << "pi" << std::endl;
        }
        std::cout << "Et("<< i << ")*Em(" << i << ")^{-1}:" << std::endl;
        std::cout << (Et.at(i)*Em.at(i).inverse()) << std::endl;
    }
    std::cout << "Here we compute Et_avg = sum 1/8*Et(i):" << std::endl;
    std::cout << "Et_avg:" << std::endl;
    std::cout << Et_avg << std::endl;
    ComplexEigenSolver<MatrixXcd> eigensolver_Et(Et_avg);
    std::cout << "The eigenvalues of Et_avg are:\n" << eigensolver_Et.eigenvalues() << std::endl;

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Here we compute H_avg from Et_avg" << std::endl;
    MatrixXcd H1a_avg;
    MatrixXcd H2a_avg;
    MatrixXcd H3a_avg;
    MatrixXcd H1b_avg;
    MatrixXcd H2b_avg;
    MatrixXcd H3b_avg;
    H1a_avg = Hermitian(Phase_Choice_0, Et_avg);
    H2a_avg = Hermitian(Phase_Choice_1, Et_avg);
    H3a_avg = Hermitian(Phase_Choice_2, Et_avg);
    H1b_avg = Hermitian(-Phase_Choice_0, Et_avg);
    H2b_avg = Hermitian(-Phase_Choice_1, Et_avg);
    H3b_avg = Hermitian(-Phase_Choice_2, Et_avg);
    std::cout << "Here we compute U_avg from H_avg." << std::endl;
    MatrixXcd U_avg;
    U_avg = BuildUmatrix(H1a_avg,H1b_avg,H2a_avg,H2b_avg,H3a_avg,H3b_avg);

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Here we compute Et_adp = U_avg*D*U_avg^{-1}. And we show that Et_adp*Em(i)^{-1} is close to unity." << std::endl;
    MatrixXcd Et_adp;
    Et_adp = U_avg*D*U_avg.inverse();
    //std::cout << "Et_adp:" << std::endl;
    //std::cout << Et_adp << std::endl;
    for (int i=0;i<N_panels;i++) {
        std::cout << "Et_adp*Em(" << i << ")^{-1}:" << std::endl;
        std::cout << (Et_adp*Em.at(i).inverse()) << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    //std::cout << "Here we compute Et_final = Et_adp + Et(i)-Et_avg" << std::endl;
    std::cout << "Here we compute Et_final = Et_adp" << std::endl;
    std::vector<MatrixXcd> Et_final;
    std::vector<MatrixXcd> EG;
    for (int i=0;i<N_panels;i++) {
        Et_final.push_back(MatrixXcd(6,6));
        EG.push_back(MatrixXcd(6,6));
        Et_final.at(i) = Et_adp;
        EG.at(i) = Em.at(i)-Et_adp;
        //std::cout << "Et_final(" << i << "):" << std::endl;
        //std::cout << Et_final.at(i) << std::endl;
        //std::cout << "Et_final("<< i << ")*Em(" << i << ")^{-1}:" << std::endl;
        //std::cout << Et_final.at(i)*Em.at(i).inverse() << std::endl;
        //ComplexEigenSolver<MatrixXcd> eigensolver_Et_final(Et_final.at(i));
        //std::cout << "The eigenvalues of Et_final(i) = Et_adp are:\n" << eigensolver_Et_final.eigenvalues() << std::endl;
        //ComplexEigenSolver<MatrixXcd> eigensolver_EG(EG.at(i));
        //std::cout << "The eigenvalues of EG(i) = Em(i)-Et_adp are:\n" << eigensolver_EG.eigenvalues() << std::endl;
    }
    VectorXd Et_final_eigenvalue_real(6);
    VectorXd Et_final_eigenvalue_imag(6);
    for (int i=0;i<N_panels;i++) {
        ComplexEigenSolver<MatrixXcd> eigensolver_Et_final(Et_final.at(i));
        Et_final_eigenvalue_real = eigensolver_Et_final.eigenvalues().real();
        Et_final_eigenvalue_imag = eigensolver_Et_final.eigenvalues().imag();
        std::cout << "The eigenvalues of Et_final("<<i<<") are:" << std::endl;
        for (int j=0;j<6;j++)
        {
            std::cout << pow(Et_final_eigenvalue_real(j)*Et_final_eigenvalue_real(j)+Et_final_eigenvalue_imag(j)*Et_final_eigenvalue_imag(j),0.5) << ", " << atan2(Et_final_eigenvalue_imag(j), Et_final_eigenvalue_real(j))/M_PI << "pi" << std::endl;
        }
    }


    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    std::cout << "Here we compute MRt(i) = Et_final(i)*Em(i)^{-1}*MR(i)." << std::endl;
    std::vector<MatrixXcd> MRt;
    for (int i=0;i<N_panels;i++) {
        MRt.push_back(MatrixXcd(6,6));
        MRt.at(i) = (Et_final.at(i)*Em.at(i).inverse())*MR.at(i);
        //std::cout << "MRt(" << i << ") = (Et_final(i)*Em(i)^{-1})*MR:" << std::endl;
        //std::cout << MRt.at(i) << std::endl;
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    MatrixXcd Q(N_panels*6,N_panels*6);
    MatrixXcd Q_reduced((N_panels-1)*6,(N_panels-1)*6);
    MatrixXcd G(N_panels*6,N_panels*6);
    MatrixXcd G_reduced((N_panels-1)*6,(N_panels-1)*6);
    MatrixXcd M_iterative((N_panels-1)*6,(N_panels-1)*6);
    MatrixXcd dV((N_panels-1)*6,1);
    MatrixXcd Q_inv_dV((N_panels-1)*6,1);
    std::vector<VectorXd> dS;
    VectorXd dS_final((N_panels-1)*6);
    VectorXd dS_8panels(N_panels*6);
    VectorXd noise_initial((N_panels-1)*6);
    VectorXd noise_final((N_panels-1)*6);
    VectorXd actuator(N_panels*6);

    G = MatrixXcd::Zero(N_panels*6,N_panels*6);
    for (int i=0;i<N_panels;i++) {
        G.block(i*6,i*6,6,6) = EG.at(i);
    }
    G_reduced = G.block(6,6,(N_panels-1)*6,(N_panels-1)*6);

    std::cout << "Here we compute Q matrix using Et_final(i)" << std::endl;
    Q = BuildQmatrix(Et_final,-1);
    //for (int i=0;i<N_panels;i++) {
        //std::cout << "Block of Q matrix (" << i << "," << i << "):" << std::endl;
        //std::cout << Q.block(i*6,i*6,6,6) << std::endl;
        //for (int j=0;j<N_panels;j++) {
        //    std::cout << "Block of Q matrix (" << i << "," << j << "):" << std::endl;
        //    std::cout << Q.block(i*6,j*6,6,6) << std::endl;
        //}
    //}
    ComplexEigenSolver<MatrixXcd> eigensolver_Q(Q);
    //std::cout << "The eigenvalues of 48x48 Q are:\n" << eigensolver_Q.eigenvalues() << std::endl;
    dV = BuilddVmatrix(Et_final,-1);
    Q_reduced = Q.block(6,6,(N_panels-1)*6,(N_panels-1)*6);
    ComplexEigenSolver<MatrixXcd> eigensolver_Q_reduced(Q_reduced);
    std::cout << "The eigenvalues of Q_reduced (after removing the first 6 columns and rows) are:\n" << eigensolver_Q_reduced.eigenvalues() << std::endl;
    //for (int i=0;i<42;i++) {
    //    std::cout << "Eigenvector of Q_reduced, eigenvalue = " << eigensolver_Q_reduced.eigenvalues()(i) << std::endl;
    //    std::cout << eigensolver_Q_reduced.eigenvectors().col(i) << std::endl;
    //}

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Here we compute dS vector using dS = sum e(i)*dV^{T}*e(i)/lambda(i)" << std::endl;
    dS.push_back(VectorXd((N_panels-1)*6));
    dS.at(0) = BuilddSmatrix(Q_reduced, dV, -1, (N_panels-1)*6);

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    int n_iteration = 100;
    M_iterative = Q_reduced.inverse()*G_reduced;
    ComplexEigenSolver<MatrixXcd> eigensolver_M_iterative(M_iterative);
    std::cout << "The eigenvalues of M_iterative (after removing the first 6 columns and rows) are:\n" << eigensolver_M_iterative.eigenvalues() << std::endl;
    for (int i=0;i<n_iteration;i++) {
        dS.push_back(VectorXd((N_panels-1)*6));
        dS.at(i+1) = BuilddSmatrixIterative(M_iterative, dS.at(i), -1, (N_panels-1)*6);
    }
    //for (int j=0;j<N_panels-1;j++) {
    //    for (int i=0;i<n_iteration;i++) {
    //        std::cout << i << "th iteration, dS(" << j << "):" << std::endl;
    //        std::cout << dS.at(i).segment(j*6,6).transpose() << std::endl;
    //    }
    //}
    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    dS_final = VectorXd::Zero((N_panels-1)*6);
    for (int i=0;i<n_iteration;i++) {
        dS_final += dS.at(i);
    }
    dS_8panels = VectorXd::Zero(N_panels*6);
    dS_8panels.segment(6,(N_panels-1)*6) = dS.at(0);
    //dS_8panels.segment(6,42) = dS_final;
    
    for (int i=0;i<N_panels;i++) {
        Em.at(i) = MR.at(i)*ML.at(i).inverse();
    }
    noise_initial = BuildNoiseVector(dS_8panels, Em); 

    for (int i=0;i<N_panels;i++) {
        std::cout << "dS_initial(" << i << "):" << std::endl;
        std::cout << dS_8panels.segment(i*6,6).transpose() << std::endl;
        std::cout << "Noise_initial = Em(i)*dS_initial(i) + dS_initial(i+1) + dSigma(i) (i=" << i << "):" << std::endl;
        std::cout << noise_initial.segment(i*6,6).transpose() << std::endl;
        std::cout << "dL(" << i << ") = ML(i).inverse()*dS(i) (total):" << std::endl;
        std::cout << ML.at(i).inverse()*dS_8panels.segment(i*6,6) << std::endl;
        //if (i!=0) {
        //std::cout << "dL(" << i << ") = ML(i).inverse()*dS(i) (each step):" << std::endl;
        //std::cout << 1./(double(i))*ML.at(i).inverse()*dS_8panels.segment(i*6,6) << std::endl;
        //}
        //actuator = ActuatorVector(i);
        //Std::cout << "old L(" << i << "):" << std::endl;
        //Std::cout << actuator << std::endl;
        //std::cout << "new L(" << i << ") = old L(i) + dL(i) (total):" << std::endl;
        //std::cout << actuator+ML.at(i).inverse()*dS_8panels.segment(i*6,6) << std::endl;
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Here we compute dS vector using dS = -Q^{-1}*dV, dimention = 42" << std::endl;
    Q_inv_dV = Q_reduced.inverse()*dV;
    for (int i=0;i<N_panels-1;i++) {
        std::cout << "Q_reduced^{-1}*dV  = -dS(" << i+1 << "):" << std::endl;
        std::cout << Q_inv_dV.block(i*6,0,6,1) << std::endl;
    }

    std::ofstream myfile;
    myfile.open ("Q_eigenvector.txt");
    for (int i=0;i<(N_panels-1)*6;i++) {
        for (int j=0;j<(N_panels-1)*6;j++) {
            myfile << eigensolver_Q_reduced.eigenvectors().col(i)(j).real() << "\n";
            //myfile << eigensolver_M_iterative.eigenvectors().col(i)(j).real() << "\n";
        }
    }
    myfile.close();



}
