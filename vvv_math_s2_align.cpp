
#include <iostream>
#include <math.h>
#include <complex>
#include "string"
#include "fstream"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;
//using Eigen::MatrixXd;

int N_panels = 8;

MatrixXcd fillMatrix(std::string name)
{
    MatrixXcd matrix(6,6);
    std::ifstream file;
    file.open ("S1_response_matrix.txt");
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
    MatrixXcd Q(48,48);
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
        if (i!=7) {
            Q.block(i*6,i*6,6,6) = I+E.at(j).transpose()*E.at(j);
            Q.block(i*6,(i+1)*6,6,6) = E.at(j).transpose();
            Q.block((i+1)*6,i*6,6,6) = E.at(j);
        }
        else {
            Q.block(7*6,7*6,6,6) = I+E.at(j).transpose()*E.at(j);
            Q.block(7*6,0,6,6) = E.at(j).transpose();
            Q.block(0,7*6,6,6) = E.at(j);
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
VectorXd BuildNoiseVector(VectorXd dS, std::vector<MatrixXcd> Et) 
{
    int n_panels = 8;
    VectorXd noise(48);
    noise = VectorXd::Zero(48);
    VectorXd dSigma(48);
    dSigma = VectorXd::Zero(48);
    dSigma.segment(42,6) = BigDeltaVector();
    for (int i=0;i<n_panels;i++) {
        if (i<n_panels-1) {
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
    MatrixXcd W(48,48);
    MatrixXcd W_reduced(42,48);
    MatrixXcd dV(42,1);
    MatrixXcd dSigma(48,1);
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
        if (i!=7) {
            W.block(i*6,i*6,6,6) = E.at(j).transpose();
            W.block((i+1)*6,i*6,6,6) = I;
        }
        else {
            W.block(7*6,7*6,6,6) = E.at(j).transpose();
            W.block(0,7*6,6,6) = I;
        }
    }
    W_reduced = W.block(6,0,42,48);
    dSigma = MatrixXcd::Zero(48,1);
    dSigma.block(42,0,6,1) = BigDeltaVector();
    std::cout << "Done calculating W and dSigma matrix." << std::endl;
    dV = W_reduced*dSigma;
    std::cout << "Done calculating dV matrix." << std::endl;
    return dV;
}
MatrixXcd BuildUmatrix(MatrixXcd H1, MatrixXcd H2, MatrixXcd H3)
{
    MatrixXcd U(6,6);

    ComplexEigenSolver<MatrixXcd> eigensolver1(H1);
    ComplexEigenSolver<MatrixXcd> eigensolver2(H2);
    ComplexEigenSolver<MatrixXcd> eigensolver3(H3);
    
    U.col(0) = eigensolver1.eigenvectors().col(0);
    U.col(1) = eigensolver1.eigenvectors().col(1);
    U.col(2) = eigensolver2.eigenvectors().col(0);
    U.col(3) = eigensolver3.eigenvectors().col(0);
    U.col(4) = eigensolver2.eigenvectors().col(1);
    U.col(5) = eigensolver3.eigenvectors().col(1);

    return U;
}
int main()
{
    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    //std::vector<MatrixXcd> MR;
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(0) = fillMatrix("M1R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(1) = fillMatrix("M2R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(2) = fillMatrix("M3R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(3) = fillMatrix("M4R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(4) = fillMatrix("M5R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(5) = fillMatrix("M6R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(6) = fillMatrix("M7R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(7) = fillMatrix("M8R");

    //std::vector<MatrixXcd> ML;
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(0) = fillMatrix("M1L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(1) = fillMatrix("M2L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(2) = fillMatrix("M3L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(3) = fillMatrix("M4L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(4) = fillMatrix("M5L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(5) = fillMatrix("M6L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(6) = fillMatrix("M7L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(7) = fillMatrix("M8L");

    //std::vector<MatrixXcd> MR;
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(0) = fillMatrix("M1R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(1) = fillMatrix("M2R");
    //MR.push_back(MatrixXcd(6,6));
    //MR.at(2) = fillMatrix("M5R");

    //std::vector<MatrixXcd> ML;
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(0) = fillMatrix("M1L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(1) = fillMatrix("M2L");
    //ML.push_back(MatrixXcd(6,6));
    //ML.at(2) = fillMatrix("M5L");

    std::vector<MatrixXcd> MR;
    MR.push_back(MatrixXcd(6,6));
    MR.at(0) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(1) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(2) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(3) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(4) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(5) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(6) = fillMatrix("M2R");
    MR.push_back(MatrixXcd(6,6));
    MR.at(7) = fillMatrix("M2R");

    std::vector<MatrixXcd> ML;
    ML.push_back(MatrixXcd(6,6));
    ML.at(0) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(1) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(2) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(3) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(4) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(5) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(6) = fillMatrix("M2L");
    ML.push_back(MatrixXcd(6,6));
    ML.at(7) = fillMatrix("M2L");

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    std::cout << "Here we compute Em(i) = MR(i)*ML(i)^{-1}" << std::endl;
    std::vector<MatrixXcd> Em;
    for (int i=0;i<N_panels;i++) {
        Em.push_back(MatrixXcd(6,6));
        Em.at(i) = MR.at(i)*ML.at(i).inverse();
        //std::cout << "Em(" << i << "):" << std::endl;
        //std::cout << Em.at(i) << std::endl;
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    //std::cout << "Here we compute H_avg = 1/8 * sum (H(i))" << std::endl;
    //MatrixXcd Havg1(6,6);
    //MatrixXcd Havg2(6,6);
    //MatrixXcd Havg3(6,6);
    //Havg1 = AvgHermitian(M_PI, Em);
    //Havg2 = AvgHermitian(3.*M_PI/4., Em);
    //Havg3 = AvgHermitian(-3.*M_PI/4., Em);
    //ComplexEigenSolver<MatrixXcd> eigensolver1(Havg1);
    //std::cout << "The eigenvalues of avg H(pi) are:\n" << eigensolver1.eigenvalues() << std::endl;
    //ComplexEigenSolver<MatrixXcd> eigensolver2(Havg2);
    //std::cout << "The eigenvalues of avg H(3pi/4) are:\n" << eigensolver2.eigenvalues() << std::endl;
    //ComplexEigenSolver<MatrixXcd> eigensolver3(Havg3);
    //std::cout << "The eigenvalues of avg H(-3pi/4) are:\n" << eigensolver3.eigenvalues() << std::endl;

    std::vector<MatrixXcd> H1;
    std::vector<MatrixXcd> H2;
    std::vector<MatrixXcd> H3;
    for (int i=0;i<N_panels;i++) {
        H1.push_back(MatrixXcd(6,6));
        H2.push_back(MatrixXcd(6,6));
        H3.push_back(MatrixXcd(6,6));
        H1.at(i) = Hermitian(M_PI, Em.at(i));
        H2.at(i) = Hermitian(3.*M_PI/4., Em.at(i));
        H3.at(i) = Hermitian(-3.*M_PI/4., Em.at(i));
        //std::cout << "Eigenvalues of panel " << i << " H matrix:" << std::endl;
        //std::cout << "Panel " << i << " H(pi) matrix:" << std::endl;
        //std::cout << H1.at(i) << std::endl;
        ComplexEigenSolver<MatrixXcd> eigensolver_1(H1.at(i));
        //std::cout << "The eigenvalues of H(pi) are:\n" << eigensolver_1.eigenvalues() << std::endl;
        //std::cout << "Panel " << i << " H(3pi/4) matrix:" << std::endl;
        //std::cout << H2.at(i) << std::endl;
        ComplexEigenSolver<MatrixXcd> eigensolver_2(H2.at(i));
        //std::cout << "The eigenvalues of H(3pi/4) are:\n" << eigensolver_2.eigenvalues() << std::endl;
        //std::cout << "Panel " << i << " H(-3pi/4) matrix:" << std::endl;
        //std::cout << H3.at(i) << std::endl;
        ComplexEigenSolver<MatrixXcd> eigensolver_3(H3.at(i));
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
        U.at(i) = BuildUmatrix(H1.at(i), H2.at(i), H3.at(i));
        //std::cout << "U(" << i << "):" << std::endl;
        //std::cout << U.at(i) << std::endl;
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

    const std::complex<double> If(0.0, 1.0);
    MatrixXcd D(6,6);
    D << -1.0 + 0.0 * If, 0., 0., 0., 0., 0.,
         0., -1.0 + 0.0 * If, 0., 0., 0., 0.,
         0., 0., cos(3*M_PI/4.) + sin(3*M_PI/4.) * If, 0., 0., 0.,
         0., 0., 0., cos(-3*M_PI/4.) + sin(-3*M_PI/4.) * If, 0., 0.,
         0., 0., 0., 0., cos(3*M_PI/4.) + sin(3*M_PI/4.) * If, 0.,
         0., 0., 0., 0., 0., cos(-3*M_PI/4.) + sin(-3*M_PI/4.) * If;
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
    std::cout << "Here we compute Et_avg = sum 1/8*Et(i):" << std::endl;
    std::cout << "Et_avg:" << std::endl;
    std::cout << Et_avg << std::endl;
    ComplexEigenSolver<MatrixXcd> eigensolver_Et(Et_avg);
    std::cout << "The eigenvalues of Et_avg are:\n" << eigensolver_Et.eigenvalues() << std::endl;

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Here we compute H_avg from Et_avg" << std::endl;
    MatrixXcd H1_avg;
    MatrixXcd H2_avg;
    MatrixXcd H3_avg;
    H1_avg = Hermitian(M_PI, Et_avg);
    H2_avg = Hermitian(3.*M_PI/4., Et_avg);
    H3_avg = Hermitian(-3.*M_PI/4., Et_avg);
    std::cout << "Here we compute U_avg from H_avg." << std::endl;
    MatrixXcd U_avg;
    U_avg = BuildUmatrix(H1_avg, H2_avg, H3_avg);

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Here we compute Et_adp = U_avg*D*U_avg^{-1}. And we show that Et_adp*Em(i)^{-1} is close to unity." << std::endl;
    MatrixXcd Et_adp;
    Et_adp = U_avg*D*U_avg.inverse();
    //std::cout << "Et_adp:" << std::endl;
    //std::cout << Et_adp << std::endl;
    //for (int i=0;i<N_panels;i++) {
    //    std::cout << "Et_adp*Em(" << i << ")^{-1}:" << std::endl;
    //    std::cout << Et_adp*Em.at(i).inverse() << std::endl;
    //}
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
        //ComplexEigenSolver<MatrixXcd> eigensolver_EG(EG.at(i));
        //std::cout << "The eigenvalues of EG(i) = Em(i)-Et_adp are:\n" << eigensolver_EG.eigenvalues() << std::endl;
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
    MatrixXcd Q(48,48);
    MatrixXcd Q_reduced(42,42);
    MatrixXcd G(48,48);
    MatrixXcd G_reduced(42,42);
    MatrixXcd M_iterative(42,42);
    MatrixXcd dV(42,1);
    MatrixXcd Q_inv_dV(42,1);
    std::vector<VectorXd> dS;
    VectorXd dS_final(42);
    VectorXd dS_8panels(48);
    VectorXd noise_initial(42);
    VectorXd noise_final(42);

    G = MatrixXcd::Zero(48,48);
    for (int i=0;i<N_panels;i++) {
        G.block(i*6,i*6,6,6) = EG.at(i);
    }
    G_reduced = G.block(6,6,42,42);

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
    Q_reduced = Q.block(6,6,42,42);
    ComplexEigenSolver<MatrixXcd> eigensolver_Q_reduced(Q_reduced);
    std::cout << "The eigenvalues of Q_reduced (after removing the first 6 columns and rows) are:\n" << eigensolver_Q_reduced.eigenvalues() << std::endl;
    //for (int i=0;i<42;i++) {
    //    std::cout << "Eigenvector of Q_reduced, eigenvalue = " << eigensolver_Q_reduced.eigenvalues()(i) << std::endl;
    //    std::cout << eigensolver_Q_reduced.eigenvectors().col(i) << std::endl;
    //}

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Here we compute dS vector using dS = sum e(i)*dV^{T}*e(i)/lambda(i)" << std::endl;
    dS.push_back(VectorXd(42));
    dS.at(0) = BuilddSmatrix(Q_reduced, dV, -1, 42);

    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    int n_iteration = 100;
    M_iterative = Q_reduced.inverse()*G_reduced;
    ComplexEigenSolver<MatrixXcd> eigensolver_M_iterative(M_iterative);
    std::cout << "The eigenvalues of M_iterative (after removing the first 6 columns and rows) are:\n" << eigensolver_M_iterative.eigenvalues() << std::endl;
    for (int i=0;i<n_iteration;i++) {
        dS.push_back(VectorXd(42));
        dS.at(i+1) = BuilddSmatrixIterative(M_iterative, dS.at(i), -1, 42);
    }
    //for (int j=0;j<N_panels-1;j++) {
    //    for (int i=0;i<n_iteration;i++) {
    //        std::cout << i << "th iteration, dS(" << j << "):" << std::endl;
    //        std::cout << dS.at(i).segment(j*6,6).transpose() << std::endl;
    //    }
    //}
    std::cout << "--------------------------------------------------------------------------------------------------------------------" << std::endl;
    dS_final = VectorXd::Zero(42);
    for (int i=0;i<n_iteration;i++) {
        dS_final += dS.at(i);
    }
    dS_8panels = VectorXd::Zero(48);
    dS_8panels.segment(6,42) = dS.at(0);
    //dS_8panels.segment(6,42) = dS_final;
    noise_initial = BuildNoiseVector(dS_8panels, Em); 
    for (int i=0;i<N_panels;i++) {
        std::cout << "dS_initial(" << i << "):" << std::endl;
        std::cout << dS_8panels.segment(i*6,6).transpose() << std::endl;
        std::cout << "Noise_initial = Em(i)*dS_initial(i) + dS_initial(i+1) + dSigma(i) (i=" << i << "):" << std::endl;
        std::cout << noise_initial.segment(i*6,6).transpose() << std::endl;
        std::cout << "dL(" << i << ") = ML(i).inverse()*dS(i):" << std::endl;
        std::cout << ML.at(i).inverse()*dS_8panels.segment(i*6,6) << std::endl;
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
    for (int i=0;i<42;i++) {
        for (int j=0;j<42;j++) {
            myfile << eigensolver_Q_reduced.eigenvectors().col(i)(j).real() << "\n";
            //myfile << eigensolver_M_iterative.eigenvectors().col(i)(j).real() << "\n";
        }
    }
    myfile.close();



}
