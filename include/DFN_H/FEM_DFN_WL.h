#pragma once
#include "../Math_WL_H/Math_WL.h"
#include "Eigen/Dense"
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{
class FEM_DFN
{
public:
    std::vector<std::vector<std::pair<Vector3d, Vector3d>>> Flux_rate;
    std::vector<std::vector<std::pair<Vector3d, Vector3d>>> Flux_rate_2D;
    ///<store the coordinates (x,y,z) that represent an element (usually it is the center of the element), and the flux rate vector at this point

    //std::vector<std::vector<RowVector6d>> pd_N_pd_x;
    // also the gridient

    //std::vector<std::vector<RowVector6d>> pd_N_pd_y;

public:
    FEM_DFN(DFN::DFN_mesh DFN_mesh, DFN::Domain dom);
    void Assemble_matrix(DFN::DFN_mesh DFN_mesh, const double Kper, double *K_overall, double *F_overall);
    void Matlab_FEM_head_result(string FileKey, DFN::DFN_mesh DFN_mesh, double *X_overall, DFN::Domain dom);
    void Calculate_flux_rate(DFN::DFN_mesh DFN_mesh, DFN::Domain dom, double *X_overall);
    void Verify_in_and_out_flux(DFN::DFN_mesh DFN_mesh, double *X_overall);
};

//*******************************************
inline FEM_DFN::FEM_DFN(DFN::DFN_mesh DFN_mesh, DFN::Domain dom)
{
    double k_coe = 1.0; // permeability of each fracture, relating to aperture
    double *K_overall = new double[DFN_mesh.Matrix_dimesions * DFN_mesh.Matrix_dimesions];
    if (K_overall == NULL)
    {
        cout << "Error! Cannot alloc to matrix 'K_overall'!\n";
        exit(0);
    }
    else
    {
        for (size_t i = 0; i < DFN_mesh.Matrix_dimesions * DFN_mesh.Matrix_dimesions; ++i)
        {
            K_overall[i] = 0;
        }
    }

    double *F_overall = new double[DFN_mesh.Matrix_dimesions];
    if (F_overall == NULL)
    {
        cout << "Error! Cannot alloc to matrix 'F_overall'!\n";
        exit(0);
    }
    else
    {
        for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
        {
            F_overall[i] = 0;
        }
    }

    Assemble_matrix(DFN_mesh, k_coe, K_overall, F_overall);

    /*
    cout << "\nK_overall;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        for (size_t j = 0; j < DFN_mesh.Matrix_dimesions; ++j)
        {
            size_t Idx = i * DFN_mesh.Matrix_dimesions + j;

            cout << K_overall[Idx] << "\t";
            size_t Idy = j * DFN_mesh.Matrix_dimesions + i;
            if (abs(K_overall[Idx] - K_overall[Idy]) > 0.0001)
            {
                cout << "Error! \n";
                exit(0);
            }
        }
        cout << "|" << endl;
    }*/

    /*
    cout << "\nF_overall;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        cout << F_overall[i] << endl;
    }*/

    // K x = F
    int n = DFN_mesh.Matrix_dimesions; // dimensions of coefficient matrix K (2D)
    int m = 1;                         // dimensions of F matrix
    //double a = 1, b = -1.00;

    //dgemv_("N", &n, &n, &a, K_overall, &n, F_overall, &m, &b, X_overall, &m);
    int *ipiv = new int[n];
    if (ipiv == NULL)
    {
        cout << "Error! Cannot alloc to matrix 'F_overall'!\n";
        exit(0);
    }
    else
    {
        for (size_t i = 0; i < (size_t)n; ++i)
        {
            ipiv[i] = 0;
        }
    }
    int info;
    std::cout << "start solving matrix;\n";
    dgesv_(&n, &m, K_overall, &n, ipiv, F_overall, &n, &info);
    std::cout << "finish solving matrix;\n";
    // now F_overall is the solution X
    if (info > 0)
    {
        cout << "The solution could not be computed!!\n";
        exit(0);
    }

    /*cout << "\nn " << n << endl;
    cout << "a " << a << endl;
    cout << "b " << b << endl;*/

    /*
    cout << "\nX_overall with BLAS;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        cout << X_overall[i] << endl;
    }
    */

    /*
    cout << "\n22 K_overall;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        for (size_t j = 0; j < DFN_mesh.Matrix_dimesions; ++j)
        {
            size_t Idx = i * DFN_mesh.Matrix_dimesions + j;

            cout << K_overall[Idx] << "\t";
            size_t Idy = j * DFN_mesh.Matrix_dimesions + i;
            if (abs(K_overall[Idx] - K_overall[Idy]) > 0.0001)
            {
                cout << "Error! \n";
                exit(0);
            }
        }
        cout << "|" << endl;
    }

    cout << "\n22 F_overall;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        cout << F_overall[i] << endl;
    }
    */
    //dgemv_("N", DFN_mesh.Matrix_dimesions, DFN_mesh.Matrix_dimesions, 1, K_overall, /*LDA*/ DFN_mesh.Matrix_dimesions, X_overall, DFN_mesh.Matrix_dimesions /*INCX*/, 1, F_overall, DFN_mesh.Matrix_dimesions);
    Calculate_flux_rate(DFN_mesh, dom, F_overall);
    Matlab_FEM_head_result("tdfn_color_head.m", DFN_mesh, F_overall, dom);
    Verify_in_and_out_flux(DFN_mesh, F_overall);
    delete[] ipiv;
    ipiv = NULL;
    delete[] F_overall;
    F_overall = NULL;
    delete[] K_overall;
    K_overall = NULL;
};
//*******************************************

inline void FEM_DFN::Assemble_matrix(DFN::DFN_mesh DFN_mesh, const double Kper, double *K_overall, double *F_overall)
{
    /*
    Eigen::MatrixXd K, F;
    K = Eigen::MatrixXd::Zero(DFN_mesh.Matrix_dimesions, DFN_mesh.Matrix_dimesions);
    F = Eigen::MatrixXd::Zero(DFN_mesh.Matrix_dimesions, 1);*/
    //pd_N_pd_x.resize(DFN_mesh.JM.size());
    //pd_N_pd_y.resize(DFN_mesh.JM.size());
    for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
    {
        //pd_N_pd_x[i].resize(DFN_mesh.JM[i].size()); // the fracture NO. i
        //pd_N_pd_y[i].resize(DFN_mesh.JM[i].size());
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            Eigen::RowVectorXd w, xi, eta;
            w = Eigen::RowVectorXd::Zero(6);
            w << 0.1713244923791700, 0.3607615730481380, 0.4679139345726910, 0.1713244923791700, 0.3607615730481380, 0.4679139345726910;
            xi = Eigen::RowVectorXd::Zero(6);
            xi << 0.9324695142031520, 0.6612093864662640, 0.2386191860831960, -0.9324695142031520, -0.6612093864662640, -0.2386191860831960;
            eta = xi;

            Eigen::MatrixXd Ke;
            Ke = Eigen::MatrixXd::Zero(6, 6);

            for (size_t ik = 0; ik < (size_t)Ke.rows(); ++ik)
            {
                for (size_t jk = 0; jk < (size_t)Ke.cols(); ++jk)
                {
                    Eigen::VectorXd pd_N_over_pd_xi, pd_N_over_pd_eta;
                    pd_N_over_pd_xi = Eigen::VectorXd::Zero(6);
                    pd_N_over_pd_eta = pd_N_over_pd_xi;

                    double N_natural_a_1 = 4 * xi(ik) + 4 * eta(jk) - 3;
                    double N_natural_a_2 = 4 * xi(ik) - 1;
                    double N_natural_a_3 = 0;
                    double N_natural_a_4 = 4 - 8 * xi(ik) - 4 * eta(jk);
                    double N_natural_a_5 = 4 * eta(jk);
                    double N_natural_a_6 = -4 * eta(jk);

                    //1 - 4 - 2 -5 -3 - 6
                    pd_N_over_pd_xi << N_natural_a_1, N_natural_a_4, N_natural_a_2, N_natural_a_5, N_natural_a_3, N_natural_a_6;

                    double N_natural_b_1 = 4 * xi(ik) + 4 * eta(jk) - 3;
                    double N_natural_b_2 = 0;
                    double N_natural_b_3 = 4 * eta(jk) - 1;
                    double N_natural_b_4 = -4 * xi(ik);
                    double N_natural_b_5 = 4 * xi(ik);
                    double N_natural_b_6 = 4 - 8 * eta(jk) - 4 * xi(ik);
                    //1 - 4 - 2 -5 -3 - 6
                    pd_N_over_pd_eta << N_natural_b_1, N_natural_b_4, N_natural_b_2, N_natural_b_5, N_natural_b_3, N_natural_b_6;

                    Eigen::VectorXd pd_x_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_xi, pd_y_over_pd_eta;
                    pd_x_over_pd_xi = Eigen::VectorXd::Zero(6);
                    pd_x_over_pd_eta = pd_x_over_pd_xi;
                    pd_y_over_pd_xi = pd_x_over_pd_xi;
                    pd_y_over_pd_eta = pd_x_over_pd_xi;

                    Eigen::MatrixXd JXYe_x, JXYe_y;
                    JXYe_x = Eigen::MatrixXd::Zero(6, 1);
                    JXYe_y = Eigen::MatrixXd::Zero(6, 1);

                    for (size_t kk = 0; kk < (size_t)JXYe_x.rows(); ++kk)
                    {
                        JXYe_x(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](kk)](0);
                        JXYe_y(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](kk)](1);
                    }

                    pd_x_over_pd_xi = pd_N_over_pd_xi.transpose() * JXYe_x;
                    pd_x_over_pd_eta = pd_N_over_pd_eta.transpose() * JXYe_x;

                    pd_y_over_pd_xi = pd_N_over_pd_xi.transpose() * JXYe_y;
                    pd_y_over_pd_eta = pd_N_over_pd_eta.transpose() * JXYe_y;

                    Eigen::MatrixXd Jacobi;
                    Jacobi = Eigen::MatrixXd::Zero(2, 2);
                    Jacobi << pd_x_over_pd_xi, pd_y_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_eta;

                    Eigen::MatrixXd tem_o;
                    tem_o = Eigen::MatrixXd::Zero(2, pd_N_over_pd_xi.rows());
                    for (size_t op = 0; op < (size_t)tem_o.cols(); ++op)
                    {
                        tem_o(0, op) = pd_N_over_pd_xi(op);
                        tem_o(1, op) = pd_N_over_pd_eta(op);
                    }

                    Eigen::MatrixXd AAA;
                    AAA = Jacobi.inverse() * tem_o;

                    Eigen::VectorXd pd_N_over_pd_x, pd_N_over_pd_y;
                    pd_N_over_pd_x = Eigen::VectorXd::Zero(pd_N_over_pd_xi.rows());
                    pd_N_over_pd_y = Eigen::VectorXd::Zero(pd_N_over_pd_xi.rows());
                    for (size_t op = 0; op < (size_t)pd_N_over_pd_x.rows(); ++op)
                    {
                        pd_N_over_pd_x(op) = AAA(0, op);
                        pd_N_over_pd_y(op) = AAA(1, op);
                    }
                    /*
                    for (size_t op = 0; op < (size_t)pd_N_pd_x[i][j].rows(); ++op)
                    {
                        pd_N_pd_x[i][j](op, 0) = pd_N_over_pd_x(op);
                        pd_N_pd_y[i][j](op, 0) = pd_N_over_pd_y(op);
                    }*/
                    Ke = Ke + Kper * w(ik) * w(jk) * Jacobi.determinant() * (pd_N_over_pd_x * pd_N_over_pd_x.transpose() + pd_N_over_pd_y * pd_N_over_pd_y.transpose());
                };
            };

            //cout << "element NO " << j << endl;
            for (size_t m = 0; m < (size_t)Ke.rows(); ++m)
            {
                for (size_t n = 0; n < (size_t)Ke.cols(); ++n)
                {
                    //i is the fracture ID
                    //j is the element ID
                    // both m and n  are node IDs
                    //cout << Ke(m, n) << "\t";
                    size_t Node_m = DFN_mesh.JM[i][j](m);
                    size_t Node_n = DFN_mesh.JM[i][j](n);
                    if (DFN_mesh.Coe_Matr_guide[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](m)](0);
                    }

                    if (DFN_mesh.Coe_Matr_guide[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](n)](0);
                    }

                    size_t Idx = Node_m * DFN_mesh.Matrix_dimesions + Node_n;
                    K_overall[Idx] = K_overall[Idx] + Ke(m, n);

                    //K(Node_m, Node_n) = K(Node_m, Node_n) + Ke(m, n);
                    /*
                    cout << "ele: " << j << ", "
                         << "local-> , m = " << m << ", n = " << n << endl;
                    cout << "ele: " << j << ", "
                         << "global->, i = " << DFN_mesh.JM[i][j](m) << ", j = " << DFN_mesh.JM[i][j](n) << endl;
                    cout << "Idx = " << Idx << endl;
                    cout << endl;

                    if (K_overall[Idx] != K(Node_m, Node_n))
                        exit(0);*/
                }
                //cout << endl;
            }
            //cout << endl;
        }
        /*
        cout << "LKKKKK: \n";
        for (size_t i = 0; i < K.rows(); ++i)
        {
            for (size_t j = 0; j < K.cols(); ++j)
            {
                cout << K(i, j) << "\t";
            }
            cout << endl;
        }*/
        /*
        std::map<std::pair<size_t, size_t>, Vector5d>::iterator it_1;
        it_1 = DFN_mesh.JB_2[i].begin();
        while(it_1 != DFN_mesh.JB_2[i].end())
        {
            cout << "element: " << it_1->first.first << ", edge: " << it_1->first.second << endl;
            it_1++;
        }
        cout << endl;*/
        for (size_t ti = 0; ti < DFN_mesh.JB_2[i].size(); ++ti)
        {
            //std::vector<std::map<std::pair<size_t, size_t>, Vector5d>> DFN_mesh.JB_2;
            std::map<std::pair<size_t, size_t>, Vector5d>::iterator it;
            it = DFN_mesh.JB_2[i].begin();
            for (size_t j = 0; j < ti; ++j)
                it++;
            if (it == DFN_mesh.JB_2[i].end())
            {
                cout << "Error happened when tranversing 'DFN_mesh.JB_2'!\n";
                exit(0);
            }

            int II = it->first.first;                   //DFN_mesh.JB_2[i](ti, 0);
            int bondary_side_number = it->first.second; //DFN_mesh.JB_2(ti, 1);

            Eigen::VectorXd Fe;
            Fe = Eigen::VectorXd::Zero(6);
            Eigen::RowVectorXd Jx, Jy;
            Jx = Eigen::RowVectorXd::Zero(6);
            Jy = Jx;
            for (size_t j = 0; j < (size_t)Fe.size(); ++j)
            {
                Jx(j) = DFN_mesh.JXY[i][DFN_mesh.JM[i][II](j)](0);
                Jy(j) = DFN_mesh.JXY[i][DFN_mesh.JM[i][II](j)](1);
            }

            if (bondary_side_number == 0)
            {
                double q1 = it->second(2);
                double q2 = it->second(3);
                double q3 = it->second(4);
                double L = pow(pow(Jx(0) - Jx(2), 2) + pow(Jy(0) - Jy(2), 2), 0.5);
                L = L * 0.5;

                Fe(0) = Kper * (L * (4 * q1 + 2 * q2 - q3)) / 30;
                Fe(1) = Kper * (L * (q1 + 8 * q2 + q3)) / 15;
                Fe(2) = Kper * (L * (2 * q2 - q1 + 4 * q3)) / 30;
            }
            else if (bondary_side_number == 1)
            {
                double q1 = it->second(2);
                double q2 = it->second(3);
                double q3 = it->second(4);
                double L = pow(pow(Jx(2) - Jx(4), 2) + pow(Jy(2) - Jy(4), 2), 0.5);
                L = L * 0.5;

                Fe(2) = Kper * (L * (4 * q1 + 2 * q2 - q3)) / 30;
                Fe(3) = Kper * (L * (q1 + 8 * q2 + q3)) / 15;
                Fe(4) = Kper * (L * (2 * q2 - q1 + 4 * q3)) / 30;
            }
            else if (bondary_side_number == 2)
            {
                double q1 = it->second(2);
                double q2 = it->second(3);
                double q3 = it->second(4);
                double L = pow(pow(Jx(4) - Jx(0), 2) + pow(Jy(4) - Jy(0), 2), 0.5);
                L = L * 0.5;

                Fe(4) = Kper * (L * (4 * q1 + 2 * q2 - q3)) / 30;
                Fe(5) = Kper * (L * (q1 + 8 * q2 + q3)) / 15;
                Fe(0) = Kper * (L * (2 * q2 - q1 + 4 * q3)) / 30;
            }

            //cout << "Fe: " << ti << endl;
            for (size_t s = 0; s < 6; ++s)
            {
                size_t Node_s = DFN_mesh.JM[i][II](s);
                if (DFN_mesh.Coe_Matr_guide[i][Node_s](0) == -1)
                {
                    size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_s](1); // the i th frac
                    size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_s](2); // the j th node
                    Node_s = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                }
                else
                {
                    Node_s = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][II](s)](0);
                }
                //cout << Fe(s, 0) << endl;
                F_overall[Node_s] = F_overall[Node_s] + Fe(s);
                //F(JM(II - 1, s) - 1, 0) = F(JM(II - 1, s) - 1, 0) + Fe(s, 0);
                //F(Node_s, 0) = F(Node_s, 0) + Fe(s, 0);
            }
        }

        double dashu = 1e16; // large number
        std::map<size_t, double>::iterator it_s;
        it_s = DFN_mesh.JB_1[i].begin();
        while (it_s != DFN_mesh.JB_1[i].end())
        {
            size_t NODE_s = it_s->first;
            double BC_1 = it_s->second;

            size_t IDX;
            double Guide_1 = DFN_mesh.Coe_Matr_guide[i][NODE_s](0);
            if (Guide_1 == -1)
            {
                size_t frac_i = DFN_mesh.Coe_Matr_guide[i][NODE_s](1);
                size_t node_j = DFN_mesh.Coe_Matr_guide[i][NODE_s](2);
                IDX = DFN_mesh.Coe_Matr_guide[frac_i][node_j](0);
            }
            else
            {
                IDX = DFN_mesh.Coe_Matr_guide[i][NODE_s](0);
            }

            size_t IDX_A = IDX * DFN_mesh.Matrix_dimesions + IDX;
            K_overall[IDX_A] = K_overall[IDX_A] * dashu;
            F_overall[IDX] = BC_1 * K_overall[IDX_A];
            it_s++;

            //K(IDX, IDX) = K(IDX, IDX) * dashu;
            //F(IDX) = BC_1 * K(IDX, IDX);
        }
    };
    /*
    Eigen::VectorXd head;
    head = Eigen::VectorXd::Zero(F.size());
    head = K.inverse() * F;
    cout << "X with Eigen:\n";
    cout << head << endl;*/

    ///----should be removed
    /*
    Eigen::MatrixXd K_inverse;
    K_inverse = K.inverse();

    for (size_t i = 0; i < K_inverse.rows(); ++i)
    {
        for (size_t j = 0; j < K_inverse.cols(); ++j)
        {
            size_t Idx = i * DFN_mesh.Matrix_dimesions + j;
            K_overall[Idx] = K_inverse(i, j);
        }
    }
    */
    //inverse_BLAS(K_overall, DFN_mesh.Matrix_dimesions);
};

inline void FEM_DFN::Matlab_FEM_head_result(string FileKey, DFN::DFN_mesh DFN_mesh, double *X_overall, DFN::Domain dom)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n\n";

    size_t NO_frac = DFN_mesh.JXY.size();
    for (size_t i = 0; i < NO_frac; ++i)
    {
        oss << "JXY_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JXY_3D[i].size(); ++j)
        {
            oss << round(DFN_mesh.JXY_3D[i][j](0), 3) << ", ";
            oss << round(DFN_mesh.JXY_3D[i][j](1), 3) << ", ";
            oss << round(DFN_mesh.JXY_3D[i][j](2), 3) << "; ";
        }
        oss << "];\n";

        oss << "JXY_2D_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JXY_3D[i].size(); ++j)
        {
            oss << round(DFN_mesh.JXY[i][j](0), 3) << ", ";
            oss << round(DFN_mesh.JXY[i][j](1), 3) << ", ";
            oss << round(DFN_mesh.JXY[i][j](2), 3) << "; ";
        }
        oss << "];\n";

        oss << "JM_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            oss << DFN_mesh.JM[i][j](0) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](1) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](2) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](3) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](4) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](5) + 1 << "; ";
        }
        oss << "];\n";

        oss << "Data_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.Coe_Matr_guide[i].size(); ++j)
        {
            size_t Node_s = 0;
            if (DFN_mesh.Coe_Matr_guide[i][j](0) == -1)
            {
                size_t i_frac = DFN_mesh.Coe_Matr_guide[i][j](1);
                size_t j_node = DFN_mesh.Coe_Matr_guide[i][j](2);
                Node_s = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
            }
            else
            {
                Node_s = DFN_mesh.Coe_Matr_guide[i][j](0);
            }
            oss << round(X_overall[Node_s], 4) << "; ";
        }
        oss << "];\n";
        oss << "figure(1)\n";
        oss << "P_" << i + 1 << " = patch('Vertices', JXY_" << i + 1 << ", 'Faces', JM_" << i + 1 << ", 'FaceVertexCData', Data_" << i + 1 << ", 'FaceColor', 'interp', 'EdgeAlpha', 1);\n";
        oss << "hold on;\n";
    }
    oss << "view(3);\n";
    oss << "colorbar;\n";

    oss << "hold on;\n";
    //Plotting the model domain
    oss << "figure(1)\n";
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << dom.Surfaces[i].Verts[j](0) << " " << dom.Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](1) << " " << dom.Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](2) << " " << dom.Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on; hold on;\n";
        }
    }
    double xmin_1 = dom.Model_domain(4), xmax_1 = dom.Model_domain(5);
    double ymin_1 = dom.Model_domain(2), ymax_1 = dom.Model_domain(3);
    double zmin_1 = dom.Model_domain(1), zmax_1 = dom.Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN head distribution');\nhold on;\n";

    oss << endl;
    oss << "%%*****flux rate vector***\n";
    // std::vector<std::vector<std::pair<Vector3d, Vector3d>>> Flux_rate;
    oss << "figure(1)\n";
    for (size_t i = 0; i < Flux_rate.size(); ++i)
    {
        oss << "q_vector_" << i + 1 << "=[";
        for (size_t j = 0; j < Flux_rate[i].size(); ++j)
        {
            for (size_t yj = 0; yj < 3; ++yj)
                oss << Flux_rate[i][j].first(yj) << ", ";
            for (size_t yj = 0; yj < 3; ++yj)
            {
                oss << Flux_rate[i][j].second(yj);
                if (yj == 2)
                    oss << "; ";
                else
                    oss << ", ";
            }
        }
        oss << "];\n";
        oss << "quiver3(q_vector_" << i + 1 << "(:,1), q_vector_" << i + 1 << "(:,2), q_vector_" << i + 1 << "(:,3), q_vector_" << i + 1 << "(:,4), q_vector_" << i + 1 << "(:,5), q_vector_" << i + 1 << "(:,6));\n\n";
        oss << "hold on;\n";
    };
    oss << "%show frac pressure contour and flow rate vectors in 2D\n";

    for (size_t i = 0; i < NO_frac; ++i)
    {
        oss << "figure(" << i + 2 << ")\n";
        oss << "x_tr" << i + 1 << " = JXY_2D_" << i + 1 << "(:, 1);\n";
        oss << "y_tr" << i + 1 << " = JXY_2D_" << i + 1 << "(:, 2);\n";
        oss << "z_tr" << i + 1 << " = Data_" << i + 1 << "(:);\n";

        oss << "nx" << i + 1 << " = 500;\n";
        oss << "ny" << i + 1 << " = 500;\n";
        oss << "[X" << i + 1 << ",Y" << i + 1 << "]";
        oss << " = meshgrid(linspace(min(x_tr" << i + 1 << "),max(x_tr" << i + 1 << "),nx" << i + 1 << "),linspace(min(y_tr" << i + 1 << "),max(y_tr" << i + 1 << "),ny" << i + 1 << ")) ;\n";
        oss << "Z" << i + 1 << " =griddata(x_tr" << i + 1 << ",y_tr" << i + 1 << ",z_tr" << i + 1 << ",X" << i + 1 << ",Y" << i + 1 << ") ;\n";
        oss << "contour(X" << i + 1 << ",Y" << i + 1 << ",Z" << i + 1 << ") ;\n";
        oss << "hold on;\n";

        oss << "S_" << i + 1 << "= patch('Vertices', JXY_2D_" << i + 1 << ", 'Faces', JM_" << i + 1 << ", 'FaceVertexCData', Data_" << i + 1 << ", 'FaceColor', 'interp', 'EdgeAlpha', 0.2, 'facealpha', 0);\n";
        oss << "hold on;\n";
        oss << "q_vector_2d" << i + 1 << "=[";
        for (size_t j = 0; j < Flux_rate_2D[i].size(); ++j)
        {
            for (size_t yj = 0; yj < 3; ++yj)
                oss << Flux_rate_2D[i][j].first(yj) << ", ";
            for (size_t yj = 0; yj < 3; ++yj)
            {
                oss << Flux_rate_2D[i][j].second(yj);
                if (yj == 2)
                    oss << "; ";
                else
                    oss << ", ";
            }
        }
        oss << "];\n";
        oss << "quiver3(q_vector_2d" << i + 1 << "(:,1), q_vector_2d" << i + 1 << "(:,2), q_vector_2d" << i + 1 << "(:,3), q_vector_2d" << i + 1 << "(:,4), q_vector_2d" << i + 1 << "(:,5), q_vector_2d" << i + 1 << "(:,6));\n\n";
        oss << "hold on;\n";
        oss << "xlabel('x (m)');\nylabel('y (m)');\ntitle('head contour and flow rate vector (Fig. " << i + 1 << ")');\nhold on;\n";
        if (i == 0 && NO_frac != 1)
        {
            oss << "display('are you sure that you wanna show all fractures?');\n";
            oss << "pause;\n";
        }
    }
    oss.close();
};

inline void FEM_DFN::Calculate_flux_rate(DFN::DFN_mesh DFN_mesh, DFN::Domain dom, double *X_overall)
{
    if (Flux_rate.size() != 0 && Flux_rate_2D.size() != 0)
    {
        cout << "Error! The array 'Flux_rate' is not initialized!\n";
        exit(0);
    }
    Flux_rate.resize(DFN_mesh.JM.size());
    Flux_rate_2D.resize(DFN_mesh.JM.size());
    for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
    {
        Flux_rate[i].resize(DFN_mesh.JM[i].size());
        Flux_rate_2D[i].resize(DFN_mesh.JM[i].size());
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            double K_coe = 1; // permeability
            Vector3d Center;  ///< The centroid of a triangle, also known as 'center of gravity'
            Vector3d pnt0 = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](0)], pnt2 = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](2)], pnt4 = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](4)];
            Vector3d pnt1 = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](1)], pnt3 = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](3)], pnt5 = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](5)];

            Eigen::VectorXd h_e;
            h_e = Eigen::VectorXd::Zero(6);

            //----extract pressure values
            for (size_t yj = 0; yj < (size_t)h_e.size(); ++yj)
            {
                /*
                size_t kf = yj;
                if (yj == 1)
                    kf = 2;
                else if (yj == 2)
                    kf = 4;*/
                size_t node_ID = DFN_mesh.JM[i][j](yj);
                int Id_guide = DFN_mesh.Coe_Matr_guide[i][node_ID](0);
                if (Id_guide == -1)
                {
                    size_t i_frac = DFN_mesh.Coe_Matr_guide[i][node_ID](1);
                    size_t j_node = DFN_mesh.Coe_Matr_guide[i][node_ID](2);
                    Id_guide = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    if (Id_guide == -1)
                    {
                        cout << "Error! find wrong repetitive point!\n";
                        exit(0);
                    }
                    h_e[yj] = X_overall[Id_guide];
                }
                else
                {
                    h_e[yj] = X_overall[Id_guide];
                }
            }
            //----

            //Center = (pnt0 + pnt2 + pnt4) / 3.;
            double qx, qy;

            Eigen::RowVectorXd p_N_p_x, p_N_p_y;
            p_N_p_x = Eigen::RowVectorXd::Zero(6);
            //p_N_p_x = Eigen::RowVectorXd::Zero(3);
            p_N_p_y = p_N_p_x;

            //------------------------
            //
            /*
            double b1, b2, b3, c1, c2, c3, A_area, edg1, edg2, edg3, p;
            edg1 = (pnt0 - pnt2).norm();
            edg2 = (pnt2 - pnt4).norm();
            edg3 = (pnt4 - pnt0).norm();
            p = (edg1 + edg2 + edg3) / 2.;

            A_area = pow(p * (p - edg1) * (p - edg2) * (p - edg3), 0.5);
            b1 = 1. / (2. * A_area) * (pnt2(1) - pnt4(1));
            b2 = 1. / (2. * A_area) * (pnt4(1) - pnt0(1));
            b3 = 1. / (2. * A_area) * (pnt0(1) - pnt2(1));
            p_N_p_x << b1, b2, b3;

            c1 = 1. / (2. * A_area) * (pnt4(2) - pnt2(2));
            c2 = 1. / (2. * A_area) * (pnt0(2) - pnt4(2));
            c3 = 1. / (2. * A_area) * (pnt2(2) - pnt0(2));
            p_N_p_y << c1, c2, c3;
            */

            /*
            double x, y;
            double xN1 = pnt0(0), xN2 = pnt2(0), xN3 = pnt4(0), xN4 = pnt1(0), xN5 = pnt3(0), xN6 = pnt5(0);
            double yN1 = pnt0(1), yN2 = pnt2(1), yN3 = pnt4(1), yN4 = pnt1(1), yN5 = pnt3(1), yN6 = pnt5(1);
            x = Center(0);
            y = Center(1);
            double yN46 = yN4 - yN6;
            double yN23 = yN2 - yN3;
            double xN23 = xN2 - xN3;
            double xN13 = xN1 - xN3;
            double yN13 = yN1 - yN3;
            double xN16 = xN1 - xN6;
            double xN46 = xN4 - xN6;
            double yN16 = yN1 - yN6;
            double yN54 = yN5 - yN4;
            double yN31 = yN3 - yN1;
            double xN31 = xN3 - xN1;
            double yN21 = yN2 - yN1;
            double xN24 = xN2 - xN4;
            double xN21 = xN2 - xN1;
            double xN54 = xN5 - xN4;
            double yN24 = yN2 - yN4;
            double yN56 = yN5 - yN6;
            double xN36 = xN3 - xN6;
            double xN56 = xN5 - xN6;
            double yN36 = yN3 - yN6;
            double yN43 = yN4 - yN3;
            double xN43 = xN4 - xN3;
            double xN41 = xN4 - xN1;
            double yN41 = yN4 - yN1;
            double yN51 = yN5 - yN1;
            double xN51 = xN5 - xN1;
            double yN61 = yN6 - yN1;
            double xN61 = xN6 - xN1;
            double yN63 = yN6 - yN3;
            double xN63 = xN6 - xN3;
            p_N_p_x(0) = (yN46 * (yN23 * (x - xN3) - xN23 * (y - yN3))) / ((xN13 * yN23 - xN23 * yN13) * (xN16 * yN46 - xN46 * yN16)) + (yN23 * (yN46 * (x - xN6) - xN46 * (y - yN6))) / ((xN13 * yN23 - xN23 * yN13) * (xN16 * yN46 - xN46 * yN16));
            p_N_p_y(0) = -(xN46 * (yN23 * (x - xN3) - xN23 * (y - yN3))) / ((xN13 * yN23 - xN23 * yN13) * (xN16 * yN46 - xN46 * yN16)) - (xN23 * (yN46 * (x - xN6) - xN46 * (y - yN6))) / ((xN13 * yN23 - xN23 * yN13) * (xN16 * yN46 - xN46 * yN16));

            p_N_p_x(1) = (yN54 * (yN31 * (x - xN1) - xN31 * (y - yN1))) / ((xN21 * yN31 - xN31 * yN21) * (xN24 * yN54 - xN54 * yN24)) + (yN31 * (yN54 * (x - xN4) - xN54 * (y - yN4))) / ((xN21 * yN31 - xN31 * yN21) * (xN24 * yN54 - xN54 * yN24));
            p_N_p_y(1) = -(xN54 * (yN31 * (x - xN1) - xN31 * (y - yN1))) / ((xN21 * yN31 - xN31 * yN21) * (xN24 * yN54 - xN54 * yN24)) - (xN31 * (yN54 * (x - xN4) - xN54 * (y - yN4))) / ((xN21 * yN31 - xN31 * yN21) * (xN24 * yN54 - xN54 * yN24));

            p_N_p_x(2) = -(yN56 * (yN21 * (x - xN1) - xN21 * (y - yN1))) / ((xN21 * yN31 - xN31 * yN21) * (xN36 * yN56 - xN56 * yN36)) - (yN21 * (yN56 * (x - xN6) - xN56 * (y - yN6))) / ((xN21 * yN31 - xN31 * yN21) * (xN36 * yN56 - xN56 * yN36));
            p_N_p_y(2) = (xN56 * (yN21 * (x - xN1) - xN21 * (y - yN1))) / ((xN21 * yN31 - xN31 * yN21) * (xN36 * yN56 - xN56 * yN36)) + (xN21 * (yN56 * (x - xN6) - xN56 * (y - yN6))) / ((xN21 * yN31 - xN31 * yN21) * (xN36 * yN56 - xN56 * yN36));

            p_N_p_x(3) = -(yN31 * (yN23 * (x - xN3) - xN23 * (y - yN3))) / ((xN23 * yN43 - xN43 * yN23) * (xN41 + yN31 - xN31 * yN41)) - (yN23 * (yN31 * (x - xN1) - xN31 * (y - yN1))) / ((xN23 * yN43 - xN43 * yN23) * (xN41 + yN31 - xN31 * yN41));
            p_N_p_y(3) = (xN31 * (yN23 * (x - xN3) - xN23 * (y - yN3))) / ((xN23 * yN43 - xN43 * yN23) * (xN41 + yN31 - xN31 * yN41)) + (xN23 * (yN31 * (x - xN1) - xN31 * (y - yN1))) / ((xN23 * yN43 - xN43 * yN23) * (xN41 + yN31 - xN31 * yN41));

            p_N_p_x(4) = (yN54 * (yN31 * (x - xN1) - xN31 * (y - yN1))) / ((xN21 * yN51 - xN51 * yN21) * (xN31 * yN51 - xN51 * yN31)) + (yN31 * (yN54 * (x - xN4) - xN54 * (y - yN4))) / ((xN21 * yN51 - xN51 * yN21) * (xN31 * yN51 - xN51 * yN31));
            p_N_p_y(4) = -(xN54 * (yN31 * (x - xN1) - xN31 * (y - yN1))) / ((xN21 * yN51 - xN51 * yN21) * (xN31 * yN51 - xN51 * yN31)) - (xN31 * (yN54 * (x - xN4) - xN54 * (y - yN4))) / ((xN21 * yN51 - xN51 * yN21) * (xN31 * yN51 - xN51 * yN31));

            p_N_p_x(5) = (yN23 * (yN21 * (x - xN1) - xN21 * (y - yN1))) / ((xN21 * yN61 - xN61 * yN21) * (xN23 * yN63 - xN63 * yN23)) + (yN21 * (yN23 * (x - xN3) - xN23 * (y - yN3))) / ((xN21 * yN61 - xN61 * yN21) * (xN23 * yN63 - xN63 * yN23));
            p_N_p_y(5) = -(xN23 * (yN21 * (x - xN1) - xN21 * (y - yN1))) / ((xN21 * yN61 - xN61 * yN21) * (xN23 * yN63 - xN63 * yN23)) - (xN21 * (yN23 * (x - xN3) - xN23 * (y - yN3))) / ((xN21 * yN61 - xN61 * yN21) * (xN23 * yN63 - xN63 * yN23));
            for (size_t ig = 0; ig < 6; ++ig)
            {
                int y = 0;
                if (p_N_p_x(ig) > 1)
                {
                    cout << "element no: " << i << ", the dN_over_dx(" << ig << ") > 1;\n";
                    y++;
                }
                if (p_N_p_y(ig) > 1)
                {
                    cout << "element no: " << i << ", the dN_over_dy(" << ig << ") > 1;\n";
                    y++;
                }
                if (y != 0)
                    cout << endl;
            }*/

            /*
            double a_i, a_j, a_k, b_i, b_j, b_k, c_i, c_j, c_k, A, edg1, edg2, edg3, p, x, y, x_i, x_j, x_k, y_i, y_j, y_k;
            edg1 = (pnt0 - pnt2).norm();
            edg2 = (pnt2 - pnt4).norm();
            edg3 = (pnt4 - pnt0).norm();
            p = (edg1 + edg2 + edg3) / 2.;
            x = Center(0);
            y = Center(1);

            A = pow(p * (p - edg1) * (p - edg2) * (p - edg3), 0.5);
            x_i = pnt0(0);
            x_j = pnt2(0);
            x_k = pnt4(0);

            y_i = pnt0(1);
            y_j = pnt2(1);
            y_k = pnt4(1);

            a_i = x_j * y_k - x_k * y_j;
            a_j = x_k * y_i - x_i * y_k;
            a_k = x_i * y_j - x_j * y_i;

            b_i = y_j - y_k;
            b_j = y_k - y_i;
            b_k = y_i - y_j;

            c_i = x_k - x_j;
            c_j = x_i - x_k;
            c_k = x_j - x_i;

            p_N_p_x(0) = (b_i * (a_i + b_i * x + c_i * y)) / (2 * pow(A, 2)) + (b_i * ((a_i + b_i * x + c_i * y) / A - 1)) / (2 * A);
            p_N_p_y(0) = (c_i * (a_i + b_i * x + c_i * y)) / (2 * pow(A, 2)) + (c_i * ((a_i + b_i * x + c_i * y) / A - 1)) / (2 * A);

            p_N_p_x(1) = (b_j * (a_j + b_j * x + c_j * y)) / (2 * pow(A, 2)) + (b_j * ((a_j + b_j * x + c_j * y) / A - 1)) / (2 * A);
            p_N_p_y(1) = (c_j * (a_j + b_j * x + c_j * y)) / (2 * pow(A, 2)) + (c_j * ((a_j + b_j * x + c_j * y) / A - 1)) / (2 * A);

            p_N_p_x(2) = (b_k * (a_k + b_k * x + c_k * y)) / (2 * pow(A, 2)) + (b_k * ((a_k + b_k * x + c_k * y) / A - 1)) / (2 * A);
            p_N_p_y(2) = (c_k * (a_k + b_k * x + c_k * y)) / (2 * pow(A, 2)) + (c_k * ((a_k + b_k * x + c_k * y) / A - 1)) / (2 * A);

            p_N_p_x(3) = (b_j * (a_i + b_i * x + c_i * y)) / pow(A, 2) + (b_i * (a_j + b_j * x + c_j * y)) / pow(A, 2);
            p_N_p_y(3) = (c_j * (a_i + b_i * x + c_i * y)) / pow(A, 2) + (c_i * (a_j + b_j * x + c_j * y)) / pow(A, 2);

            p_N_p_x(4) = (b_k * (a_j + b_j * x + c_j * y)) / pow(A, 2) + (b_j * (a_k + b_k * x + c_k * y)) / pow(A, 2);
            p_N_p_y(4) = (c_k * (a_j + b_j * x + c_j * y)) / pow(A, 2) + (c_j * (a_k + b_k * x + c_k * y)) / pow(A, 2);

            p_N_p_x(5) = (b_k * (a_i + b_i * x + c_i * y)) / pow(A, 2) + (b_i * (a_k + b_k * x + c_k * y)) / pow(A, 2);
            p_N_p_y(5) = (c_k * (a_i + b_i * x + c_i * y)) / pow(A, 2) + (c_i * (a_k + b_k * x + c_k * y)) / pow(A, 2);*/

            double xi, eta;
            Vector2d A, B, C, D;
            A << 0, 0;
            B << 1, 0;
            C << 0, 1;
            D = (A + B + C) / 3;
            xi = D(0);
            eta = D(1);

            double N1 = (1 - xi - eta) * (1 - 2 * xi - 2 * eta);
            double N2 = xi * (2 * xi - 1);
            double N3 = eta * (2 * eta - 1);
            double N4 = 4 * xi * (1 - xi - eta);
            double N5 = 4 * xi * eta;
            double N6 = 4 * eta * (1 - xi - eta);

            Eigen::VectorXd x_e, y_e, N;
            x_e = Eigen::VectorXd::Zero(6);
            y_e = x_e;
            N = x_e;
            x_e << pnt0(0), pnt1(0), pnt2(0), pnt3(0), pnt4(0), pnt5(0);
            y_e << pnt0(1), pnt1(1), pnt2(1), pnt3(1), pnt4(1), pnt5(1);
            //1 - 4 - 2 -5 -3 - 6
            N << N1, N4, N2, N5, N3, N6;

            Center(0) = x_e.dot(N);
            Center(1) = y_e.dot(N);
            Center(2) = 0;

            Eigen::VectorXd pd_N_over_pd_xi, pd_N_over_pd_eta;
            pd_N_over_pd_xi = Eigen::VectorXd::Zero(6);
            pd_N_over_pd_eta = pd_N_over_pd_xi;

            double N_natural_a_1 = 4 * xi + 4 * eta - 3;
            double N_natural_a_2 = 4 * xi - 1;
            double N_natural_a_3 = 0;
            double N_natural_a_4 = 4 - 8 * xi - 4 * eta;
            double N_natural_a_5 = 4 * eta;
            double N_natural_a_6 = -4 * eta;

            //1 - 4 - 2 -5 -3 - 6
            pd_N_over_pd_xi << N_natural_a_1, N_natural_a_4, N_natural_a_2, N_natural_a_5, N_natural_a_3, N_natural_a_6;

            double N_natural_b_1 = 4 * xi + 4 * eta - 3;
            double N_natural_b_2 = 0;
            double N_natural_b_3 = 4 * eta - 1;
            double N_natural_b_4 = -4 * xi;
            double N_natural_b_5 = 4 * xi;
            double N_natural_b_6 = 4 - 8 * eta - 4 * xi;
            //1 - 4 - 2 -5 -3 - 6
            pd_N_over_pd_eta << N_natural_b_1, N_natural_b_4, N_natural_b_2, N_natural_b_5, N_natural_b_3, N_natural_b_6;

            Eigen::VectorXd pd_x_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_xi, pd_y_over_pd_eta;
            pd_x_over_pd_xi = Eigen::VectorXd::Zero(6);
            pd_x_over_pd_eta = pd_x_over_pd_xi;
            pd_y_over_pd_xi = pd_x_over_pd_xi;
            pd_y_over_pd_eta = pd_x_over_pd_xi;

            pd_x_over_pd_xi = pd_N_over_pd_xi.transpose() * x_e;
            pd_x_over_pd_eta = pd_N_over_pd_eta.transpose() * x_e;

            pd_y_over_pd_xi = pd_N_over_pd_xi.transpose() * y_e;
            pd_y_over_pd_eta = pd_N_over_pd_eta.transpose() * y_e;

            Eigen::MatrixXd Jacobi;
            Jacobi = Eigen::MatrixXd::Zero(2, 2);
            Jacobi << pd_x_over_pd_xi, pd_y_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_eta;

            Eigen::MatrixXd tem_o;
            tem_o = Eigen::MatrixXd::Zero(2, pd_N_over_pd_xi.rows());
            for (size_t op = 0; op < (size_t)tem_o.cols(); ++op)
            {
                tem_o(0, op) = pd_N_over_pd_xi(op);
                tem_o(1, op) = pd_N_over_pd_eta(op);
            }

            Eigen::MatrixXd AAA;
            AAA = Jacobi.inverse() * tem_o;

            for (size_t op = 0; op < 6; ++op)
            {
                p_N_p_x(op) = AAA(0, op);
                p_N_p_y(op) = AAA(1, op);
            }

            //------------------------

            qx = -K_coe * p_N_p_x * h_e; // qx = -delta h / delta x * k,  or say, qx = - ix * k
            qy = -K_coe * p_N_p_y * h_e;

            //------

            Vector3d q_overall;
            q_overall << qx, qy, 0;
            Flux_rate_2D[i][j].first = Center;
            Flux_rate_2D[i][j].second = q_overall;
            //cout << "Element NO: " << j << ", q: " << pow(pow(q_overall(0), 2) + pow(q_overall(1), 2), 0.5) << ", " << q_overall << "; ";

            //-----rotation to 3D
            double R_angle_temp1 = DFN_mesh.Rota_angle[i].first(0);
            Vector3d Frac_center;
            Frac_center << DFN_mesh.Rota_angle[i].first(1), DFN_mesh.Rota_angle[i].first(2), DFN_mesh.Rota_angle[i].first(3);
            Vector3d temp3;
            temp3 << DFN_mesh.Rota_angle[i].second(3), DFN_mesh.Rota_angle[i].second(4), DFN_mesh.Rota_angle[i].second(5);
            Vector3d Normal_frac;
            Normal_frac << DFN_mesh.Rota_angle[i].second(0), DFN_mesh.Rota_angle[i].second(1), DFN_mesh.Rota_angle[i].second(2);

            if (abs(Normal_frac(0)) < 0.0001 && abs(Normal_frac(1)) < 0.0001 && abs(Normal_frac(2)) < 0.0001)
            {
                //nothing to do
                Center += Frac_center;
            }
            else
            {
                Quaternion_t Q_axis_1;
                NormalizeRotation(R_angle_temp1, temp3, Q_axis_1);

                Vector3d temp4;
                Rotation(q_overall, Q_axis_1, temp4);
                q_overall = temp4;

                Vector3d temp5;
                Rotation(Center, Q_axis_1, temp5);
                Center = temp5 + Frac_center;
            }
            //cout << "after rotation, q: " << pow(pow(q_overall(0), 2) + pow(q_overall(1), 2) + pow(q_overall(2), 2), 0.5) << ", " << q_overall << endl;
            Flux_rate[i][j].first = Center;
            Flux_rate[i][j].second = q_overall;
        }
    }
};

inline void FEM_DFN::Verify_in_and_out_flux(DFN::DFN_mesh DFN_mesh, double *X_overall)
{
    double Q_in = 0, Q_out = 0;
    double q_in = 0, q_out = 0;
    double length_in = 0, length_out = 0;
    cout << "in: \n";
    for (size_t i = 0; i < DFN_mesh.Inlet.size(); ++i)
    {
        std::map<std::pair<size_t, size_t>, double>::iterator inlet_A = DFN_mesh.Inlet[i].begin();

        while (inlet_A != DFN_mesh.Inlet[i].end())
        {
            size_t i_frac = inlet_A->first.first;
            size_t j_ele = inlet_A->first.second;
            double edge_length = inlet_A->second;
            Vector3d q_vector = Flux_rate[i_frac][j_ele].second;
            double q = pow(pow(q_vector(0), 2) + pow(q_vector(1), 2) + pow(q_vector(2), 2), 0.5);
            double Q_ele = q * edge_length;
            length_in += edge_length;
            Q_in += Q_ele;
            q_in += q;
            cout << "element no: " << j_ele << ", Q: " << Q_ele << "; flux rate:" << q << "\n";
            inlet_A++;
        }
    }
    cout << "\n\nout: \n";
    for (size_t i = 0; i < DFN_mesh.Outlet.size(); ++i)
    {
        std::map<std::pair<size_t, size_t>, double>::iterator outlet_A = DFN_mesh.Outlet[i].begin();
        while (outlet_A != DFN_mesh.Outlet[i].end())
        {
            size_t i_frac = outlet_A->first.first;
            size_t j_ele = outlet_A->first.second;
            double edge_length = outlet_A->second;
            Vector3d q_vector = Flux_rate[i_frac][j_ele].second;
            double q = pow(pow(q_vector(0), 2) + pow(q_vector(1), 2) + pow(q_vector(2), 2), 0.5);
            double Q_ele = q * edge_length;
            length_out += edge_length;
            Q_out += Q_ele;
            q_out += q;
            cout << "element no: " << j_ele << ", Q: " << Q_ele << "; flux rate:" << q << "\n";
            outlet_A++;
        }
    }
    std::cout << "\n\nQ_in: " << Q_in << ", trace length: " << length_in << "; "
              << "Q_out: " << Q_out << ", trace length: " << length_out << "; " << std::endl;
    std::cout << "q_in: " << q_in << ", q_out: " << q_out << std::endl;
};

}; // namespace DFN