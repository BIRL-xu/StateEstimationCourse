//
// Created by baul on 6/13/20.
//

#include "myICP/ICP.h"
#include <Eigen/Eigen>
#include <cmath>

namespace myICP
{
    bool ICP::solve() {
        // Step1: 计算两个点云的重心
        if(m_refPoints.points.empty() || m_secPoints.points.empty())
        {
            ROS_ERROR("Empty points cloud.");
            return false;
        }
        // 取数量少的点个数作为匹配的数量
        const std::uint32_t pointsNum =
                (m_refPoints.points.size() <= m_secPoints.points.size())?\
                m_refPoints.points.size():m_secPoints.points.size();
        ROS_INFO("Matching points number: %d" ,pointsNum);

        Eigen::Vector3d p(0.0, 0.0, 0.0);
        Eigen::Vector3d y(0.0, 0.0, 0.0);
        for(std::uint32_t i = 0; i < pointsNum; ++i)
        {
            p.x() += m_refPoints.points[i].x / pointsNum;
            p.y() += m_refPoints.points[i].y / pointsNum;
            p.z() += m_refPoints.points[i].z / pointsNum;

            y.x() += m_secPoints.points[i].x / pointsNum;
            y.y() += m_secPoints.points[i].y / pointsNum;
            y.z() += m_secPoints.points[i].z / pointsNum;
        }
        std::cout << "p: \n" << p << std::endl << std::endl;
        std::cout << "y: \n" << y << std::endl << std::endl;
        // Step2: Construct M and I
        // P^
        Eigen::Matrix3d phat = constructPhat(p.x(), p.y(), p.z());
        std::cout << "phat: \n" << phat << std::endl << std::endl;

        // I
        Eigen::Matrix3d I;
        I << 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0;

        for(std::uint32_t i = 0; i < pointsNum; ++i)
        {
            auto temp = constructPhat(
                    m_refPoints.points[i].x - p.x(),
                    m_refPoints.points[i].y - p.y(),
                    m_refPoints.points[i].z - p.z());
            I -= temp * temp;
        }
        I /= pointsNum;
        std::cout << "I: \n" << I << std::endl << std::endl;
        // M
        Eigen::MatrixXd A(6, 6);
        A.setIdentity();
        A.block<3, 3>(0, 3) = phat;
/*        A.setIdentity();
        A.block<3, 3>(3, 0) = phat;*/

        Eigen::MatrixXd B(6, 6);
        B.setIdentity();
        B.block<3, 3>(3, 3) = I;
        Eigen::MatrixXd M(6, 6);
        M = A.transpose() * B * A;
//        M = A * B * A.transpose();
        std::cout << "M: \n" << M << std::endl << std::endl;

        // Step3: 迭代计算扰动e
        // 构造W
        Eigen::Matrix3d W;
        W << 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0,
             0.0, 0.0, 0.0;
        for(std::uint32_t i = 0; i < pointsNum; ++i)
        {
            Eigen::Vector3d yBias(
                    m_secPoints.points[i].x - y.x(),
                    m_secPoints.points[i].y - y.y(),
                    m_secPoints.points[i].z - y.z());
            Eigen::Vector3d pBias(
                    m_refPoints.points[i].x - p.x(),
                    m_refPoints.points[i].y - p.y(),
                    m_refPoints.points[i].z - p.z());
            Eigen::Matrix3d temp = yBias*pBias.transpose();
            W += temp;
        }
        W /= pointsNum;
        std::cout << "W: \n" << W << std::endl << std::endl;
        // 初始估计
        Eigen::Vector3d t(y.x() - p.x(),
                y.y() - p.y(),
                y.z() - p.z());
/*        Eigen::Vector3d t(y.x(),
                          y.y(),
                          y.z());*/
        std::cout << "initial-t: \n" << t << std::endl << std::endl;
//        Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + constructPhat(0, 0, 1);
        double initialAngle = 5 * M_PI / 180.0;
        Eigen::Matrix3d R;
        R << cos(initialAngle), -sin(initialAngle), 0,
             sin(initialAngle), cos(initialAngle), 0,
             0, 0, 1;
        // 齐次矩阵
        Eigen::Matrix4d T;
        T.block<3, 3>(0, 0) = R;
        T.block<3, 1>(0, 3) = -R * t;
        T.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
        T(3, 3) = 1;
        std::cout << "T: \n" << T << std::endl << std::endl;
        Eigen::Matrix3d b1 = constructPhat(1, 0, 0) * R * W.transpose();
        Eigen::Matrix3d b2 = constructPhat(0, 1, 0) * R * W.transpose();
        Eigen::Matrix3d b3 = constructPhat(0, 0, 1) * R * W.transpose();
        Eigen::Vector3d b(b1.trace(), b2.trace(), b3.trace());
        std::cout << "b: \n" << b << std::endl << std::endl;
        Eigen::VectorXd a(6, 1);
        a.block<3, 1>(0, 0) = y - R * (p - t);
        a.block<3, 1>(3, 0) = b - constructPhat(y.x(),y.y(),y.z()) * R * (p - t);
        std::cout << "a: \n" << a << std::endl << std::endl;
        // 初始扰动量
        Eigen::VectorXd e(6, 1);
        std::cout << "M.inverse(): \n" << M.inverse() << std::endl << std::endl;
        std::cout << "ad(T): \n" << adjointMatrix(T) << std::endl << std::endl;

        e = adjointMatrix(T) * M.inverse() * adjointMatrix(T).transpose() * a;

        std::cout << "initial-e: \n" << e << std::endl << std::endl;
        std::cout << "e-norm: \n" << e.norm() << std::endl << std::endl;

        //　迭代次数
        std::int32_t iterLimit = 30;
        std::int32_t iterNum = 0;
        // 最小扰动，即||e||_2 < minE时停止迭代
        const double minE = 1E-3;
        while(iterLimit && (e.norm() > minE))
        {
            // 计算扰动的李代数
            Eigen::Vector3d translation = e.block<3, 1>(0, 0);
            std::cout << "translation: \n" << translation << std::endl << std::endl;
            Eigen::Vector3d rotation = e.block<3, 1>(3, 0);
            std::cout << "rotation: \n" << rotation << std::endl << std::endl;

            // ehat
            Eigen::Matrix4d ehat = getEhat(translation, rotation);

            std::cout << "ehat: \n" << ehat << std::endl << std::endl;

            Eigen::Matrix4d I4;
            I4.setIdentity();
            std::cout << "I4: \n" << I4 << std::endl << std::endl;
            T = (I4 + ehat) * T;
            // 扰动e的旋转矩阵
/*            Eigen::Matrix3d C = Eigen::Matrix3d::Identity() + constructPhat(rotation.x(), rotation.y(), rotation.z());
            std::cout << "C: \n" << C << std::endl << std::endl;
            // 左雅克比矩阵
            Eigen::Matrix3d J = Eigen::Matrix3d::Identity() + 0.5*constructPhat(rotation.x(), rotation.y(), rotation.z());
            std::cout << "J: \n" << J << std::endl << std::endl;
            // 扰动e对应的李代数SE(3)
            Eigen::Matrix4d Te;
            Te.block<3, 3>(0, 0) = C;
            Te.block<3, 1>(0, 3) = J * translation;
            Te.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
            Te(3, 3) = 1;
            std::cout << "Te: \n" << Te << std::endl << std::endl;
            // 更新T
            T = Te * T;*/

            std::cout << "T: \n" << T << std::endl << std::endl;
            // 更新a
            R = (T.block<3, 3>(0, 0) + T.block<3, 3>(0, 0).transpose()) / 2.0;
            std::cout << "R: \n" << R << std::endl << std::endl;
            t = -R.inverse() * T.block<3, 1>(0, 3);
            std::cout << "t: \n" << t << std::endl << std::endl;

            b1 = constructPhat(1, 0, 0) * R * W.transpose();
            b2 = constructPhat(0, 1, 0) * R * W.transpose();
            b3 = constructPhat(0, 0, 1) * R * W.transpose();
            b << b1.trace(), b2.trace(), b3.trace();
            std::cout << "b: \n" << b << std::endl << std::endl;

            a.block<3, 1>(0, 0) = y - R * (p - t);
            a.block<3, 1>(3, 0) = b - constructPhat(y.x(),y.y(),y.z()) * R * (p - t);
            std::cout << "a: \n" << a << std::endl << std::endl;

            //更新扰动e
            e = adjointMatrix(T) * M.inverse() * adjointMatrix(T).transpose() * a;

            --iterLimit;
            ++iterNum;
            std::cout << "\tAfter iteration-T: \n" << T << std::endl << std::endl;
        }

        ROS_INFO("Iteration finished:");
        ROS_INFO("\tIteration times: %d", iterNum);
        std::cout << "\tAfter iteration-e: \n" << e << std::endl << std::endl;
        ROS_INFO("\tNormal of e: %.5f", e.norm());
        m_R = T.block<3, 3>(0, 0);
        m_t = T.block<3, 1>(0, 3);
        return true;
    }

    bool ICP::solve2() {
        // Step1: 计算两个点云的重心
        if(m_refPoints.points.empty() || m_secPoints.points.empty())
        {
            ROS_ERROR("Empty points cloud.");
            return false;
        }
        // 取数量少的点个数作为匹配的数量
        const std::uint32_t pointsNum =
                (m_refPoints.points.size() <= m_secPoints.points.size())?\
                m_refPoints.points.size():m_secPoints.points.size();
        ROS_INFO("Matching points number: %d" ,pointsNum);

        Eigen::Vector3d p(0.0, 0.0, 0.0);
        Eigen::Vector3d y(0.0, 0.0, 0.0);
        Eigen::Vector3d sumP = Eigen::Vector3d::Zero();
        Eigen::Vector3d sumy = Eigen::Vector3d::Zero();
        for(std::uint32_t i = 0; i < pointsNum; ++i)
        {
            p.x() += m_refPoints.points[i].x / pointsNum;
            p.y() += m_refPoints.points[i].y / pointsNum;
            p.z() += m_refPoints.points[i].z / pointsNum;

            y.x() += m_secPoints.points[i].x / pointsNum;
            y.y() += m_secPoints.points[i].y / pointsNum;
            y.z() += m_secPoints.points[i].z / pointsNum;
        }

        std::cout << "p: \n" << p << std::endl << std::endl;
        std::cout << "y: \n" << y << std::endl << std::endl;
        // Step2: Construct I, W, b
        // I
        Eigen::Matrix3d I;
        I << 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0;
        Eigen::Matrix3d W;
        W << 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0;

        for(std::uint32_t i = 0; i < pointsNum; ++i)
        {
            Eigen::Matrix3d temp = constructPhat(
                    m_refPoints.points[i].x - p.x(),
                    m_refPoints.points[i].y - p.y(),
                    m_refPoints.points[i].z - p.z());
            I -= temp * temp;

            Eigen::Vector3d yBias(
                    m_secPoints.points[i].x - y.x(),
                    m_secPoints.points[i].y - y.y(),
                    m_secPoints.points[i].z - y.z());
            Eigen::Vector3d pBias(
                    m_refPoints.points[i].x - p.x(),
                    m_refPoints.points[i].y - p.y(),
                    m_refPoints.points[i].z - p.z());
            Eigen::Matrix3d bias = yBias * pBias.transpose();
            W += bias;

        }
        I /= pointsNum;
        std::cout << "I: \n" << I << std::endl << std::endl;
        W /= pointsNum;
        std::cout << "W: \n" << W << std::endl << std::endl;
        // b
        // 初始估计
        Eigen::Matrix3d R;
/*        Eigen::Vector3d rotation(0.01, 0.01, 0.01);
        R = Eigen::Matrix3d::Identity() + constructPhat(rotation);*/

        const double initialAngle = 5 * M_PI / 180.0;
        R << cos(initialAngle), -sin(initialAngle), 0,
                sin(initialAngle), cos(initialAngle), 0,
                0, 0, 1;

        Eigen::Matrix3d b1 = constructPhat(1, 0, 0) * R * W.transpose();
        Eigen::Matrix3d b2 = constructPhat(0, 1, 0) * R * W.transpose();
        Eigen::Matrix3d b3 = constructPhat(0, 0, 1) * R * W.transpose();
        Eigen::Vector3d b(b1.trace(), b2.trace(), b3.trace());

        // 扰动phi
        Eigen::Vector3d  phi;
        phi = R * I.inverse() * R.transpose() * b;
        std::cout << "initial-phi: \n" << phi << std::endl << std::endl;
        //　迭代次数
        std::int32_t iterLimit = 50;
        std::int32_t iterNum = 0;
        // 最小扰动，即||e||_2 < minE时停止迭代
        const double minE = 1E-6;
        while(iterLimit && (phi.norm() > minE))
        {
            R = (Eigen::Matrix3d::Identity() + constructPhat(phi)) * R;
            Eigen::Matrix3d b1 = constructPhat(1, 0, 0) * R * W.transpose();
            Eigen::Matrix3d b2 = constructPhat(0, 1, 0) * R * W.transpose();
            Eigen::Matrix3d b3 = constructPhat(0, 0, 1) * R * W.transpose();
            b << b1.trace(), b2.trace(), b3.trace();

            phi = R * I.inverse() * R.transpose() * b;
            --iterLimit;
            iterNum++;
        }
        std::cout << "Iteration finished:" << std::endl;
        std::cout << "\tIteration times:" << iterNum << std::endl;
        std::cout << "\tphi: \n" << phi << std::endl << std::endl;
        std::cout << "\tphi-norm: \n" << phi.norm() << std::endl << std::endl;
        std::cout << "\tR: \n" << R << std::endl << std::endl;
        // 更新结果
        m_R = R;
        m_t = p - R.transpose() * y;
        std::cout << "\tt: \n" << m_t << std::endl << std::endl;

        return true;
    }

    pcl::PointCloud<pcl::PointXYZ> ICP::transformToRefFrame() {
        pcl::PointCloud<pcl::PointXYZ> out;
        out.header = m_secPoints.header;
        out.height = m_secPoints.height;
        out.width = m_secPoints.width;
        out.is_dense = m_secPoints.is_dense;

        // 齐次转换矩阵
        Eigen::Matrix4d T;
        T.setIdentity();
        T.block<3, 3>(0, 0) = m_R;
        T.block<3, 1>(0, 3) = m_t;
        std::cout << "T: \n" << T << std::endl << std::endl;
        for(auto &rawPoint : m_secPoints.points)
        {
            Eigen::Vector4d eigenPoint(rawPoint.x, rawPoint.y, rawPoint.z, 1);
            eigenPoint = T.inverse() * eigenPoint;

            pcl::PointXYZ pclPoint;
            pclPoint.x = eigenPoint.head<3>().x();
            pclPoint.y = eigenPoint.head<3>().y();
            pclPoint.z = eigenPoint.head<3>().z();
            out.points.push_back(pclPoint);
        }

        return out;
    }

} // namespace myICP