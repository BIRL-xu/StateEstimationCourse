//
// Created by baul on 6/13/20.
//

#ifndef SRC_ICP_H
#define SRC_ICP_H

#include<pcl/point_cloud.h>
#include<pcl_conversions/pcl_conversions.h>

namespace myICP
{
    class ICP{
    public:
        ICP(){

        }

        void setPointsCloud(pcl::PointCloud<pcl::PointXYZ> firstPoints,
                            pcl::PointCloud<pcl::PointXYZ> secondPoints)
        {
            m_refPoints = firstPoints;
            m_secPoints = secondPoints;
        }

        /*变换矩阵解法*/
        bool solve();

        /*先求迭代求解旋转矩阵，再直接求解平移量*/
        bool solve2();

        // 将第二帧点云转换到第一帧点云坐标系下以验证所求(R,t)阵是否准确。
        pcl::PointCloud<pcl::PointXYZ> transformToRefFrame();
    private:
        Eigen::Vector3d getMeanPoint();

        // 由三个数据构造叉乘矩阵
        Eigen::Matrix3d constructPhat(const double x, const double y, const double z)
        {
            Eigen::Matrix3d phat;
            phat << 0, -z, y,
                    z, 0, -x,
                    -y, x, 0;
            return phat;
        }

        /*由Eigen三维矢量够着叉乘矩阵*/
        Eigen::Matrix3d constructPhat(const Eigen::Vector3d& vec)
        {
            return constructPhat(vec.x(), vec.y(), vec.z());
        }

        /*构造六维矢量的矩阵--变换矩阵解法使用*/
        Eigen::Matrix4d getEhat(const Eigen::Vector3d& trans, const Eigen::Vector3d& rot)
        {
            Eigen::Matrix4d eHat;
            eHat.block<3, 3>(0, 0) = constructPhat(rot.x(), rot.y(), rot.z());
            eHat.block<3, 1>(0, 3) = trans;
            eHat.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
            eHat(3, 3) = 0;
            return eHat;
        }

        // 4x4齐次矩阵的伴随矩阵(6x6)
        Eigen::MatrixXd adjointMatrix(const Eigen::Matrix4d& T)
        {
            Eigen::MatrixXd adT(6, 6);
            Eigen::Matrix3d C = T.block<3, 3>(0, 0);
            Eigen::Vector3d t = T.block<3, 1>(0, 3);
            adT.block<3, 3>(0, 0) = C;
            adT.block<3, 3>(0, 3) = constructPhat(t.x(), t.y(), t.z()) * C;
            adT.block<3, 3>(3, 0) = Eigen::Matrix3d::Zero();
            adT.block<3, 3>(3, 3) = C;
            return adT;
        }

        /*《机器人学状态估计》式8.83，　导致程序太慢，不使用*/
        Eigen::Vector3d calculate_b(const Eigen::Vector3d& p,
                const Eigen::Vector3d& y,
                const Eigen::Matrix3d& C,
                const std::uint32_t iterNum)
        {
            Eigen::Vector3d b = Eigen::Vector3d::Zero();
            for(std::uint32_t i = 0; i < iterNum; ++i)
            {
                Eigen::Vector3d pBias(m_refPoints.points[i].x - p.x(), m_refPoints.points[i].y - p.y(), m_refPoints.points[i].z - p.z());
                Eigen::Vector3d yBias(m_secPoints.points[i].x - y.x(), m_secPoints.points[i].y - y.y(), m_secPoints.points[i].z - y.z());
                b -= constructPhat(yBias) * C * pBias;
            }

            b /= iterNum;
            return b;
        }

    private:
        pcl::PointCloud<pcl::PointXYZ> m_refPoints;
        pcl::PointCloud<pcl::PointXYZ> m_secPoints;

        // 最终R, t阵
        Eigen::Matrix3d m_R = Eigen::Matrix3d::Identity();
        Eigen::Vector3d m_t = Eigen::Vector3d::Zero();
    };
}

#endif //SRC_ICP_H
