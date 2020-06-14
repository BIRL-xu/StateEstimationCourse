//
// Created by baul on 6/13/20.
//

#include "myICP/ICP.h"
#include <ros/ros.h>
#include<pcl/point_cloud.h>
#include<pcl_conversions/pcl_conversions.h>
#include<sensor_msgs/PointCloud2.h>
#include<pcl/io/pcd_io.h>
#include <memory>

int main (int argc, char **argv)
{
    ros::init (argc, argv, "myICP");
    ros::NodeHandle nh;
    ros::Publisher pcl_pub1 = nh.advertise<sensor_msgs::PointCloud2> ("points1", 1);
    ros::Publisher pcl_pub2 = nh.advertise<sensor_msgs::PointCloud2> ("points2", 1);
    ros::Publisher pcl_pub3 = nh.advertise<sensor_msgs::PointCloud2> ("transformedPoints", 1);
    pcl::PointCloud<pcl::PointXYZ> cloudFirst;
    pcl::PointCloud<pcl::PointXYZ> cloudSecond;

    pcl::io::loadPCDFile ("/home/baul/Shenlan/StateEstimationCourse/7/PCDdata/first.pcd", cloudFirst);
    pcl::io::loadPCDFile ("/home/baul/Shenlan/StateEstimationCourse/7/PCDdata/second.pcd", cloudSecond);

    //Convert the cloud to ROS message
    sensor_msgs::PointCloud2 outputFirst;
    sensor_msgs::PointCloud2 outputSecond;
    pcl::toROSMsg(cloudFirst, outputFirst);
    outputFirst.header.frame_id = "lidar";
    pcl::toROSMsg(cloudSecond, outputSecond);
    outputSecond.header.frame_id = "lidar";

    myICP::ICP doICP;
    doICP.setPointsCloud(cloudFirst, cloudSecond);
//    doICP.solve();
    doICP.solve2();
    sensor_msgs::PointCloud2 transformedSecondPoints;
    pcl::toROSMsg(doICP.transformToRefFrame(), transformedSecondPoints);
    transformedSecondPoints.header.frame_id = "lidar";
    ros::Rate loop_rate(1);
    while (ros::ok())
    {
        pcl_pub1.publish(outputFirst);
        pcl_pub2.publish(outputSecond);
        pcl_pub3.publish(transformedSecondPoints);
        ros::spinOnce();
        loop_rate.sleep();
    }
    return 0;
}