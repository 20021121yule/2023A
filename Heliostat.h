#ifndef HELIOSTAT_H
#define HELIOSTAT_H

#include <vector>

#include "angles_sun.h"
#include "Eigen/Dense"

class Heliostat : public angles_sun {
public:
    double cord_x, cord_y, _theta; // 定日镜坐标和角度（范围 0-2*pi）
    int _circle; // 所在圈数
    double _length, _width, _height; // 固定参数：长、宽、高
    double square; // 面积
    double alpha_heliostat = 0, gamma_heliostat = 0; // 定日镜动态参数（弧度制）

    double cos_efficient = 0; // 余弦效率
    double shadow_ratio = 0; // 阴影效率
    double trunc = 0; // 截断效率

    double optical_efficiency = 0; // 该定日镜光学效率
    double sum_energy_per_heliostat = 0; // 该定日镜在所有时间内总的产生的能量

    Heliostat(const double x, const double y, const double theta, const int circle, const double l, const double w,
              const double h, const double st_calculated,
              const int days_calculated) : angles_sun(st_calculated, days_calculated), cord_x(x), cord_y(y),
                                           _theta(theta), _circle(circle), _length(l), _width(w), _height(h) {
        surface_coordinates = {
            {-l / 2, w / 2, 0}, {l / 2, w / 2, 0},
            {l / 2, -w / 2, 0}, {-l / 2, -w / 2, 0}
        }; // 初始化该定日镜镜面坐标系下4个端点的坐标
        square = l * w; // 计算该定日镜面积
    }

    // 找到最临近的定日镜
    // 在保存的所有的定日镜序列里面找
    // 返回指向最临近定日镜的指针
    [[nodiscard]] Heliostat *find_closest_heliostat(const std::vector<Heliostat *> &heliostats) const;

    // 对该定日镜做所有的计算
    // 需要外界给入参数：指向最临近的定日镜的指针
    // 直接调用即可计算
    void calculate_operations_on_heliostat(const Heliostat *closest_heliostat) {
        calculate_central_norm_and_field_coordinates();
        calculate_projection_and_reflection_coordinates_in_field(closest_heliostat);
        calculate_trunc_and_shadow_ratio();
        calculate_optical_efficiency();
        calculate_energy_per_heliostat_on_sum_and_average();
    }

private:
    Eigen::Vector3d central_norm; // 中心点法向量

    std::vector<Eigen::Vector3d> field_coordinates; // 该定日镜镜场坐标系下4个端点的坐标
    std::vector<Eigen::Vector3d> surface_coordinates; // 该定日镜镜面坐标系下4个端点的坐标
    std::vector<Eigen::Vector3d> projection_coordinates = {}; // 前一个定日镜顺着阳光在此定日镜上的投影点坐标（面坐标系）
    std::vector<Eigen::Vector3d> reflection_coordinates = {}; // 前一个定日镜逆着阳光在此定日镜上的反射点坐标（面坐标系）

    // 吸收塔以及集热器参数
    // 默认其在中心
    double Height_collector = 80.0; // 吸收塔高度
    double Height_tower = 8.0; // 集热器高度
    double diameter = 7.0; // 集热器直径

    // 遮挡效率
    // 计算定日镜场面积
    static constexpr double r_sq = 337.032 * 337.032 + (-8.21) * (-8.21); // 最大半径
    static constexpr double s_field = M_PI * r_sq; // 面积
    const double block = (80 + 8) * 7 / (tan(alpha_sun) * s_field); // block效率

    // 私有函数
    // 计算中心法向量（central_norm）以及四个端点在镜场坐标系下的坐标
    void calculate_central_norm_and_field_coordinates();

    // 计算在镜场坐标系下的投射（projection）和反射（reflection）坐标
    void calculate_projection_and_reflection_coordinates_in_field(const Heliostat *adjacent_heliostat);

    // 计算截断（trunc）和阴影（shadow）效率
    void calculate_trunc_and_shadow_ratio();

    // 计算光学效率
    void calculate_optical_efficiency();

    // 计算总能量和平均能量
    void calculate_energy_per_heliostat_on_sum_and_average();

    // 镜场/镜面坐标系相互转换的函数
    [[nodiscard]] Eigen::Vector3d field_trans_to_surf(const Eigen::Vector3d &cords) const; // 镜场坐标系转化成镜面坐标系

    // 生成定日镜面上在镜场坐标系下的离散坐标序列
    // 返回一系列坐标，这些坐标是一个序列。
    std::vector<Eigen::Vector3d> generate_meshgrid_cords_in_field();

    // 几何辅助函数
    // 根据围成封闭图形的4个点，按照顺序返回其4条边的方程(系数)
    // coordinates是一系列点的坐标，这里是4个
    // 返回这4个点围成的四边形的4条方程
    static std::vector<Eigen::Vector4d> calculate_line_equations(const std::vector<Eigen::Vector3d> &quad);
    // 判定一个点(coordinate)是否在围成的封闭图形当中
    // coordinate是点的坐标，polygon_line_equations是封闭图形边缘的直线方程
    // 返回bool值，0:不在四边形当中，1:在四边形当中
    static bool is_cord_in_polygon(const Eigen::Vector3d &cords, const std::vector<Eigen::Vector4d> &polygon);

    // 计算以vec方向为轴的旋转矩阵
    static Eigen::Matrix3d calculate_rotation_matrix_along_vec(const Eigen::Vector3d &vec);
};


#endif // HELIOSTAT_H
