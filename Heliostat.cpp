#include "Heliostat.h"
#include <vector>
#include <cmath>
#include <stdexcept>
using namespace std;

// 构造函数已在头文件中内联定义

// 找到最接近的定日镜
Heliostat* Heliostat::find_closest_heliostat(const std::vector<Heliostat*>& heliostats) const {
    double min_diff = INFINITY;
    Heliostat* closest = nullptr;
    for (auto& heliostat : heliostats) {
        if (heliostat->_circle == this->_circle - 1) {
            if (const double diff = fabs(heliostat->_theta - this->_theta); diff < min_diff) {
                min_diff = diff;
                closest = heliostat;
            }
        }
    }
    return closest;
}

// 计算中心法向量和镜场坐标系下四个顶点
void Heliostat::calculate_central_norm_and_field_coordinates() {
    const double dist = std::sqrt(cord_x * cord_x + cord_y * cord_y + (Height_collector - _height) * (Height_collector - _height));
    const Eigen::Vector3d e_in(-std::cos(alpha_sun) * std::sin(gamma_sun), -std::cos(alpha_sun) * std::cos(gamma_sun), -std::sin(alpha_sun));
    const Eigen::Vector3d e_out(-cord_x / dist, -cord_y / dist, (Height_collector - _height) / dist);
    central_norm = (e_out - e_in).normalized();
    alpha_heliostat = std::acos(central_norm[2]);
    gamma_heliostat = std::atan2(central_norm[0], central_norm[1]);

    const double norm_xy = std::sqrt(central_norm[0] * central_norm[0] + central_norm[1] * central_norm[1]);
    const double norm_z = std::sqrt(1 - central_norm[2] * central_norm[2]);

    for (int ii = 0; ii < 4; ++ii) {
        const double sign_x = (ii == 1 || ii == 2) ? 1.0 : -1.0;
        const double sign_y = (ii == 0 || ii == 3) ? 1.0 : -1.0;
        const double sign_x1 = (ii == 0 || ii == 1) ? 1.0 : -1.0;
        const double x = cord_x + sign_x1 * (0.5 * _width * central_norm[2] * central_norm[0] / norm_xy)
                         + sign_x * (0.5 * _length * central_norm[1] / norm_xy);
        const double y = cord_y + sign_x1 * (0.5 * _width * central_norm[2] * central_norm[1] / norm_xy)
                         + sign_y * (0.5 * _length * central_norm[0] / norm_xy);
        const double z = _height + (ii < 2 ? -1 : 1) * (0.5 * _width * norm_z);
        field_coordinates.emplace_back(x, y, z);
    }
}

// 计算投影和反射点
void Heliostat::calculate_projection_and_reflection_coordinates_in_field(const Heliostat* adjacent_heliostat) {
    if (adjacent_heliostat == nullptr) return;
    const Eigen::Vector3d e_in(-std::cos(alpha_sun) * std::sin(gamma_sun), -std::cos(alpha_sun) * std::cos(gamma_sun), -std::sin(alpha_sun));
    Eigen::Vector3d e_out = (e_in - 2 * (e_in.dot(central_norm)) * central_norm).normalized();

    for (int ii = 0; ii < 4; ++ii) {
        const Eigen::Vector3d field_point = adjacent_heliostat->field_coordinates[ii];
        const Eigen::Vector3d point_to_heliostat = (field_point - Eigen::Vector3d(cord_x, cord_y, _height)).eval();

        const double denominator1 = central_norm.dot(e_in);
        const double denominator2 = central_norm.dot(e_out);
        if (fabs(denominator1) < 1e-6 || fabs(denominator2) < 1e-6) {
            throw std::runtime_error("Error: Denominator too close to zero!");
        }

        double t1 = -point_to_heliostat.dot(central_norm) / denominator1;
        double t2 = -point_to_heliostat.dot(central_norm) / denominator2;

        const Eigen::Vector3d projection_point = field_point + t1 * e_in;
        const Eigen::Vector3d reflection_point = field_point + t2 * e_out;

        projection_coordinates.push_back(field_trans_to_surf(projection_point));
        reflection_coordinates.push_back(field_trans_to_surf(reflection_point));
    }
}

// 计算截断效率和遮挡比例
void Heliostat::calculate_trunc_and_shadow_ratio() {
    const std::vector<Eigen::Vector3d> grid_points = generate_meshgrid_cords_in_field();
    const auto lineEquations1 = calculate_line_equations(projection_coordinates);
    const auto lineEquations2 = calculate_line_equations(reflection_coordinates);
    const Eigen::Vector3d e_in(-std::cos(alpha_sun) * std::sin(gamma_sun), -std::cos(alpha_sun) * std::cos(gamma_sun), -std::sin(alpha_sun));
    const Eigen::Matrix3d RotationMatrix = calculate_rotation_matrix_along_vec(e_in);

    int trunc_points = 0;
    int shadowed_points = 0;
    constexpr double alpha_max = 4.65 * 1e-3;
    constexpr double delta_1 = 2 * 1e-3;
    constexpr int n1 = static_cast<int>(alpha_max / delta_1);
    constexpr double gamma_max = 2 * M_PI;
    constexpr double delta_2 = 2.0;
    constexpr int n2 = static_cast<int>(gamma_max / delta_2);

    for (int ii = 0; ii <= n1; ++ii) {
        for (int jj = 0; jj <= n2; ++jj) {
            const double alpha_current = ii * delta_1;
            const double gamma_current = jj * delta_2;
            const Eigen::Vector3d e_in_shift(-sin(alpha_current) * cos(gamma_current),
                                             -sin(alpha_current) * sin(gamma_current), -cos(alpha_current));
            Eigen::Vector3d e_in_curr = (RotationMatrix * e_in_shift).normalized();
            Eigen::Vector3d e_out = (e_in_curr - 2 * (e_in_curr.dot(central_norm)) * central_norm).normalized();

            for (const auto& grid_vec : grid_points) {
                const Eigen::Vector3d grid_point_2d = field_trans_to_surf(grid_vec);

                if (_circle > 1) {
                    if (!is_cord_in_polygon(grid_point_2d, lineEquations1) &&
                        !is_cord_in_polygon(grid_point_2d, lineEquations2)) {
                        const double A = e_out.head<2>().squaredNorm();
                        const double B = 2 * e_out.head<2>().dot(grid_vec.head<2>());
                        const double C = grid_vec.head<2>().squaredNorm() - diameter * diameter / 4;
                        if (const double delta = B * B - 4 * A * C; delta > 0) {
                            const double t2 = (-B - std::sqrt(delta)) / (2 * A);
                            const double z_reflected = grid_vec[2] + t2 * e_out[2];
                            if (z_reflected >= Height_collector - 0.5 * Height_tower &&
                                z_reflected <= Height_collector + 0.5 * Height_tower) {
                                trunc_points++;
                            }
                        }
                    } else {
                        shadowed_points++;
                    }
                } else {
                    const double A = e_out.head<2>().squaredNorm();
                    const double B = 2 * e_out.head<2>().dot(grid_vec.head<2>());
                    const double C = grid_vec.head<2>().squaredNorm() - diameter * diameter / 4;
                    if (const double delta = B * B - 4 * A * C; delta > 0) {
                        const double t2 = (-B - std::sqrt(delta)) / (2 * A);
                        const double z_reflected = grid_vec[2] + t2 * e_out[2];
                        if (z_reflected >= Height_collector - 0.5 * Height_tower &&
                            z_reflected <= Height_collector + 0.5 * Height_tower) {
                            trunc_points++;
                        }
                    }
                }
            }
        }
    }

    const double total = static_cast<double>(grid_points.size()) * (n1 + 1) * (n2 + 1);
    trunc = static_cast<double>(trunc_points) / total;
    shadow_ratio = static_cast<double>(shadowed_points) / total + block;
}

// 计算光学效率
void Heliostat::calculate_optical_efficiency() {
    constexpr double ref_efficient = 0.92;
    constexpr double atm_efficient = 0.99;
    cos_efficient = std::cos(alpha_sun) * std::sin(gamma_sun) * central_norm[0] +
                    std::cos(alpha_sun) * std::cos(gamma_sun) * central_norm[1] +
                    std::sin(alpha_sun) * central_norm[2];
    optical_efficiency = ref_efficient * atm_efficient * cos_efficient * (1 - shadow_ratio) * trunc;
}

// 计算能量总量
void Heliostat::calculate_energy_per_heliostat_on_sum_and_average() {
    static constexpr double G_0 = 1.366;
    static constexpr double H = 3.0;
    static constexpr double a = 0.4237 - 0.00821 * (6.0 - H) * (6.0 - H);
    static constexpr double b = 0.5055 + 0.00595 * (6.5 - H) * (6.5 - H);
    static constexpr double c = 0.2711 + 0.01858 * (2.5 - H) * (2.5 - H);
    const double DNI = G_0 * (a + b * std::exp(-c / std::sin(alpha_sun)));
    sum_energy_per_heliostat = DNI * square * optical_efficiency;
}

// 镜场坐标系 → 镜面坐标系
Eigen::Vector3d Heliostat::field_trans_to_surf(const Eigen::Vector3d& cords) const {
    const double norm_xy = std::sqrt(central_norm[0] * central_norm[0] + central_norm[1] * central_norm[1]);
    const double norm_z = std::sqrt(1 - central_norm[2] * central_norm[2]);
    const Eigen::Vector3d U(central_norm[0] / norm_xy, central_norm[1] / norm_xy, 0.0);
    Eigen::Vector3d W(-_width * central_norm[2] * central_norm[0] / norm_xy + _length * central_norm[1] / norm_xy,
                      -_width * central_norm[2] * central_norm[1] / norm_xy - _length * central_norm[0] / norm_xy,
                      _width * norm_z);
    W.normalize();
    const Eigen::Vector3d P_diff = (cords - Eigen::Vector3d(cord_x, cord_y, _height)).eval();
    return {P_diff.dot(U), P_diff.dot(W), 0.0};
}

// 网格点生成
std::vector<Eigen::Vector3d> Heliostat::generate_meshgrid_cords_in_field() {
    constexpr double delta = 0.1;
    const int n1 = static_cast<int>(_length / delta);
    const int n2 = static_cast<int>(std::abs((_width * central_norm[2])) / delta);
    std::vector<Eigen::Vector3d> grid_points;
    Eigen::Matrix3d R_z;
    R_z << std::cos(gamma_heliostat), std::sin(gamma_heliostat), 0,
          -std::sin(gamma_heliostat), std::cos(gamma_heliostat), 0,
           0, 0, 1;

    for (int ii = 0; ii <= n1; ++ii) {
        for (int jj = 0; jj <= n2; ++jj) {
            const double x_local = -_length / 2 + ii * delta;
            const double y_local = -std::abs((_width * central_norm[2])) / 2 + jj * delta;
            Eigen::Vector3d local_cords_xy(x_local, y_local, 0);
            const Eigen::Vector3d global_cords_xy = (R_z * local_cords_xy + Eigen::Vector3d(cord_x, cord_y, 0.0)).eval();
            const double z = _height - (central_norm[0] * (global_cords_xy.x() - cord_x)
                                        + central_norm[1] * (global_cords_xy.y() - cord_y)) / central_norm[2];
            grid_points.emplace_back(global_cords_xy.x(), global_cords_xy.y(), z);
        }
    }
    return grid_points;
}

// 辅助函数：生成四边形的直线方程
std::vector<Eigen::Vector4d> Heliostat::calculate_line_equations(const std::vector<Eigen::Vector3d>& quad) {
    if (quad.empty()) return {};
    std::vector<Eigen::Vector4d> lineEquations;
    for (size_t i = 0; i < quad.size(); ++i) {
        const auto& P1 = quad[i];
        const auto& P2 = quad[(i + 1) % quad.size()];
        const Eigen::Vector3d norm_curr = (P2 - P1).cross(Eigen::Vector3d(0, 0, 1));
        const double D = -norm_curr.dot(P1);
        lineEquations.emplace_back(norm_curr.x(), norm_curr.y(), norm_curr.z(), D);
    }
    return lineEquations;
}

// 判断一个点是否在封闭区域内
bool Heliostat::is_cord_in_polygon(const Eigen::Vector3d& cords, const std::vector<Eigen::Vector4d>& polygon) {
    for (const auto& line : polygon) {
        if (const double result = line.head<3>().dot(cords) + line[3]; result > 0) {
            return false;
        }
    }
    return true;
}

// 计算旋转矩阵（Rodriguez）
Eigen::Matrix3d Heliostat::calculate_rotation_matrix_along_vec(const Eigen::Vector3d& vec) {
    Eigen::Vector3d v = vec.normalized();
    double x = v[0], y = v[1], z = v[2];
    Eigen::Vector3d u(-y, x, 0.0);
    if (u.norm() > 1e-6) u.normalize();
    else return Eigen::Matrix3d::Identity();

    double cos_theta = z;
    double sin_theta = std::sqrt(1 - z * z);
    Eigen::Matrix3d u_cross;
    u_cross << 0.0, -u[2], u[1],
               u[2], 0.0, -u[0],
              -u[1], u[0], 0.0;
    return (Eigen::Matrix3d::Identity() + sin_theta * u_cross + (1 - cos_theta) * u_cross * u_cross).eval();
}
