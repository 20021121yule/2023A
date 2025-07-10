//
// Created by 余乐 on 24-12-11.
//

#ifndef ANGLES_SUN_H
#define ANGLES_SUN_H

#include <cmath>

class angles_sun {
public:
    double alpha_sun; // 太阳高度角(角度制)
    double gamma_sun; // 太阳方位角(角度制)

    // 构造函数
    angles_sun(const double ST_input, const int D_input) : ST(ST_input), D(D_input) {
        // 按照顺序计算成员
        delta_s = delta_sun();
        omega_s = omega_sun();
        alpha_sun = calculate_alpha_sun();
        gamma_sun = calculate_gamma_sun();
    };

private:
    double delta_s; // 太阳赤纬角
    double phi_s = 39.4 * M_PI / 180; // 当地纬度
    double omega_s; // 太阳时角
    double ST; // 当地时间
    int D; // 天数

    [[nodiscard]] inline double delta_sun() const;
    [[nodiscard]] inline double omega_sun() const;
    [[nodiscard]] inline double calculate_alpha_sun() const;
    [[nodiscard]] inline double calculate_gamma_sun() const;
};

inline double angles_sun::delta_sun() const {
    return asin(sin(2 * M_PI * D / 365) * sin(2 * M_PI * 23.45 / 360));// 太阳赤纬角（黄赤交角）
}

inline double angles_sun::omega_sun() const{
    return M_PI * (ST - 12) / 12;
}

inline double angles_sun::calculate_alpha_sun() const{
    return asin(cos(delta_s) * cos(phi_s) * cos(omega_s) + sin(delta_s) * sin(phi_s));
}

inline double angles_sun::calculate_gamma_sun() const {
    double cos_gamma = (sin(delta_s) - sin(alpha_sun) * sin(phi_s)) / (cos(alpha_sun) * cos(phi_s));
    if (cos_gamma < -1) cos_gamma = -1; // 防止超出范围
    if (cos_gamma > 1) cos_gamma = 1;
    if (omega_s > 0) return (2 * M_PI - acos(cos_gamma));
    return acos(cos_gamma);
}


#endif //ANGLES_SUN_H
