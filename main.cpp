#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <format>
#include "Heliostat.h"

// 全局参数
double global_sum_cos_eff = 0.0, global_sum_s_b_eff = 0.0, global_sum_trunc_eff = 0.0, global_sum_optical_eff = 0.0; // 全局效率参数
int global_count_cos_eff = 0, global_count_s_b_eff = 0, global_count_trunc_eff = 0, global_count_optical_eff = 0; // 全局计数参数（计算平均值）
double sum_energy_total = 0.0; // 总能量参数
double average_energy_total_per_heliostat = 0.0; // 单个定日镜的总能量参数

void global_eff_parameters_update(const Heliostat *heliostat) {
    // 将 heliostat 计算后的结果来更新全局参数
    // 更新全局各种效率总和
    global_sum_cos_eff +=  heliostat->cos_efficient;
    global_sum_s_b_eff += 1 - heliostat->shadow_ratio;
    global_sum_trunc_eff += heliostat->trunc;
    global_sum_optical_eff += heliostat->optical_efficiency;
    sum_energy_total += heliostat->sum_energy_per_heliostat;

    // 计数
    global_count_cos_eff++;
    global_count_s_b_eff++;
    global_count_trunc_eff++;
    global_count_optical_eff++;
}

std::vector<Heliostat *> read_and_calculate_heliostats_at_day_and_st(const std::string &file_path,
                                                                     const int &days_calculated,
                                                                     const double &st_calculated) {
    std::vector<Heliostat *> heliostats; // 用指针储存
    // 加载数据文件
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Error: 无法打开文件!");
    }
    std::string line;
    std::getline(file, line); // 跳过标题行
    while (std::getline(file, line)) {
        // 检查并去除可能的 \r
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        std::istringstream iss(line);
        if (std::string x_str, y_str, theta_str, circle_str;
            getline(iss, x_str, ',') && getline(iss, y_str, ',') && getline(iss, theta_str, ',') && getline(
                iss, circle_str, ',')) {
            const double x = stod(x_str);
            const double y = stod(y_str);
            const double theta = stod(theta_str);
            const int circle = stoi(circle_str);
            heliostats.emplace_back(new Heliostat(x, y, theta, circle, 6.0, 6.0, 4.0, st_calculated,
                                                 days_calculated));
        }
    }
    file.close();

    for ( auto const &heliostat:heliostats) {
        const Heliostat *closest_heliostat = nullptr;
        if (heliostat->_circle > 1) { closest_heliostat = heliostat->find_closest_heliostat(heliostats); }
        heliostat->calculate_operations_on_heliostat(closest_heliostat);
        global_eff_parameters_update(heliostat);
    }

    return heliostats;
}

int main() {
    // 初始化太阳角
    const std::vector days_calculated = {-58, -27, 0, 30, 60, 91, 121, 152, 183, 213, 244, 274};
    const std::vector st_calculated = {9.0, 10.5, 12.0, 13.5, 15.0};

    const std::string file_path = "/Users/yule/Desktop/2023_A_完整代码/A题/附件_结果.csv"; // 文件路径

    for (auto const &d: days_calculated) {
        for (auto const &s: st_calculated) {
            // 生成在 d 和 s 下的一组定日镜
            const std::vector<Heliostat *> heliostats_curr =
                         read_and_calculate_heliostats_at_day_and_st(file_path, d, s);
        }
    }

    // 全局平均值
    const double global_average_cos_eff_ratio = global_sum_cos_eff / global_count_cos_eff;
    const double global_average_s_b_eff_ratio = global_sum_s_b_eff / global_count_s_b_eff;
    const double global_average_trunc_eff_ratio = global_sum_trunc_eff / global_count_trunc_eff;
    const double global_average_optical_eff_ratio = global_sum_optical_eff / global_count_optical_eff;

    std::cout << std::format("The cos efficient is: {:.4f}", global_average_cos_eff_ratio) << std::endl;
    std::cout << std::format("The s&b efficient is: {:.4f}", global_average_s_b_eff_ratio) << std::endl;
    std::cout << std::format("The trunc efficient is: {:.4f}", global_average_trunc_eff_ratio) << std::endl;
    std::cout << std::format("The optical efficient is: {:.4f}", global_average_optical_eff_ratio) << std::endl;
    std::cout << std::format("The total energy per month is: {:.4f}", sum_energy_total / (1e3 * 5.0 * 12.0)) << " MW" << std::endl;

    return 0;
}
