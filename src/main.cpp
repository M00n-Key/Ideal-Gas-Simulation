
#include <iostream>
#include <cmath>

#include <random>
#include <ctime>
#include <string>
#include <algorithm>

#include <fstream>
#include <filesystem>

#include "Vector.hpp"

const float deltatime = 1e-6;
// Половина размера ячейки, в которой происходит симуляция
// [-length, length]x[-length, length]x[-length, length]
const float length = 1.e-3;
const float bolzman_constant = 1.38 * 1.e-23;

Vector force(const Vector& source, const Vector& target) {
    static const float epsilon = 120 * bolzman_constant;
    static const float sigma = 0.34 * 1.e-9;

    // На частицу может действовать не только сама частица, но и ее "образ" при периодическом отображении
    // Поэтому надо учитывать всех ближайших соседей. Большая часть из них будет отсечена по расстоянию (см. ниже)
    std::vector<Vector> possible_source = {
        source, source + Vector(length, 0, 0), source + Vector(length, 0, 0), 
        source + Vector(0, length, 0), source + Vector(0, -length, 0), 
        source + Vector(0, 0, length), source + Vector(0, 0, -length), 
        source + Vector(length, length, 0), source + Vector(length, -length, 0), 
        source + Vector(-length, length, 0), source + Vector(-length, -length, 0), 
        source + Vector(length, 0, length), source + Vector(length, 0, -length),
        source + Vector(-length, 0, length), source + Vector(-length, 0, -length),
        source + Vector(0, length, length), source + Vector(0, length, -length),
        source + Vector(0, -length, length), source + Vector(0, -length, -length),
        source + Vector(length, length, length), source + Vector(length, length, -length),
        source + Vector(length, -length, length), source + Vector(length, -length, -length),
        source + Vector(-length, length, length), source + Vector(-length, length, -length),
        source + Vector(-length, -length, length), source + Vector(-length, -length, -length),
    };

    std::transform(possible_source.cbegin(), possible_source.cend(), possible_source.begin(), [target](const Vector vec){return vec - target;});

    // Обрубаем потенциал на значении r = 2.5 sigma
    std::remove_if(possible_source.begin(), possible_source.end(), [](const Vector& vec){return vec.magnitude() > 2.5 * sigma;});

    std::vector<float> distances(possible_source.size());
    std::transform(possible_source.cbegin(), possible_source.cend(), distances.begin(), [](const Vector& vec){return vec.magnitude();});

    for (size_t i = 0; i < possible_source.size(); i++) {
        possible_source[i] = (48 * epsilon / sigma / sigma) * (powf(sigma / distances[i], 14) - powf(sigma / distances[i], 8) / 2) * possible_source[i];
    }

    return std::accumulate(possible_source.cbegin(), possible_source.cend(), zero, [](const Vector& left, const Vector& right){return left + right;});
}

void clamp_to_box(Vector& vec) {
    // Ограничивает по x
    while (vec.x > length) {
        vec.x -= length;
    }
    while (vec.x < -length) {
        vec.x += length;
    }
    // Ограничивает по y
    while (vec.y > length) {
        vec.y -= length;
    }
    while (vec.y < -length) {
        vec.y += length;
    } 
    // Ограничивает по z
    while (vec.z > length) {
        vec.z -= length;
    }
    while (vec.z < -length) {
        vec.z += length;
    }
}

int main() {
    const float mass = 18 * 1.66 * 1.e-27;

    const int number_of_particles = 1000;

    Vector positions[number_of_particles];
    Vector velocities[number_of_particles];

    Vector next_positions[number_of_particles];
    Vector next_velocities[number_of_particles];
    
    Vector forces[number_of_particles];

    int seed = time(0);

    std::default_random_engine reng(seed);
    std::uniform_real_distribution<float> dist(-length, length); // распределение положения частиц
    float max_possible_velocity = sqrtf(2 * bolzman_constant * 300 / mass);
    std::uniform_real_distribution<float> vel_dist(-1, 1); // распределение скоростей частиц
    std::normal_distribution norm_dist(max_possible_velocity, max_possible_velocity / 2); // распределение модуля вектора

    std::cout << "Initialize points" << std::endl;

    // Случайная инициализация точек и скоростей
    for (int i = 0; i < number_of_particles; i++) {
        float x_init = dist(reng), y_init = dist(reng), z_init = dist(reng);
        float phi_init = vel_dist(reng), theta_init = vel_dist(reng);

        positions[i].x = x_init;
        positions[i].y = y_init;
        positions[i].z = z_init;

        velocities[i].x = cos(phi_init) * cos(theta_init);
        velocities[i].y = sin(phi_init) * cos(theta_init);
        velocities[i].z = sin(theta_init);

        velocities[i] = norm_dist(reng) * velocities[i];
    }

    std::cout << "Done" << std::endl;

    for (int step = 0; step < 10; step++) {

        std::cout << "Step " << step << ":" << std::endl;

        std::cout << "Calculate forces..." << std::endl;

        // считаем силы, действующие на частицы
        for (int i = 0; i < number_of_particles; i++) {
            forces[i] = zero;
            for (int j = 0; j < number_of_particles; j++) {
                if (i == j) continue;

                forces[i] += force(positions[i], positions[j]);
            }
        }

        std::cout << std::endl << "Done." << std::endl;

        // алгоритм Верле

        // считаем новые позиции

        std::cout << "Calculate positions..." << std::endl;
        for (int i = 0; i < number_of_particles; i++) {
            next_positions[i] = positions[i] + deltatime * velocities[i] + deltatime * deltatime / 2 / mass * forces[i];
            clamp_to_box(next_positions[i]);
        }
        std::cout << std::endl << "Done." << std::endl;

        std::cout << "Calculate new forces...";
        // снова считаем силы
        Vector new_forces[number_of_particles];
        for (int i = 0; i < number_of_particles; i++) {
            new_forces[i] = zero;
            for (int j = 0; j < number_of_particles; j++) {
                if (i == j) continue;

                new_forces[i] += force(next_positions[i], next_positions[j]);
            }
        }
        std::cout << std::endl << "Done." << std::endl;

        // считаем скорости
        std::cout << "Calculate velocities...";
        for (int i = 0; i < number_of_particles; i++) {
            velocities[i] = velocities[i] + deltatime / 2 * (forces[i] + new_forces[i]);
        }
        std::cout << "Done." << std::endl;

        std::cout << "Swap positions";
        for (int i = 0; i < number_of_particles; i++) {
            std::swap(next_positions[i], positions[i]);
        }
        std::cout << std::endl << "Done." << std::endl;

        // Записываем в файл
        std::cout << "Write a file..."; 
        std::filesystem::create_directory("output");
        std::ofstream file(std::string("output/step") + std::to_string(step) + ".csv");
        if (!file.is_open()) continue;
        for (int i = 0; i < number_of_particles; i++) {
            file << positions[i].x << " " << positions[i].y << " " << positions[i].z << " "
                << velocities[i].x << " " << velocities[i].y << " " << velocities[i].z << std::endl;
        }
        std::cout << std::endl << "Done." << std::endl;
    }

    return 0;
}

