
#include <iostream>
#include <cmath>

#include <random>
#include <ctime>
#include <string>
#include <algorithm>

#include <fstream>
#include <filesystem>

#include "Vector.hpp"

// Половина размера ячейки, в которой происходит симуляция
// [-length, length]x[-length, length]x[-length, length]
const float length = 1.e-2;
const float bolzman_constant = 1.38 * 1.e-23;

// Параметры для аргона
const float epsilon = 120 * bolzman_constant;
const float sigma = 0.34 * 1.e-9;

const float mass = 18 * 1.66 * 1.e-27;

const float deltatime = 1e-5;

Vector force(const Vector& source, const Vector& target) {
    
    Vector sum_force = zero;

    // На частицу может действовать не только сама частица, но и ее "образ" при периодическом отображении
    // Поэтому надо учитывать всех ближайших соседей. Большая часть из них будет отсечена по расстоянию 
    for (int i = 0; i < 27; i++) {
        int digits[3] = {i / 9, i / 3 % 3, i % 3};
        for (int j = 0; j < 3; j++) {
            if (digits[j] == 2) digits[j] = -1;
        }
        
        Vector possible_source = source + length * Vector(digits[0], digits[1], digits[2]);

        float distance = (possible_source - target).magnitude();

        if (distance > 2.5 * sigma) continue;

        sum_force += (48 * epsilon / sigma / sigma) * (powf(sigma / distance, 14) - powf(sigma / distance, 8) / 2) * (target - possible_source);
    }

    return sum_force;
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
    std::uniform_real_distribution<float> vel_dist(-max_possible_velocity, max_possible_velocity); // распределение скоростей частиц

    std::cout << "Initialize points" << std::endl;

    // Случайная инициализация точек и скоростей
    for (int i = 0; i < number_of_particles; i++) {
        float x_init = dist(reng), y_init = dist(reng), z_init = dist(reng);
        float vx_init = vel_dist(reng), vy_init = vel_dist(reng), vz_init = vel_dist(reng);

        positions[i].x = x_init;
        positions[i].y = y_init;
        positions[i].z = z_init;

        velocities[i].x = vx_init;
        velocities[i].y = vy_init;
        velocities[i].z = vz_init;
    }

    std::cout << "Done" << std::endl;
    
    // Записываем в файл
    std::cout << "Write a file..."; 
    std::filesystem::create_directory("output");
    std::ofstream file(std::string("output/step") + std::to_string(0) + ".csv");
    if (!file.is_open()) return -1;
    for (int i = 0; i < number_of_particles; i++) {
        file << positions[i].x << " " << positions[i].y << " " << positions[i].z << " "
            << velocities[i].x << " " << velocities[i].y << " " << velocities[i].z << std::endl;
    }
    std::cout << "Done." << std::endl;


    for (int step = 1; step <= 10; step++) {

        std::cout << std::endl << "Step " << step << ":" << std::endl;

        std::cout << "Calculate forces...";

        // считаем силы, действующие на частицы
        for (int i = 0; i < number_of_particles; i++) {
            forces[i] = zero;
            for (int j = 0; j < number_of_particles; j++) {
                if (i == j) continue;

                forces[i] += force(positions[i], positions[j]);
            }
        }

        std::cout << "Done." << std::endl;

        // алгоритм Верле

        // считаем новые позиции

        std::cout << "Calculate positions...";
        for (int i = 0; i < number_of_particles; i++) {
            next_positions[i] = positions[i] + deltatime * velocities[i] + deltatime * deltatime / 2 / mass * forces[i];
            clamp_to_box(next_positions[i]);
        }
        std::cout << "Done." << std::endl;

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
        std::cout << "Done." << std::endl;

        // считаем скорости
        std::cout << "Calculate velocities...";
        for (int i = 0; i < number_of_particles; i++) {
            next_velocities[i] = velocities[i] + deltatime / 2 * (forces[i] + new_forces[i]);
        }
        std::cout << "Done." << std::endl;

        std::cout << "Swap positions and velocities...";
        for (int i = 0; i < number_of_particles; i++) {
            std::swap(next_positions[i], positions[i]);
            std::swap(next_velocities[i], velocities[i]);
        }
        std::cout << "Done." << std::endl;

        // Записываем в файл
        std::cout << "Write a file..."; 
        std::filesystem::create_directory("output");
        std::ofstream file(std::string("output/step") + std::to_string(step) + ".csv");
        if (!file.is_open()) continue;
        for (int i = 0; i < number_of_particles; i++) {
            file << positions[i].x << " " << positions[i].y << " " << positions[i].z << " "
                << velocities[i].x << " " << velocities[i].y << " " << velocities[i].z << std::endl;
        }
        std::cout << "Done." << std::endl;
    }

    return 0;
}

