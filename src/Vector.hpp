#pragma once

#include <fstream>
#include <cmath>
#include <utility>

struct Vector {
    Vector() : Vector(0, 0, 0) {}
    Vector(float x, float y, float z) : x(x), y(y), z(z) {}
    Vector(const Vector& vec) : x(vec.x), y(vec.y), z(vec.z) {}
    Vector(Vector&& vec) {
        std::swap(x, vec.x);
        std::swap(y, vec.y);
        std::swap(z, vec.z);
    }

    Vector& operator=(const Vector& vec) {
        x = vec.x;
        y = vec.y;
        z = vec.z;
        return *this;
    }

    Vector& operator+=(const Vector& vec) {
        this->x += vec.x;
        this->y += vec.y;
        this->z += vec.z;
        return *this;
    }

    float magnitude() const {return sqrtf(x*x + y*y + z*z);}

    float x, y, z;
};

const Vector zero(0, 0, 0);

Vector operator+(const Vector& left, const Vector& right) {
    return Vector{left.x + right.x, left.y + right.y, left.z + right.z};
}

Vector operator*(float factor, const Vector& vec) {
    return Vector{factor * vec.x, factor * vec.y, factor * vec.z};
}

Vector operator-(const Vector& left, const Vector& right) {
    return Vector{left.x - right.x, left.y - right.y, left.z - right.z};
}

std::ostream& operator<<(std::ostream& out, const Vector& vec) {
    return out << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
}