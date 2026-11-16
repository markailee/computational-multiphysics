#pragma once

#include <algorithm>
#include <cmath>
#include <istream>
#include <ostream>

namespace mp {

struct Vector2 {
    double x = 0.0;
    double y = 0.0;

    double operator[](int i) const { return i == 0 ? x : y; }
    double& operator[](int i) { return i == 0 ? x : y; }

    Vector2& operator+=(const Vector2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2& operator-=(const Vector2& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2& operator*=(double value) {
        x *= value;
        y *= value;
        return *this;
    }

    Vector2& operator/=(double value) {
        x /= value;
        y /= value;
        return *this;
    }
};

inline Vector2 operator+(Vector2 a, const Vector2& b) { return a += b; }
inline Vector2 operator-(Vector2 a, const Vector2& b) { return a -= b; }
inline Vector2 operator-(const Vector2& a) { return {-a.x, -a.y}; }
inline Vector2 operator*(Vector2 a, double value) { return a *= value; }
inline Vector2 operator*(double value, Vector2 a) { return a *= value; }
inline Vector2 operator/(Vector2 a, double value) { return a /= value; }
inline double dot(const Vector2& a, const Vector2& b) { return a.x * b.x + a.y * b.y; }
inline double cross(const Vector2& a, const Vector2& b) { return a.x * b.y - a.y * b.x; }
inline double normSquared(const Vector2& a) { return dot(a, a); }
inline double norm(const Vector2& a) { return std::sqrt(normSquared(a)); }
inline Vector2 normalized(const Vector2& a) { return a / norm(a); }

inline std::ostream& operator<<(std::ostream& stream, const Vector2& value) {
    return stream << value.x << ' ' << value.y;
}

inline std::istream& operator>>(std::istream& stream, Vector2& value) {
    return stream >> value.x >> value.y;
}

} // namespace mp
