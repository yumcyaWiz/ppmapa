#ifndef _CAMERA_H
#define _CAMERA_H
#include <spdlog/spdlog.h>

#include <cmath>

#include "core.h"

class Camera {
 protected:
  Vec3f position;
  Vec3f forward;
  Vec3f right;
  Vec3f up;

 public:
  Camera(const Vec3f& position, const Vec3f& forward)
      : position(position), forward(forward) {
    right = normalize(cross(forward, Vec3f(0, 1, 0)));
    up = normalize(cross(right, forward));

    spdlog::info("[Camera] position: ({}, {}, {})", position[0], position[1],
                 position[2]);
    spdlog::info("[Camera] forward: ({}, {}, {})", forward[0], forward[1],
                 forward[2]);
    spdlog::info("[Camera] right: ({}, {}, {})", right[0], right[1], right[2]);
    spdlog::info("[Camera] up: ({}, {}, {})", up[0], up[1], up[2]);
  }

  // sample ray emitting from the given sensor coordinate
  // NOTE: uv: [-1, -1] x [1, 1], sensor coordinate
  virtual bool sampleRay(const Vec2f& uv, Ray& ray, float& pdf) const = 0;
};

// pinhole camera
class PinholeCamera : public Camera {
 private:
  float FOV;
  float focal_length;

 public:
  PinholeCamera(const Vec3f& position, const Vec3f& forward,
                float FOV = 0.5f * PI)
      : Camera(position, forward) {
    // compute focal length from FOV
    focal_length = 1.0f / std::tan(0.5f * FOV);
    spdlog::info("[PinholeCamera] focal_length: {}", focal_length);
  }

  bool sampleRay(const Vec2f& uv, Ray& ray, float& pdf) const override {
    const Vec3f pinholePos = position + focal_length * forward;
    const Vec3f sensorPos = position + uv[0] * right + uv[1] * up;
    ray = Ray(sensorPos, normalize(pinholePos - sensorPos));
    pdf = 1.0f;
    return true;
  }
};

#endif