#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include <omp.h>

#include <optional>

#include "core.h"
#include "photon_map.h"
#include "scene.h"

class Integrator {
 public:
  // compute radiance coming from the given ray
  virtual Vec3f integrate(const Ray& ray, const Scene& scene,
                          Sampler& sampler) const = 0;

  // compute cosine term
  // NOTE: need to account for the asymmetry of BSDF when photon tracing
  // https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation#x3-Non-symmetryDuetoShadingNormals
  // Veach, Eric. Robust Monte Carlo methods for light transport simulation.
  // Stanford University, 1998. Section 5.3
  static float cosTerm(const Vec3f& wo, const Vec3f& wi,
                       const SurfaceInfo& surfaceInfo,
                       const TransportDirection& transport_dir) {
    const float wi_ns = dot(wi, surfaceInfo.shadingNormal);
    const float wi_ng = dot(wi, surfaceInfo.geometricNormal);
    const float wo_ns = dot(wo, surfaceInfo.shadingNormal);
    const float wo_ng = dot(wo, surfaceInfo.geometricNormal);

    // prevent light leaks
    if (wi_ng * wi_ns <= 0 || wo_ng * wo_ns <= 0) {
      return 0;
    }

    if (transport_dir == TransportDirection::FROM_CAMERA) {
      return std::abs(wi_ns);
    } else if (transport_dir == TransportDirection::FROM_LIGHT) {
      return std::abs(wo_ns) * std::abs(wi_ng) / std::abs(wo_ng);
    } else {
      spdlog::error("invalid transport direction");
      std::exit(EXIT_FAILURE);
    }
  }
};

// implementation of path tracing
// NOTE: for reference purpose
class PathTracing : public Integrator {
 private:
  const int maxDepth;

 public:
  PathTracing(int maxDepth = 100) : maxDepth(maxDepth) {}

  void build(const Scene& scene, Sampler& sampler) override {}

  Vec3f integrate(const Ray& ray_in, const Scene& scene,
                  Sampler& sampler) const override {
    Vec3f radiance(0);
    Ray ray = ray_in;
    Vec3f throughput(1, 1, 1);

    for (int k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // russian roulette
        if (k > 0) {
          const float russian_roulette_prob = std::min(
              std::max(throughput[0], std::max(throughput[1], throughput[2])),
              1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob) {
            break;
          }
          throughput /= russian_roulette_prob;
        }

        // Le
        if (info.hitPrimitive->hasAreaLight()) {
          radiance += throughput *
                      info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        // sample direction by BxDF
        Vec3f dir;
        float pdf_dir;
        Vec3f f = info.hitPrimitive->sampleBxDF(
            -ray.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA,
            sampler, dir, pdf_dir);

        // update throughput and ray
        throughput *= f *
                      cosTerm(-ray.direction, dir, info.surfaceInfo,
                              TransportDirection::FROM_CAMERA) /
                      pdf_dir;
        ray = Ray(info.surfaceInfo.position, dir);
      } else {
        break;
      }
    }

    return radiance;
  }
};

class SPPM : public Integrator {
 private:
  // number of iterations
  const int nIteration;
  // number of photons in each iteration
  const int nPhotons;
  // initial global radius
  const float initialRadius;

  float globalRadius;
  PhotonMap photonMap;

  // compute reflected radiance with photon map
  Vec3f computeRadianceWithPhotonMap(const Vec3f& wo,
                                     const IntersectInfo& info) const {
    // get nearby photons
    float max_dist2;
    const std::vector<int> photon_indices =
        globalPhotonMap.queryKNearestPhotons(info.surfaceInfo.position,
                                             nEstimationGlobal, max_dist2);

    Vec3f Lo;
    for (const int photon_idx : photon_indices) {
      const Photon& photon = globalPhotonMap.getIthPhoton(photon_idx);
      const Vec3f f = info.hitPrimitive->evaluateBxDF(
          wo, photon.wi, info.surfaceInfo, TransportDirection::FROM_CAMERA);
      Lo += f * photon.throughput;
    }
    if (photon_indices.size() > 0) {
      Lo /= (nPhotonsGlobal * PI * max_dist2);
    }
    return Lo;
  }

  // sample initial ray from light and compute initial throughput
  Ray sampleRayFromLight(const Scene& scene, Sampler& sampler,
                         Vec3f& throughput) {
    // sample light
    float light_choose_pdf;
    const std::shared_ptr<Light> light =
        scene.sampleLight(sampler, light_choose_pdf);

    // sample point on light
    float light_pos_pdf;
    const SurfaceInfo light_surf = light->samplePoint(sampler, light_pos_pdf);

    // sample direction on light
    float light_dir_pdf;
    const Vec3f dir =
        light->sampleDirection(light_surf, sampler, light_dir_pdf);

    // spawn ray
    Ray ray(light_surf.position, dir);
    throughput = light->Le(light_surf, dir) /
                 (light_choose_pdf * light_pos_pdf * light_dir_pdf) *
                 std::abs(dot(dir, light_surf.shadingNormal));

    return ray;
  }

  void buildPhotonMap(uint64_t seed) {
    std::vector<Photon> photons;

    // init sampler for each thread
    std::vector<std::unique_ptr<Sampler>> samplers(omp_get_max_threads());
    for (int i = 0; i < samplers.size(); ++i) {
      samplers[i] = sampler.clone();
      samplers[i]->setSeed(seed * (i + 1));
    }

    // build photon map
    spdlog::info("tracing photons...");
#pragma omp parallel for
    for (int i = 0; i < nPhotonsGlobal; ++i) {
      auto& sampler_per_thread = *samplers[omp_get_thread_num()];

      // sample initial ray from light and set initial throughput
      Vec3f throughput;
      Ray ray = sampleRayFromLight(scene, sampler_per_thread, throughput);

      // trace photons
      // whener hitting diffuse surface, add photon to the photon array
      // recursively tracing photon with russian roulette
      for (int k = 0; k < maxDepth; ++k) {
        if (std::isnan(throughput[0]) || std::isnan(throughput[1]) ||
            std::isnan(throughput[2])) {
          spdlog::error("photon throughput is NaN");
          break;
        } else if (throughput[0] < 0 || throughput[1] < 0 ||
                   throughput[2] < 0) {
          spdlog::error("photon throughput is minus");
          break;
        }

        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
          if (bxdf_type == BxDFType::DIFFUSE) {
            // TODO: remove lock to get more speed
#pragma omp critical
            {
              photons.emplace_back(throughput, info.surfaceInfo.position,
                                   -ray.direction);
            }
          }

          // russian roulette
          if (k > 0) {
            const float russian_roulette_prob = std::min(
                std::max(throughput[0], std::max(throughput[1], throughput[2])),
                1.0f);
            if (sampler_per_thread.getNext1D() >= russian_roulette_prob) {
              break;
            }
            throughput /= russian_roulette_prob;
          }

          // sample direction by BxDF
          Vec3f dir;
          float pdf_dir;
          const Vec3f f = info.hitPrimitive->sampleBxDF(
              -ray.direction, info.surfaceInfo, TransportDirection::FROM_LIGHT,
              sampler_per_thread, dir, pdf_dir);

          // update throughput and ray
          throughput *= f *
                        cosTerm(-ray.direction, dir, info.surfaceInfo,
                                TransportDirection::FROM_LIGHT) /
                        pdf_dir;
          ray = Ray(info.surfaceInfo.position, dir);
        } else {
          // photon goes to the sky
          break;
        }
      }
    }

    // build photon map
    spdlog::info("building global photon map");
    globalPhotonMap.setPhotons(photons);
    globalPhotonMap.build();
  }
};

#endif