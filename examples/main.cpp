#include "camera.h"
#include "integrator.h"
#include "scene.h"

int main() {
  const int width = 512;
  const int height = 512;
  const int n_iterations = 1;
  const int n_photons = 10000;
  const int max_depth = 100;

  Image image(width, height);

  const auto camera =
      std::make_shared<Camera>(Vec3f(0, 1, 6), Vec3f(0, 0, -1), 0.25 * PI);

  // build scene
  Scene scene;
  scene.loadModel("CornellBox-Original.obj");
  scene.build();

  // render
  UniformSampler sampler;
  PathTracing integrator(camera, 10000);
  integrator.render(scene, sampler, image);

  // gamma correction
  image.gammaCorrection(2.2f);

  // output image
  image.writePPM("output.ppm");

  return 0;
}