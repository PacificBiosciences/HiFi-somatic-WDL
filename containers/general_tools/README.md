* This uses docker buildx to build for multiple architecture or single architecture.
* The issue with building for arm64 is that a lot of conda packages are not available!
  * `docker buildx create --name mybuilder --use --bootstrap`
  * `docker buildx build --platform linux/amd64 -t kpinpb/severus:v0.1 .`
