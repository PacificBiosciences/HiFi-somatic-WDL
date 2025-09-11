VERSION=v0.1.1

docker buildx build . --platform linux/amd64 -t kpinpb/tmb-calculator:${VERSION} --progress=plain
docker tag kpinpb/tmb-calculator:${VERSION} kpinpb/tmb-calculator:latest
docker push kpinpb/tmb-calculator:latest && docker push kpinpb/tmb-calculator:${VERSION}
