VERSION=v0.1

docker buildx build . --platform linux/amd64 -t kpinpb/tmb-calculator:${VERSION} --progress=plain
# Then docker run IMAGEID
# Use docker ps -l to get commit ID, example b0e1c82cba75 below
docker commit b0e1c82cba75 quay.io/pacbio/tmb_calculator
docker commit b0e1c82cba75 quay.io/pacbio/tmb_calculator:${VERSION}
docker push quay.io/pacbio/tmb_calculator && docker push quay.io/pacbio/tmb_calculator:${VERSION}
