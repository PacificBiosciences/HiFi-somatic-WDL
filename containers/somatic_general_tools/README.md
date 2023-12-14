- The following command build to quay.io

```bash
docker buildx build --platform linux/amd64 -t pacbio/somatic_general-tools:v0.1 . --progress=plain --load &> build.log
# Use docker images to find image ID, then run docker run IMAGEID, and use docker ps -l to find container ID
# for commit
docker commit 54667beb3c81 quay.io/pacbio/somatic_general_tools:v0.1
docker commit 54667beb3c81 quay.io/pacbio/somatic_general_tools
docker push quay.io/pacbio/somatic_general_tools
docker push quay.io/pacbio/somatic_general_tools:v0.1
