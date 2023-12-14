- Following commands will build the docker image to quay.io

``` bash
docker buildx build --platform linux/amd64 -t pacbio/somatic_r_tools:v0.1 . --progress=plain --load &> build.log
# Use docker images to find image ID, then run docker run IMAGEID, and use docker ps -l to find container ID
# for commit
docker commit 7a3473ef4a45 quay.io/pacbio/somatic_r_tools
docker commit 7a3473ef4a45 quay.io/pacbio/somatic_r_tools:v0.1
docker push quay.io/pacbio/somatic_r_tools
docker push quay.io/pacbio/somatic_r_tools:v0.1
```
