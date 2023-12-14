- The following command build to quay.io

```bash
docker buildx build --platform linux/amd64 -t pacbio/purple:v4.0 . --progress=plain --load &> build.log
# Use docker images to find image ID, then run docker run IMAGEID, and use docker ps -l to find container ID
# for commit
docker commit 1d82deb169aa quay.io/pacbio/purple:v4.0
docker commit 1d82deb169aa quay.io/pacbio/purple
docker push quay.io/pacbio/purple
docker push quay.io/pacbio/purple:v4.0
```
