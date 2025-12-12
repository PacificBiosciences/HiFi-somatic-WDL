- This will build longphase to quay.io

```bash
docker buildx build --platform linux/amd64 -t pacbio/longphase:v1.5.2 . --progress=plain --load &> build.log
# Use docker images to find image ID, then run docker run IMAGEID, and use docker ps -l to find container ID
# for commit
docker commit 4c4a6caff388 quay.io/pacbio/longphase:v1.5.2
docker commit 4c4a6caff388 quay.io/pacbio/longphase
docker push quay.io/pacbio/longphase
docker push quay.io/pacbio/longphase:v1.5.2
```
