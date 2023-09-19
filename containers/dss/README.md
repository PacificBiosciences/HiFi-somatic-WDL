docker build -t kpinpb/dss:v0.1 .
docker tag kpinpb/dss:v0.1 kpinpb/dss:latest
docker push kpinpb/dss:v0.1
docker push kpinpb/dss:latest
