#! /bin/bash

docker build -t wasicse/disbpred - < Dockerfile

# docker commit CONTAINERNAME  wasicse/disbpred:latest

docker push  wasicse/disbpred:latest