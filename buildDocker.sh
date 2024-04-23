#! /bin/bash

docker build -t wasicse/disbpred - < Dockerfile

# docker commit CONTAINERNAME  wasicse/esmdispred:latest

# docker push  wasicse/disbpred:latest