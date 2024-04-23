FROM debian:buster-slim
RUN mkdir -p /usr/share/man/man1 /usr/share/man/man2 && \
apt-get update && \
apt-get install -y --no-install-recommends make build-essential libssl-dev  wget curl llvm libidn11  ca-certificates-java openjdk-11-jdk git nano tcsh sudo gfortran
WORKDIR "/opt"

RUN git clone --depth=1 https://github.com/pyenv/pyenv.git .pyenv
ENV PYENV_ROOT="/opt/.pyenv" 
ENV PATH="${PYENV_ROOT}/shims:${PYENV_ROOT}/bin:${PATH}"
RUN echo 'export PYENV_ROOT=/opt/.pyenv' >> ~/.bashrc 
RUN echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc 
RUN echo 'eval "$(pyenv init -)"' >> ~/.bashrc

# ENV PYTHON_VERSION2=miniconda3-4.7.12
# RUN pyenv install ${PYTHON_VERSION2} && \
#       pyenv global ${PYTHON_VERSION2}

RUN echo "Cloning diSBPred from github..."      
RUN git clone https://github.com/wasicse/diSBPred.git && \
chmod -R 777 /opt/diSBPred
WORKDIR "/opt/diSBPred"

RUN ./install_dependencies.sh
RUN ./download_largefiles.sh
RUN ./download_dataset.sh
ENTRYPOINT [ "/bin/bash" ]


