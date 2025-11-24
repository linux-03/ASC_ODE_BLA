# Use a Linux base — this is run inside Docker, even on macOS hosts
FROM ubuntu:24.04 AS pyodide_wheel

ENV DEBIAN_FRONTEND=noninteractive \
    PYODIDE_VERSION=0.26.1 \
    EMSDK_VERSION=3.1.58 \
    PYTHON_VERSION=3.12.1 \
    EMSDK_DIR=/opt/emsdk \
    VENV_DIR=/opt/pyodide-venv \
    PROJECT_DIR=/project

# Install system packages (one RUN to reduce layers)
RUN apt-get update \
 && apt-get -y upgrade \
 && apt-get install -y --no-install-recommends \
      nodejs ccache python3-pip python3 git make pkg-config g++ \
      lbzip2 xz-utils autoconf libtool unzip xxd wget python3-venv \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Install emscripten in /opt/emsdk
WORKDIR /opt
RUN git clone https://github.com/emscripten-core/emsdk.git ${EMSDK_DIR} \
 && cd ${EMSDK_DIR} \
 && ./emsdk install ${EMSDK_VERSION} \
 && ./emsdk activate ${EMSDK_VERSION}

# Create isolated venv for pyodide build tools
RUN python3 -m venv ${VENV_DIR} \
 && ${VENV_DIR}/bin/python -m pip install --upgrade pip \
 && ${VENV_DIR}/bin/pip install pyodide-build==${PYODIDE_VERSION} --break-system-packages

# Copy project source into image
WORKDIR ${PROJECT_DIR}
COPY . ${PROJECT_DIR}

SHELL ["/bin/bash", "-lc"]

# Source emsdk and venv then run pyodide build
RUN source ${EMSDK_DIR}/emsdk_env.sh \
 && source ${VENV_DIR}/bin/activate \
 && pyodide build

# Install the generated wheel (use wildcard so you don't hardcode the filename)
WORKDIR ${PROJECT_DIR}/dist
RUN source ${VENV_DIR}/bin/activate \
 && ${VENV_DIR}/bin/pip install --no-deps --force-reinstall dist/*.whl || true \
 && python3 -c "import sys; print('PYTHON', sys.version)"

# Stage: jupyterlite build (reuse built artifacts)
FROM pyodide_wheel AS jupyterlite_build

ENV VENV_DIR=/opt/pyodide-venv
WORKDIR /project

# install requirements (use venv pip)
RUN ${VENV_DIR}/bin/pip install -r requirements.txt --break-system-packages

# clone jupyter-lite site repo (same as you had)
WORKDIR /root
RUN git clone https://github.com/triadtitans/rigid_body_interactive.git
WORKDIR /root/rigid_body_interactive

# Use venv's jupyter lite to build the site; adjust pyodide tarball if needed
RUN TZ=UTC ${VENV_DIR}/bin/jupyter lite build --pyodide https://ngsolve.org/files/pyodide-0.26.0/ngsolve_pyodide.tar.bz2

# copy pyodide artifacts from the first stage into the site dir
RUN mkdir -p /root/rigid_body_interactive/dist/static/pyodide \
 && cp -r /project/dist/* /root/rigid_body_interactive/dist/static/pyodide/ || true

# Final runtime image with only the built site
FROM ubuntu:24.04 AS server
RUN apt-get update && apt-get install -y --no-install-recommends python3 && rm -rf /var/lib/apt/lists/*
WORKDIR /srv/site
COPY --from=jupyterlite_build /root/rigid_body_interactive/dist /srv/site
EXPOSE 8000
CMD ["python3", "-m", "http.server", "8000", "--directory", "/srv/site"]
