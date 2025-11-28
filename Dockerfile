# Stage 1: build pyodide wheel
FROM ubuntu:24.04 AS pyodide_wheel

ENV DEBIAN_FRONTEND=noninteractive \
    PYODIDE_VERSION=0.26.1 \
    EMSDK_VERSION=3.1.58 \
    PYTHON_VERSION=3.12.1 \
    EMSDK_DIR=/root/emsdk \
    VENV_DIR=/root/env26 \
    VENV_DIR_SETUP=/root/venv_setup \
    PROJECT_DIR=/root/ASC_ODE_BLA

# install packages in one step, keep layers small
RUN apt-get update \
 && apt-get -y upgrade \
 && apt-get install -y --no-install-recommends \
        curl \
        nodejs \
        ccache \
        python3-pip \
        python3 \
        git \
        make \
        pkg-config \
        g++ \
        lbzip2 \
        xz-utils \
        autoconf \
        libtool \
        unzip \
        xxd \
        wget \
        python3.12-venv \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /root
# clone emsdk and install specified version
RUN git clone https://github.com/emscripten-core/emsdk.git ${EMSDK_DIR} \
 && cd ${EMSDK_DIR} \
 && ./emsdk install ${EMSDK_VERSION} \
 && ./emsdk activate ${EMSDK_VERSION}

# Create isolated venv for pyodide build tools
RUN python3 -m venv ${VENV_DIR_SETUP} \
 && ${VENV_DIR_SETUP}/bin/python -m pip install --upgrade pip \
 && ${VENV_DIR_SETUP}/bin/pip install pyodide-build==${PYODIDE_VERSION} --break-system-packages


# Create a pyodide-build venv (pyodide venv command provided by pyodide-build)
# this creates a usable python environment for building/packaging
SHELL ["/bin/bash", "-lc"]

RUN source ${VENV_DIR_SETUP}/bin/activate \
 && which pyodide \
 && pyodide --version \
 && pyodide venv ${VENV_DIR}


# copy project into image (do this after heavy installs to keep cache)
WORKDIR ${PROJECT_DIR}
COPY . ${PROJECT_DIR}

SHELL ["/bin/bash", "-lc"]

# source emsdk and venv, then run pyodide build
RUN source ${EMSDK_DIR}/emsdk_env.sh \
 && source ${VENV_DIR}/bin/activate \
 && pyodide build

# go to dist and install wheel (use wildcard to avoid hardcoded filename)
WORKDIR ${PROJECT_DIR}/dist
RUN source ${VENV_DIR}/bin/activate \
 && sh -c 'ls -1 *.whl 2>/dev/null || true' \
 && sh -c 'if ls *.whl >/dev/null 2>&1; then ${VENV_DIR}/bin/pip install --no-deps --force-reinstall *.whl || true; else echo "No wheel found"; fi' \
 && python3 -c "import sys; print('PYTHON', sys.version)"

# Stage 2: jupyter-lite build (re-use artifacts from pyodide_wheel)
FROM pyodide_wheel AS jupyterlite_build

ENV VENV_DIR=/root/env26
WORKDIR ${PROJECT_DIR}

# install jupyter-lite requirements (use venv pip)
RUN ${VENV_DIR}/bin/pip install --upgrade pip \
 && if [ -f requirements.txt ]; then ${VENV_DIR}/bin/pip install -r requirements.txt --break-system-packages; fi

# clone the jupyter-lite site repo (adjust to your site repo)
WORKDIR /root
RUN git clone https://github.com/linux-03/rigid_body_inter.git rigid_body_inter

WORKDIR /root/rigid_body_inter
# Ensure jupyter-lite is available in venv, then build site with pyodide tarball
RUN ${VENV_DIR}/bin/pip install jupyterlite-core --break-system-packages \
 && TZ=UTC ${VENV_DIR}/bin/jupyter lite build --pyodide https://ngsolve.org/files/pyodide-0.26.0/ngsolve_pyodide.tar.bz2

# copy pyodide artifacts from first stage into the site static folder if present
RUN mkdir -p /root/rigid_body_inter/dist/static/pyodide \
 && cp -r /root/ASC_ODE_BLA/dist/* /root/rigid_body_inter/dist/static/pyodide/ || true

# Stage 3: runtime - serve the built site
FROM ubuntu:24.04 AS server
RUN apt-get update && apt-get install -y --no-install-recommends python3 && rm -rf /var/lib/apt/lists/*

WORKDIR /srv/site
# copy the built jupyter-lite site from previous stage
COPY --from=jupyterlite_build /root/rigid_body_inter/dist /srv/site

EXPOSE 8000
# simple, reproducible server
CMD ["python3", "-m", "http.server", "8000", "--directory", "/srv/site"]
