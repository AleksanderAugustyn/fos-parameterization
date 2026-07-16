#!/usr/bin/env bash
# Build the manylinux2014 x86_64 wheel. Runs INSIDE
# quay.io/pypa/manylinux2014_x86_64 (CentOS 7, devtoolset-10 GCC 10.2):
#
#   docker run --rm -v "$PWD":/io -w /io quay.io/pypa/manylinux2014_x86_64 ./ci/build-wheel.sh
#
# CI passes SETUPTOOLS_SCM_PRETEND_VERSION_FOR_<DIST> via -e; local runs fall
# back to git describe. Identical local and CI by design.
set -euxo pipefail

PYBIN=/opt/python/cp312-cp312/bin
gfortran --version | head -1

# The bind-mounted repo belongs to the host user; git (setuptools_scm)
# refuses to read it otherwise.
git config --global --add safe.directory "$PWD"

rm -rf dist_raw wheelhouse

"$PYBIN/python" -m pip install --upgrade --quiet build
# nehalem == the x86-64-v2 feature floor (SSE4.2+POPCNT, no AVX) in GCC-10
# vocabulary: devtoolset-10 predates the psABI level names (GCC 11). Portable
# down to the cluster's Haswell nodes; the libs cost 0.1-0.4 ms/call, so
# portability wins over host-tuned codegen.
"$PYBIN/python" -m build --wheel --outdir dist_raw \
    -Ccmake.define.GCC_OPTS_MARCH=nehalem

# --plat makes auditwheel fail loudly if any symbol exceeds glibc 2.17.
auditwheel repair dist_raw/*.whl --plat manylinux2014_x86_64 -w wheelhouse

# Wheel-floor guard: a regression to the host tag must fail the build.
wheel_file=(wheelhouse/*.whl)
case "$(basename "${wheel_file[0]}")" in
    *manylinux_2_17_x86_64*|*manylinux2014_x86_64*)
        echo "wheel floor OK: ${wheel_file[0]}" ;;
    *)
        echo "ERROR: wheel is not manylinux2014: ${wheel_file[0]}" >&2
        exit 1 ;;
esac
