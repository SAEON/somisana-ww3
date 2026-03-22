# This Dockerfile provides instructions to build a docker image
# which forms the image for running the command line interface - cli.py in the root of this repo
# the cli allows us to call whatever python functions in this repo that we
# need for our WW3 operations
#
# The docker image is built by
# .github/workflows/build_image.yml
# and gets triggered whenever any changes are made to this repo

# start with a miniforge base image
# We're fixing the miniforge version due to dependency issues which get introduced if left open ended
FROM condaforge/miniforge3:26.1.0-0

ENV DEBIAN_FRONTEND noninteractive

RUN mkdir /somisana-ww3
WORKDIR /somisana-ww3

# Install somisana-ww3 environment into base conda environment
COPY environment.yml .
RUN mamba env update -n base -f environment.yml

# add the somisana-ww3 code and install into the base environment
ADD . /somisana-ww3
RUN pip install -e .

# Set the cli.py as the entry point
ENTRYPOINT ["python", "cli.py"]
