# Use wishartlab/phastest-docker-single as base
FROM docker.io/wishartlab/phastest-docker-single

# Root as workdir
WORKDIR /

# Install basics
RUN apt-get install unzip

# Remove old script subfolder, download needed scripts
RUN rm -r phastest-app && \
    wget https://phastest.ca/download_file/phastest-docker.zip && \
    unzip -qq phastest-docker.zip && \
    mv phastest/phastest-app-docker phastest-app && \
    rm phastest-docker.zip && \
    rm -r phastest/ && \
    rm -r __MACOSX && \
    cp -r phastest-app/sub_programs/ncbi-blast-2.3.0+/ BLAST+/

# Make files executable
RUN chmod -R 755 "/phastest-app" && \
    chmod -R 755 "/BLAST+"

# Copy modified files for input and database as parameters
COPY *.pl /phastest-app/scripts/
COPY phastest /usr/bin/phastest

# Entrypoint
CMD ["bash"]
