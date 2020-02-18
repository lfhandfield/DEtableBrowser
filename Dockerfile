FROM rocker/shiny:3.5.1

# install some R required stuff
RUN apt-get update -y --no-install-recommends \
    && apt-get -y install -f \
       zlib1g-dev \
       libssl-dev \
       libcurl4-openssl-dev \
       wget \
       libxml2-dev \
       && apt-get clean && \
       rm -rf /var/lib/apt/lists/*

ADD app app/
COPY *.R app/
RUN install2.r data.table DT devtools ggplot2

#COPY detablebrowser.conf /etc/detablebrowser/detablebrowser.conf
CMD ["/usr/bin/detablebrowser.sh"]
