FROM rocker/shiny:3.5.1
ENV DOWNLOAD_URL http://midori-browser.org/downloads/midori_0.5.11-0_amd64_.deb

RUN \
  sed -i 's/http.debian.net/ftp.ch.debian.org/' /etc/apt/sources.list && \
  apt-get update && \
  apt-get install -y \
    ca-certificates \
    fonts-dejavu \
    gnome-icon-theme \
    libcanberra-gtk-module \
    wget && \
  wget "${DOWNLOAD_URL}" -O midori.deb && \
  { dpkg -i midori.deb || true; } && \
  apt-get update && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y -f && \
  apt-get purge -y --auto-remove wget && \
  rm -rf /var/lib/apt/lists/*
  
#RUN \
#  sed -i 's/http.debian.net/ftp.ch.debian.org/' /etc/apt/sources.list && \
#  apt-get update && \
#  apt-get install -y \
#    ca-certificates \
#    fonts-dejavu \
#    gnome-icon-theme \
#    libcanberra-gtk-module \
#    wget && \
#  wget "${DOWNLOAD_URL}" -O midori.deb && \
#  { dpkg -i midori.deb || true; } && \
#  apt-get update && \
#  DEBIAN_FRONTEND=noninteractive apt-get install -y -f && \
#  apt-get purge -y --auto-remove wget && \
#  rm -rf /var/lib/apt/lists/*
#RUN \
#  sed -i 's/http.debian.net/ftp.ch.debian.org/' /etc/apt/sources.list && \
#  apt-get update -y --no-install-recommends && \
#  apt-get install -y -f \
#    zlib1g-dev \
#    libssl-dev \
#    libcurl4-openssl-dev \
#    libxml2-dev \
#    ca-certificates \
#    fonts-dejavu \
#    gnome-icon-theme \
#    libcanberra-gtk-module \
#    wget && \
#  wget "${DOWNLOAD_URL}" -O midori.deb && \
#  { dpkg -i midori.deb || true; } && \
#  apt-get update && \
#  DEBIAN_FRONTEND=noninteractive apt-get install -y -f && \
#  apt-get purge -y --auto-remove wget && \
#  apt-get clean && \
#  rm -rf /var/lib/apt/lists/* && \
#  install2.r data.table DT devtools ggplot2

RUN \
  export uid=1000 gid=1000 && \
  groupadd --gid ${gid} user && \
  useradd --uid ${uid} --gid ${gid} --create-home user
  
COPY app/*.R /srv/shiny-server/
COPY server.conf /etc/shiny-server/shiny-server.conf
  
CMD ["/usr/bin/shiny-server.sh"]
