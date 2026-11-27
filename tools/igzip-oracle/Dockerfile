# ISA-L 2.31.0 built from source for x86-64, the bytes crates/gkl-deflate/src/igzip.rs ports (decision 0050).
FROM --platform=linux/amd64 debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates git gcc make nasm autoconf automake libtool \
 && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 --branch v2.31.0 https://github.com/intel/isa-l /src/isa-l \
 && cd /src/isa-l && ./autogen.sh >/dev/null && ./configure --prefix=/usr >/dev/null && make -j8 >/dev/null && make install >/dev/null && ldconfig
COPY oracle.c /src/oracle.c
RUN gcc -O2 -o /usr/local/bin/igzip-oracle /src/oracle.c -lisal
