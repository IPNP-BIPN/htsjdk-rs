# Oracle image for htsjdk-rs golden generation.
#
# Every input is pinned: the base image by digest, every jar by sha256. Rebuilding this file
# must produce a functionally identical environment, because the meaning of every golden in
# this repository depends on it.
#
# Platform is linux/amd64 on purpose. See docs/decisions/0004-oracle-platform.md: Intel GKL
# ships x86-only natives, and they do load under Rosetta translation, so this image is a valid
# oracle even on Apple Silicon. Build and run it with --platform linux/amd64.

FROM eclipse-temurin@sha256:068a8f9ae4b74d9a20de3ecca771ba1a6437f1d7a8a8ad6deaf9dbdd2274397a

ARG HTSJDK_VERSION=4.2.0
ARG GKL_VERSION=0.8.11

RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/oracle/jars

# commons-io and commons-logging are NOT optional. GKL extracts its native library using
# org.apache.commons.io.FilenameUtils; without it the load fails with a WARNING and silently
# falls back to the JDK deflater. That failure mode is the reason OracleProbe asserts.
RUN set -eux; \
    fetch() { curl -fsSL -o "$2" "$1"; echo "$3  $2" | sha256sum -c -; }; \
    M=https://repo1.maven.org/maven2; \
    fetch $M/com/github/samtools/htsjdk/${HTSJDK_VERSION}/htsjdk-${HTSJDK_VERSION}.jar \
          htsjdk.jar 52c9eb1a568d8261767ebf888d6ebafa60911bd44e0c9242d413fbff1d1e2398; \
    fetch $M/com/intel/gkl/gkl/${GKL_VERSION}/gkl-${GKL_VERSION}.jar \
          gkl.jar 612976f6f8c1fffe85bf143447dceab5410bdaa4d828f98bcc899959f7b475ba; \
    fetch $M/org/broadinstitute/gatk-native-bindings/1.0.0/gatk-native-bindings-1.0.0.jar \
          gatk-native-bindings.jar c9f5af7134cd649e00234e9c5c8310eb3e3854feb0b5e5cbe43de94412305593; \
    fetch $M/commons-logging/commons-logging/1.2/commons-logging-1.2.jar \
          commons-logging.jar daddea1ea0be0f56978ab3006b8ac92834afeefbd9b7e4e6316fca57df0fa636; \
    fetch $M/commons-io/commons-io/2.16.1/commons-io-2.16.1.jar \
          commons-io.jar f41f7baacd716896447ace9758621f62c1c6b0a91d89acee488da26fc477c84f; \
    fetch $M/org/apache/commons/commons-jexl/2.1.1/commons-jexl-2.1.1.jar \
          commons-jexl.jar 03c9a9fae5da78ce52c0bf24467cc37355b7e23196dff4839e2c0ff018a01306; \
    fetch $M/org/apache/commons/commons-math3/3.5/commons-math3-3.5.jar \
          commons-math3.jar dcf9768726b82e8d7b8a928447b08bca7cfe2773d5120a63a41ec224207c1a29; \
    fetch $M/org/apache/commons/commons-compress/1.26.0/commons-compress-1.26.0.jar \
          commons-compress.jar 051aceb8bbcc62d0f5b2b8ac72c53767f9c59bfbd050151e65bef6f51c8ed9c9

# commons-math3 is pinned to 3.5, not to the 3.6.1 the jmath corpus job downloads for its
# FastMath columns. 3.5 is what GATK pins (`strictly '3.5'` in its build.gradle), so it is the
# version whose numbers reach a golden downstream. The two pins disagreeing is a real gap and is
# recorded in the manifest rather than papered over.
# commons-compress is NOT optional either. `FastaReferenceWriter` counts the bytes it has
# written with `org.apache.commons.compress.utils.CountingOutputStream`, and every offset in
# the .fai it writes comes from that counter, so without the jar the class does not load at
# all: a NoClassDefFoundError rather than a degraded output. htsjdk 4.2.0 pins 1.26.0.
ENV ORACLE_CP=/opt/oracle/jars/htsjdk.jar:/opt/oracle/jars/gkl.jar:/opt/oracle/jars/gatk-native-bindings.jar:/opt/oracle/jars/commons-logging.jar:/opt/oracle/jars/commons-io.jar:/opt/oracle/jars/commons-jexl.jar:/opt/oracle/jars/commons-math3.jar:/opt/oracle/jars/commons-compress.jar:/opt/oracle/classes

WORKDIR /opt/oracle
COPY OracleProbe.java .
RUN mkdir -p classes && javac -cp "$ORACLE_CP" -d classes OracleProbe.java

# Fail the build itself if the environment does not satisfy the contract, so a degraded image
# can never reach the point of producing a golden.
RUN java -cp "$ORACLE_CP" OracleProbe

WORKDIR /work
ENTRYPOINT ["/bin/sh", "-c"]
CMD ["java -cp \"$ORACLE_CP\" OracleProbe"]
