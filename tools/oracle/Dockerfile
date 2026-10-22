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

ARG HTSJDK_VERSION=5.0.0
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
          htsjdk.jar 59d862b74ec3e847ec84bc87bed34cecd1a80095495d965b9403f863e374e99c; \
    fetch $M/com/intel/gkl/gkl/${GKL_VERSION}/gkl-${GKL_VERSION}.jar \
          gkl.jar 612976f6f8c1fffe85bf143447dceab5410bdaa4d828f98bcc899959f7b475ba; \
    fetch $M/org/broadinstitute/gatk-native-bindings/1.0.0/gatk-native-bindings-1.0.0.jar \
          gatk-native-bindings.jar c9f5af7134cd649e00234e9c5c8310eb3e3854feb0b5e5cbe43de94412305593; \
    fetch $M/commons-logging/commons-logging/1.3.0/commons-logging-1.3.0.jar \
          commons-logging.jar 66d3c980470b99b0c511dad3dfc0ae7b265ec1fb144e96bc0253a8a175fd34d9; \
    fetch $M/commons-io/commons-io/2.16.1/commons-io-2.16.1.jar \
          commons-io.jar f41f7baacd716896447ace9758621f62c1c6b0a91d89acee488da26fc477c84f; \
    fetch $M/org/apache/commons/commons-jexl/2.1.1/commons-jexl-2.1.1.jar \
          commons-jexl.jar 03c9a9fae5da78ce52c0bf24467cc37355b7e23196dff4839e2c0ff018a01306; \
    fetch $M/org/apache/commons/commons-math3/3.5/commons-math3-3.5.jar \
          commons-math3.jar dcf9768726b82e8d7b8a928447b08bca7cfe2773d5120a63a41ec224207c1a29; \
    fetch $M/org/apache/commons/commons-compress/1.26.0/commons-compress-1.26.0.jar \
          commons-compress.jar 051aceb8bbcc62d0f5b2b8ac72c53767f9c59bfbd050151e65bef6f51c8ed9c9; \
    fetch $M/com/fulcrumgenomics/jlibdeflate/0.1.0/jlibdeflate-0.1.0.jar \
          jlibdeflate.jar c20407700e94307b80b8ed832be3fbc52370e6c3f826fe8b91b89d9db1910ef2; \
    fetch $M/org/xerial/snappy/snappy-java/1.1.10.5/snappy-java-1.1.10.5.jar \
          snappy-java.jar 0f3f1857ed33116583f480b4df5c0218836c47bfbc9c6221c0d73f356decf37b; \
    fetch $M/org/tukaani/xz/1.9/xz-1.9.jar \
          xz.jar 211b306cfc44f8f96df3a0a3ddaf75ba8c5289eed77d60d72f889bb855f535e5; \
    fetch $M/org/json/json/20231013/json-20231013.jar \
          json.jar 0f18192df289114e17aa1a0d0a7f8372cc9f5c7e4f7e39adcf8906fe714fa7d3

# commons-math3 is pinned to 3.5, not to the 3.6.1 the jmath corpus job downloads for its
# FastMath columns. 3.5 is what GATK pins (`strictly '3.5'` in its build.gradle), so it is the
# version whose numbers reach a golden downstream. The two pins disagreeing is a real gap and is
# recorded in the manifest rather than papered over.
# commons-compress is NOT optional either. `FastaReferenceWriter` counts the bytes it has
# written with `org.apache.commons.compress.utils.CountingOutputStream`, and every offset in
# the .fai it writes comes from that counter, so without the jar the class does not load at
# all: a NoClassDefFoundError rather than a degraded output. htsjdk 5.0.0 pins 1.26.0.
# The four jars after it are htsjdk 5.0.0's own runtime dependencies, at its build.gradle's
# versions, as is commons-logging 1.3.0. jlibdeflate carries libdeflate v1.25 for linux-x86_64,
# the deflater 5.0.0 picks by default; the manifest turns it off per run where a golden needs to.
ENV ORACLE_CP=/opt/oracle/jars/htsjdk.jar:/opt/oracle/jars/gkl.jar:/opt/oracle/jars/gatk-native-bindings.jar:/opt/oracle/jars/commons-logging.jar:/opt/oracle/jars/commons-io.jar:/opt/oracle/jars/commons-jexl.jar:/opt/oracle/jars/commons-math3.jar:/opt/oracle/jars/commons-compress.jar:/opt/oracle/jars/jlibdeflate.jar:/opt/oracle/jars/snappy-java.jar:/opt/oracle/jars/xz.jar:/opt/oracle/jars/json.jar:/opt/oracle/classes

# htsjdk 5.0.0 deflates with libdeflate unless told otherwise, and the property is read at class
# load. Setting it here reaches every JVM in the image (suites, probes, fuzzer, benchmark); a run
# that wants libdeflate passes -Dsamjdk.use_libdeflate=true, which the command line wins over.
ENV JAVA_TOOL_OPTIONS=-Dsamjdk.use_libdeflate=false

WORKDIR /opt/oracle
COPY OracleProbe.java .
RUN mkdir -p classes && javac -cp "$ORACLE_CP" -d classes OracleProbe.java

# Fail the build itself if the environment does not satisfy the contract, so a degraded image
# can never reach the point of producing a golden.
RUN java -cp "$ORACLE_CP" OracleProbe

WORKDIR /work
ENTRYPOINT ["/bin/sh", "-c"]
CMD ["java -cp \"$ORACLE_CP\" OracleProbe"]
