# Dockerfile
FROM alpine:3.20
COPY smina.static /usr/local/bin/smina
RUN chmod +x /usr/local/bin/smina
WORKDIR /work
ENTRYPOINT ["smina"]