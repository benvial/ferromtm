FROM benvial/pyonelab:latest



RUN git clone https://github.com/benvial/ferromtm.git

WORKDIR /ferromtm

RUN make env

RUN echo "conda activate ferromtm" >> ~/.bashrc

RUN make testenv && make req


ENV APPNAME ferromtm

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
