# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/ndexkinomeloader
COPY ./ /tmp/ndexkinomeloader/
RUN pip install /tmp/ndexkinomeloader

RUN rm -rf /tmp/ndexkinomeloader

RUN mkdir /ndexkinomeloader
WORKDIR /ndexkinomeloader

CMD /bin/echo "To use run: ndexloadkinome.py"
