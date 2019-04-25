FROM continuumio/anaconda3

RUN conda install -c bioconda ska
RUN pip install --upgrade --force-reinstall pip==9.0.3
RUN pip install awscli --disable-pip-version-check
RUN pip install boto3
RUN pip install --upgrade pip
ADD create_sketch.py /usr/local/bin/create_sketch.py
RUN chmod 777 /usr/local/bin/create_sketch.py

WORKDIR /scratch