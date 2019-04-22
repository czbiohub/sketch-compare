FROM continuumio/anaconda3

RUN conda install -c bioconda ska
RUN pip install --upgrade --force-reinstall pip==9.0.3
RUN pip install awscli --disable-pip-version-check
RUN pip install boto3
RUN pip install --upgrade pip
ADD prep_files.py /usr/local/bin/prep_files.py
RUN chmod 777 /usr/local/bin/prep_files.py

WORKDIR /scratch