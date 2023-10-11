FROM databio/pepatac:latest
MAINTAINER Aurelie Gabriel <aurelie.gabriel@unil.ch>

#RUN wget --no-check-certificate https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
#RUN tar xzvf sratoolkit.2.11.3-ubuntu64.tar.gz
#ENV PATH $PATH:/home/sratoolkit.2.11.3-ubuntu64/bin
RUN wget --no-check-certificate https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz
RUN tar xzvf sratoolkit.3.0.0-ubuntu64.tar.gz
ENV PATH $PATH:/home/sratoolkit.3.0.0-ubuntu64/bin


# RUN cp -r sratoolkit.3.0.0-ubuntu64/bin/* /usr/local/bin

RUN mkdir /root/.ncbi/
RUN echo $'## auto-generated configuration file - DO NOT EDIT ## \n/LIBS/GUID = "50b75afa-afc3-43ef-8cc0-2177e53eddac"\n/config/default = "false" \n/repository/user/ad/public/apps/file/volumes/flatAd = "."\n/repository/user/ad/public/apps/refseq/volumes/refseqAd = "."\n/repository/user/ad/public/apps/sra/volumes/sraAd = "."\n/repository/user/ad/public/apps/sraPileup/volumes/ad = "."\n/repository/user/ad/public/apps/sraRealign/volumes/ad = "."\n/repository/user/ad/public/apps/wgs/volumes/wgsAd = "."\n/repository/user/ad/public/root = "."\n/repository/user/default-path = "/root/ncbi"  ' > /root/.ncbi/user-settings.mkfg

RUN git clone https://github.com/databio/pepatac.git
ENV PATH $PATH:/home/pepatac/pipelines/


RUN pip install multiqc
RUN apt-get install --no-install-recommends -y build-essential
RUN apt-get update
#RUN apt-get install --no-install-recommends -y python3-dev python3-pip python3-numpy python3-setuptools python3-matplotlib
#RUN wget https://github.com/nboley/idr/archive/refs/tags/2.0.3.zip
#RUN unzip 2.0.3.zip
#RUN cd /home/idr-2.0.3 && python3 setup.py install


#RUN git clone https://github.com/rhysnewell/ChIP-R.git
#RUN cd ChIP-R && python3 setup.py install

# Install featureCounts
RUN apt-get -y install subread

#RUN  sed -i 's/:2:30:10/:2:24:10/g' pepatac/pipelines/pepatac.py
#RUN chmod +x pepatac/pipelines/pepatac.py

RUN pip3 install biopython
RUN apt-get update
RUN pip3 install sinto
RUN pip3 install pyega3

RUN apt install --no-install-recommends -y libgeos-dev
RUN Rscript -e "install.packages(c('Seurat','Signac'))"
RUN Rscript -e "BiocManager::install('liftOver')"
#RUN Rscript -e "BiocManager::install(version = '3.14')"

RUN pip3 install pyega3 --upgrade
RUN Rscript -e "BiocManager::install('edgeR')"
RUN Rscript -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38','BSgenome.Hsapiens.UCSC.hg19','BSgenome.Mmusculus.UCSC.mm10'))"
