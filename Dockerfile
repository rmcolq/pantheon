FROM condaforge/mambaforge:latest AS conda

COPY environment.yml .

RUN /opt/conda/bin/mamba env create -f /environment.yml

RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

ENV PATH=/opt/conda/envs/pantheon/bin:${HOME}/edirect:$PATH

CMD ["/bin/bash"]