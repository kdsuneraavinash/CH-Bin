# CH-Bin

## Installation

### Requirements

You will need python and build requirements (and optionally venv). In ubuntu 20.04, you can install them by,
```bash
sudo apt-get install build-essential
sudo apt-get install python3 python-is-python3
sudo apt-get install python3-dev
sudo apt-get install python3.8-venv
```

1. Install [FragGeneScan](https://sourceforge.net/projects/fraggenescan).
    ```bash
    wget -O tools/FragGeneScan1.31.tar.gz https://sourceforge.net/projects/fraggenescan/files/FragGeneScan1.31.tar.gz
    tar -xzf tools/FragGeneScan1.31.tar.gz -C tools
    cd tools/FragGeneScan1.31 && make clean && make fgs && cd ../..
    ./tools/FragGeneScan1.31/run_FragGeneScan.pl
    ```
2. Install [HMMER](http://hmmer.org/).
    ```bash
    wget -O tools/hmmer.tar.gz http://eddylab.org/software/hmmer/hmmer.tar.gz
    tar -xzf tools/hmmer.tar.gz -C tools
    cd tools/hmmer-3.3.2 && ./configure && make && cd ../..
    ./tools/hmmer-3.3.2/src/hmmalign -h
    ```
3. Install [kmer-counter](https://github.com/alexpreynolds/kmer-counter).
    ```bash
    git clone https://github.com/alexpreynolds/kmer-counter tools/kmer-counter
    cd tools/kmer-counter/ && make && cd ../..
    ./tools/kmer-counter/kmer-counter --help
    ```
4. Get marker-gene hmm file.
    ```bash
    wget -O tools/marker.hmm -nc https://raw.githubusercontent.com/sufforest/SolidBin/4c9b9ea7b8d8a0df1b772669872b69006c490e67/auxiliary/marker.hmm
    head ./tools/marker.hmm
    ```

### CH-Bin Installation

1. Install using `setup.py`. (Recommended installing in a virtual environment)
    ```bash
    python -m venv .venv
    source .venv/bin/activate
    python setup.py install
    ```
2. Run following command to run the tool.
    ```
    ch_bin --contig [CONTIG_FASTA] --coverages [COVERAGE_TSV] --out [OUT_DIR]
    ```

## Development

1. Install the requirements of the project by running
    ```bash
    pip install -r requirements.txt
    ```
2. Run the tool via:
    ```bash
    python -m ch_bin.ch_bin
    ```
3. (Optional) Build the documentation via:
    ```bash
    cd docs
    sphinx-autobuild . _build/html --port 8001
    ```

Additionally, linting and type-checking are configured to this project. You may install the git-hooks for the formatters
via,

 ```bash
pre-commit install
 ```

Type-checking can be done as:

 ```bash
mypy .
 ```
