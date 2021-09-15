# CH-Bin

## Installation

First clone the CH-Bin repository to a local directory. Note that CH-Bin only supports linux.

```bash
git clone https://github.com/kdsuneraavinash/CH-Bin
```

### Dependencies

You will need python and build requirements (and optionally venv). In ubuntu 20.04, you can install them by,

```bash
sudo apt-get install build-essential
sudo apt-get install python3 python-is-python3
sudo apt-get install python3-dev
sudo apt-get install python3.8-venv
```

Additionally, [FragGeneScan](https://sourceforge.net/projects/fraggenescan), [HMMER](http://hmmer.org/)
and [kmer-counter](https://github.com/alexpreynolds/kmer-counter) tools are required. If you have them already installed, you may provide their paths in the configuration. Otherwise, follow the below steps to install them manually. Following commands will install the required tools in the `tools` directory.

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

If you installed them in a different directory than `tools`, you need to update the configuration. Edit `config/default.ini` as follows, (If you followed the default installation commands provided above, you do not need to change the configuration.)

```ini
[COMMANDS]
FragGeneScan = <FragGeneScan Path>
Hmmer = <Hmmer Path>
KMerCounter = <kmer-counter Path>
```

### CH-Bin Installation

1. Install using `setup.py`. (Recommended installing in a virtual environment)
    ```bash
    python -m venv .venv        # Optional: creating venv
    source .venv/bin/activate   # Optional: activating venv
    python setup.py install
    ```
2. Then run following command to run the tool.
    ```bash
    ch_bin --contigs [CONTIG_FASTA] --coverages [COVERAGE_TSV] --out [OUT_DIR]
    ```

For example, to bin the sample dataset, run the following command.
 ```bash
 ch_bin --contigs test_data/five-genomes-contigs.fasta --coverages test_data/five-genomes-abundance.abund --out out
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
