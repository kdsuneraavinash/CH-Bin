# Bin-X

## Installation

### Requirements

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
   ```````

### Bin-X Installation

1. Install using `setup.py`.
   ```bash
   python setup.py install
   ```
2. Copy the `config/default.ini` file and modify it to fit your configuration.
   ```bash
   cp config/default.ini config.ini
   vim config.ini
   ```
3. Run following command to run the tool.
   ```
   bin_x
   ```

## Development

1. First install [Poetry](https://python-poetry.org/docs/).
   ```bash
   curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
   ```
2. Install the requirements of the project by running
    ```bash
    poetry install
    ```
3. Run the tool via:
    ```bash
    poetry run python -m bin_x.bin_x
    ```
4. (Optional) Build the documentation via:
    ```bash
   cd docs
    poetry run sphinx-autobuild . _build/html --port 8001
    ```

Additionally, linting and type-checking are configured to this project. You may install the git-hooks for the formatters
via,

 ```bash
poetry run pre-commit install
 ```

Type-checking can be done as:

 ```bash
poetry run mypy .
 ```
