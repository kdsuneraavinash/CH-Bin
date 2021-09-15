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

## Required File Formats

### Contig File

Contig file should be a multi-FASTA file containing the genomes in following format.

```fasta
>NODE_509_length_56_cov_70.000000
AAGGCTCTTCAGGAATAAGAGTGTAACCACCTGAAACCAACACCCCGATTCCCGGG
>NODE_510_length_56_cov_68.000000
CCAGCAGAACCCCTGGTCCTGCTAACTCGGTGTCCACTACCCGGGGTGAACCTCAC
```

### Coverage File

Coverage file should be a tab-separated file with the following format.

```tsv
NODE_1_length_1189502_cov_16.379288	16.379288
NODE_2_length_1127036_cov_16.549343	16.549343
NODE_3_length_1009819_cov_16.436396	16.436396
NODE_4_length_861895_cov_21.063754	21.063754
NODE_5_length_737013_cov_20.834031	20.834031
NODE_6_length_659011_cov_21.171279	21.171279
```

The first column should contain the contig id, as given in the contigs FASTA. Each column after that should refer to the coverages from each sample. (There should be at-least one sample) The example file above shows a coverage file with data taken from one sample.

### Coverage File

If you followed the installation document, the requirements should be installed in the `tools` directory. If so you can
skip this section and run the tool directly. If the tools were installed in a different manner, or you want to change
some internal configuration, you might want to provide a custom configuration.

#### Configuration file format

CH-Bin uses a customizable `ini` file format to store the configuration. The default configuration can be found
in `config/default.ini`.

```ini
[COMMANDS]
FragGeneScan = DIR OF run_FragGeneScan.pl
Hmmer = DIR OF Hmmer
KMerCounter = DIR OF kmer-counter
Seq2Vec = DIR OF seq2vec

[RESOURCES]
MarkersHmm = PATH OF marker.hmm

[PARAMETERS]
KmerK = INTEGER
KmerCounterTool = EITHER kmer_counter OR seq2vec
ContigLengthFilterBp = INTEGER
ScmCoverageThreshold = FLOAT BETWEEN 0 AND 1
ScmSelectPercentile = FLOAT BETWEEN 0 AND 1
SeedContigSplitLengthBp = INTEGER
AlgoNumNeighbors = INTEGER
AlgoMaxIterations = INTEGER
AlgoDistanceMetric = convex
AlgoQpSolver = EITHER quadprog OR cvxopt
InMemDistMatrix = EITHER yes OR no
```

If the requirement installations were done differently, you might want to change the `COMMANDS` section.

#### Commands/Resource Configuration

For each tool provide the command that can be used to run the specified tool. For example, for `FragGeneScan`, provide
the `run_FragGeneScan.pl` script location. Note that all the related files should be executable.
(In `FragGeneScan` case, both `FragGeneScan` and `run_FragGeneScan.pl` should be executable)

The resource section describes the resource locations. Generally, you do not need to change this section.

#### Parameters Configuration

In the parameters section, you can adjust the default tool settings. Following table shows available properties.

| Parameter               | Description                                                  |
| ----------------------- | ------------------------------------------------------------ |
| KmerK                   | K value of the kmers to count.                               |
| KmerCounterTool         | Kmer counter tool to use. (kmer_counter/seq2vec)             |
| ContigLengthFilterBp    | Threshold to filter the short contigs.                       |
| ScmCoverageThreshold    | Threshold for a hit to be considered for the seed frequency distribution. |
| ScmSelectPercentile     | Percentile to use for selecting the number of seeds. For example, 0.5 will take the median number of seeds. |
| SeedContigSplitLengthBp | Length to split the seed contigs.                            |
| AlgoNumNeighbors        | Number of neighbors to consider for polytope.                |
| AlgoMaxIterations       | Number of maximum iterations to perform.                     |
| AlgoDistanceMetric      | Polytope distance matrix (convex/affine)                     |
| AlgoQpSolver            | Quadratic programming problem solver. (quadprog/cvxopt)      |
| InMemDistMatrix         | Whether to use in-memory or in-disk distance matrix. (yes/no)    |

#### Providing Configuration

When running `ch_bin` you can provide the custom configuration file via, `-s` or `--config` parameter.

#### Using seq2vec

`seq2vec` is a fast kmer-counter tool. But due to difficulty in installing, by default, CH-Bin uses kmer-counter. However, there will be a performance improvement if `seq2vec` is used.

1. First install boost (1.72+)
    ```bash
    sudo apt-get install libboost-all-dev
    ```

2. In case that the version installed by above command is less than 1.72, you may want
   to [install from source](https://www.boost.org/doc/libs/1_76_0/more/getting_started/unix-variants.html#id20).
    ```bash
    sudo apt-get install build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev
    wget -O tools/boost_1_76_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.76.0/boost_1_76_0.tar.gz/download
    cd tools && tar xzvf boost_1_76_0.tar.gz
    cd tools/boost_1_76_0/ && ./bootstrap.sh --prefix=/usr/ && chmod +x b2 && ./b2 && sudo ./b2 install && cd ../..
    ```

3. Then download and build the `seq2vec` project.
    ```bash
    git clone https://github.com/anuradhawick/seq2vec.git tools/seq2vec
    cd tools/seq2vec/ && mkdir build && cmake . && make -j8 && cd ../..
    ./tools/seq2vec/seq2vec --help
    ```

4. Finally, set the configuration parameters as follows and provide the modified configuration file when running `ch_bin`
   .

   ```ini
   [COMMANDS]
   Seq2Vec = tools/seq2vec/seq2vec

   [PARAMETERS]
   KmerCounterTool = seq2vec
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
