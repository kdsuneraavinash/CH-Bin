# Setup

## Custom configuration

If you followed the installation document, the requirements should be installed in the `tools` directory. If so you can
skip this section and run the tool directly. If the tools were installed in a different manner, or you want ot change
some internal configuration, you might want to provide a custom configuration.

### Configuration file format

Bin-X uses a customizable `ini` file format to store the configuration. The default configuration can be found
in `config/default.ini`.

```ini
[COMMANDS]
FragGeneScan = PATH OF run_FragGeneScan.pl
HmmSearch = PATH OF hmmsearch
KMerCounter = PATH OF kmer-counter
Seq2Vec = PATH OF seq2vec

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
```

If the requirement installations were done differently, you might want to change the `COMMANDS` section.

### Commands/Resource Configuration

For each tool provide the command that can be used to run the specified tool. For example, for `FragGeneScan`, provide
the `run_FragGeneScan.pl` script location. Note that all the related files should be executable.
(In `FragGeneScan` case, both `FragGeneScan` and `run_FragGeneScan.pl` should be executable)

The resource section describes the resource locations. Generally, you do not need to change this section.

### Parameters Configuration

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

### Providing Configuration

When running `bin_x` you can provide the custom configuration file via, `-s` or `--config` parameter.

## Using seq2vec

`seq2vec` is a fast kmer-counter tool. But due to installation issues, by default, bin-x uses kmer-counter. However,
there will be a performance improvement if `seq2vec` is used.

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

4. Finally, set the configuration parameters as follows and provide the modified configuration file when running `bin_x`
   .

   ```ini
   [COMMANDS]
   Seq2Vec = tools/seq2vec/seq2vec
   
   [PARAMETERS]
   KmerCounterTool = seq2vec
   ```
