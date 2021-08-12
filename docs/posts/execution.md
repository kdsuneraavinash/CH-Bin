# Execution

The tool can be run via `bin_x` command. Alternatively you can execute it directly using `python -m bin_x.bin_x`.

```bash
$ python -m bin_x.bin_x --help
Usage: python -m bin_x.bin_x [OPTIONS]

Options:
  -i, --contigs PATH       The contig file to perform the binning operation.
                           [required]
  -c, --coverages PATH     The tab-seperated file with abundance data.
                           [required]
  -g, --ground_truth PATH  The ground truth CSV.  [required]
  -s, --config PATH        The configuration file path.
  -o, --out TEXT           The output directory for the tool.
  -t, --n_iter INTEGER     Number of Iterations to run.
  --help                   Show this message and exit.
```
