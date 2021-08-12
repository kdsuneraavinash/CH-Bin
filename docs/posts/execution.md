# Execution

The tool can be run via `bin_x` command. Alternatively you can execute it directly using `python -m bin_x.bin_x`.

```bash
$ python -m bin_x.bin_x run --help
Usage: python -m bin_x.bin_x run [OPTIONS]

Options:
  -i, --contigs PATH    The contig file to perform the binning operation.
                        [required]
  -c, --coverages PATH  The tab-seperated file with abundance data.
                        [required]
  -s, --config PATH     The configuration file path.
  -o, --out PATH        The output directory for the tool.
  --help                Show this message and exit.
```
