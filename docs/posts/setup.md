# Setup

If you followed the installation document, the requirements should be installed in the `tools` directory. If so you can
skip this section and run the tool directly. If the tools were installed in a different manner, or you want ot change
some internal configuration, you might want to provide a custom configuration.

## Configuration file format

Bin-X uses a customizable `ini` file format to store the configuration. The default configuration can be found
in `config/default.ini`.

```ini
[COMMANDS]
FragGeneScan = tools/FragGeneScan1.31/run_FragGeneScan.pl
HmmSearch = tools/hmmer-3.3.2/src/hmmsearch
KMerCounter = tools/kmer-counter/kmer-counter

[RESOURCES]
MarkersHmm = tools/marker.hmm

[PARAMETERS]
```

If the requirement installations were done differently, you might want to change the `COMMANDS` section.

## Commands Configuration

For each tool provide the command that can be used to run the specified tool. For example, for `FragGeneScan`, provide
the `run_FragGeneScan.pl` script location. Note that all the related files should be executable.
(In `FragGeneScan` case, both `FragGeneScan` and `run_FragGeneScan.pl` should be executable)

## Resources Configuration

This section describes the resource locations. Generally, you do not need to change this section.

## Parameters Configuration

In the parameters section, you can adjust the default tool settings. Following table shows available properties.

> TODO

## Providing Configuration

When running `bin_x` you can provide the custom configuration file via, `-s` or `--config` parameter.
