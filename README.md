# sourmash_plugin_extract_errors

## Motivation

Illumina's sequencing errors are unique per sample. If you are working on a pool of samples and you want to remove only the errors that are unique to each sample, you can use this plugin to create a sourmash signature of these errors then subtract it from the original signatures.


## Description

The plugin is backed with C++ for performance and it supports multi-threading. To use it, you can run the following command:

```
sourmash scripts extract_errors --ksize 31 --num-threads 16 --out errors.sig *sig
```

## Usage

```
usage:  extract_errors [-h] [-q] [-d] [--ksize KSIZE] [--num-threads NUM_THREADS] [--out OUT] sig_paths [sig_paths ...]

positional arguments:
  sig_paths

options:
  -h, --help            show this help message and exit
  -q, --quiet           suppress non-error output
  -d, --debug           provide debugging output
  --ksize KSIZE
  --num-threads NUM_THREADS
  --out OUT
```


## Installation

This plugin has a C++ extension, so you will need a C++ compiler to install it. On Linux, you can install the `build-essential` package to get a C++ compiler. On macOS, you can install the Xcode command line tools with `xcode-select --install`.

```
pip install git+https://github.com/mr-eyes/sourmash_plugin_extract_errors
```


## Technical details

We use `simdjson` for parsing the json file, and we use the `parallel_hashmap` to count the kmers in parallel.