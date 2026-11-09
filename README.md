# SOWA

Reorder DNA reads files for a better compression by all-purpose compressors (gzip, zstd..)

## Installation
```
git clone https://github.com/girunivlille/sowa
cd sowa
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Usage
```
RUSTFLAGS="-C target-cpu=native" cargo run --release -- --reads-in input_file.fasta --reads-out output_file.fasta
```
