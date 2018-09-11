# biotools

some small scripts or tools for bioinfomatics

## bamdst_parse

parse bamdst output file.

### Usage

```bash
$ ./bamdst_parse.py 
usage: bamdst_parse.py [-h] -s SORTEDBAM -r RMDUPBAM -b BED [-d BAMDST]
                       [-o TEMPDIR] [-f FLANK] [-m MAPQ] [-u UNCOVER] [-n]
                       [-p SEP]

optional arguments:
  -h, --help            show this help message and exit
  -s SORTEDBAM, --sortedbam SORTEDBAM
                        sorted bam file
  -r RMDUPBAM, --rmdupbam RMDUPBAM
                        rmdup bam file
  -b BED, --bed BED     target bed file
  -d BAMDST, --bamdst BAMDST
                        path to bamdst software.
                        default:/GPFS01/softwares/bamdst/bamdst
  -o TEMPDIR, --tempdir TEMPDIR
                        path to bamdst output dir. default=./
  -f FLANK, --flank FLANK
                        base pairs in each direction. default:100
  -m MAPQ, --mapq MAPQ  map quality cutoff value, greater or equal to the
                        value will be count. default:20
  -u UNCOVER, --uncover UNCOVER
                        region will included in uncover file if below it.
                        default:20
  -n, --noheader        print header or not, default:with header
  -p SEP, --sep SEP     delimiter, default:,
```

### Example

```bash
./bamdst_parse.py -s test.sorted.bam -r test.sorted.rmdup.bam -b target.bed > output.csv

./bamdst_parse.py -s test.sorted.bam -r test.sorted.rmdup.bam -b target.bed -p \t > output.tsv
```

### Output
