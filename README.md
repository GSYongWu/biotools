# biotools

some small scripts or tools for bioinfomatics

## bamdst_parse

parse bamdst output file.

### Usage

```shell
$ ./bamdst_parse.py
usage: bamdst_parse.py [-h] -s SORTEDBAM -r RMDUPBAM -b BED [-d BAMDST]
                       [-t TEMPDIR] [-o OUTPUT] [-f FLANK] [-m MAPQ]
                       [-u UNCOVER] [-n] [-p SEP]

optional arguments:
  -h, --help            show this help message and exit
  -s SORTEDBAM, --sortedbam SORTEDBAM
                        sorted bam file
  -r RMDUPBAM, --rmdupbam RMDUPBAM
                        rmdup bam file
  -b BED, --bed BED     target bed file
  -d BAMDST, --bamdst BAMDST
                        path to bamdst software
                        [/GPFS01/softwares/bamdst/bamdst]
  -t TEMPDIR, --tempdir TEMPDIR
                        path to bamdst output dir [./]
  -o OUTPUT, --output OUTPUT
                        output file name [stdout]
  -f FLANK, --flank FLANK
                        base pairs in each direction [100]
  -m MAPQ, --mapq MAPQ  map quality cutoff value, greater or equal to the
                        value will be count [20]
  -u UNCOVER, --uncover UNCOVER
                        region will included in uncover file if below it [20]
  -n, --noheader        print no header [False]
  -p SEP, --sep SEP     delimiter, default:,

```

### Example

```bash
./bamdst_parse.py -s test.sorted.bam -r test.sorted.rmdup.bam -b target.bed > output.csv

./bamdst_parse.py -s test.sorted.bam -r test.sorted.rmdup.bam -b target.bed -p \t > output.tsv
```

### Output

| Title                     | Description                             |
| ------------------------- | --------------------------------------- |
| #SAMPLE                   | Sample name                             |
| CLEAN_READS               | Clean reads count                       |
| CLEAN_BASES               | Clean base count(MB)                    |
| INSERT_SIZE               | Insert size                             |
| MAPQ20(%)                 | MapQ>=20 reads pct                      |
| DUPLICATE(%)              | Duplicate pct                           |
| DUPLICATE_TARGET(%)       | duplicate(target region)                |
| ON_TARGET_CORE(%)         | Ontarget in core region(base)           |
| ON_TARGET_EXT(%)          | Ontarget in extend 100bp region(base)   |
| ON_TARGET_READS_CORE(%)   | Ontarget reads pct                      |
| ON_TARGET_READS_EXT(%)    | Ontarget reads flank 100bp pct          |
| RATIO_OF_MAPPED(%)        | Ratio of mapped reads                   |
| MEAN_DEPTH                | Mean depth                              |
| MEDIAN_DEPTH              | Median depth                            |
| 1X_COVERAGE(%)            | Coverage >=1X                           |
| 20X_COVERAGE(%)           | Coverage >=20X                          |
| 50X_COVERAGE(%)           | Coverage >=50X                          |
| 100X_COVERAGE(%)          | Coverage >=100X                         |
| 200X_COVERAGE(%)          | Coverage >=200X                         |
| 500X_COVERAGE(%)          | Coverage >=500X                         |
| 10%_COVERAGE(%)           | Coverage >=(Mean depth)*10%             |
| 20%_COVERAGE(%)           | Coverage >=(Mean depth)*20%             |
| 50%_COVERAGE(%)           | Coverage >=(Mean depth)*50%             |
| MEAN_DEPTH_DEDUP          | Mean dedup depth                        |
| MEDIAN_DEPTH_DEDUP        | Median dedup depth                      |
| 1X_COVERAGE_DEDUP(%)      | Dedup coverage >=1X                     |
| 20X_COVERAGE_DEDUP(%)     | Dedup coverage >=20X                    |
| 50X_COVERAGE_DEDUP(%)     | Dedup coverage >=50X                    |
| 100X_COVERAGE_DEDUP(%)    | Dedup coverage >=100X                   |
| 200X_COVERAGE_DEDUP(%)    | Dedup coverage >=200X                   |
| 500X_COVERAGE_DEDUP(%)    | Dedup coverage >=500X                   |
| 10%MEAN_COVERAGE_DEDUP(%) | Dedup coverage >=(Mean dedup depth)*10% |
| 20%MEAN_COVERAGE_DEDUP(%) | Dedup coverage >=(Mean dedup depth)*20% |
| 50%MEAN_COVERAGE_DEDUP(%) | Dedup coverage >=(Mean dedup depth)*50% |
| CV_SCORE                  | cv score                                |