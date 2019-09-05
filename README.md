## fastq stat 
* reads num 、data volume、Q20、Q30
* Q20、Q30 of barcode and UMI 

## usage
### help info
```bash
python fastq_stat.py --help
usage: fastq_stat.py [-h] [--fastq FASTQ] [--file_type {gzip,fastq}]
                     [--barcode [BARCODE [BARCODE ...]]]
                     [--umi [UMI [UMI ...]]]

stat raw fastq data info, including total data、(barcode、UMI)Q20 percent、Q30 percent 

        python fastq_stat.py -f test.fastq -t fastq [-b 0 16 -u 16 25]

optional arguments:
  -h, --help            show this help message and exit
  --fastq FASTQ, -f FASTQ
                        full path of fastq file
  --file_type {gzip,fastq}, -t {gzip,fastq}
                        fastq file type, gzip of fastq
  --barcode [BARCODE [BARCODE ...]], -b [BARCODE [BARCODE ...]]
                        barcode index list
  --umi [UMI [UMI ...]], -u [UMI [UMI ...]]
                        umi index list
```

### detailed
```bash
python fastq_stat.py -f test.fastq.gz
    39561   5103369 4493270 3856336
```


## plan 
- [ ] parallel
