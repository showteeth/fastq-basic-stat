import os,sys,re
import argparse
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def get_qual_num(qual_li,threshold,qual_range=None):
    if not qual_range:
        filter_qual=[i for i in qual_li if i >= int(threshold)]
        filter_len=len(filter_qual)
    else:
        if len(qual_range)==1:
            filter_qual=[i for i in qual_li[int(qual_range[0]):] if i >= int(threshold)]
        else:
            filter_qual=[i for i in qual_li[int(qual_range[0]):int(qual_range[1])] if i >= int(threshold)]
        filter_len=len(filter_qual)
    return filter_len

def get_info(handle,barcode=None,UMI=None):
    read_flag=0;total_data=0;total_q20=0;total_q30=0
    barcode_q20=0;barcode_q30=0;UMI_q20=0;UMI_q30=0
    remain_q20=0;remain_q30=0
    for title, seq, qual in FastqGeneralIterator(handle):
        read_flag+=1
        qual_num=[(ord(i)-33) for i in qual]
        raw_len=len(qual_num)
        total_data+=raw_len
        q20_len=get_qual_num(qual_num,20)
        q30_len=get_qual_num(qual_num,30)
        total_q20+=q20_len
        total_q30+=q30_len
        if barcode and UMI:
            barcode_q20_len=get_qual_num(qual_num,20,barcode)
            UMI_q20_len=get_qual_num(qual_num,20,UMI)
            barcode_q30_len=get_qual_num(qual_num,30,barcode)
            UMI_q30_len=get_qual_num(qual_num,30,UMI)
            max_index=max(barcode + UMI)
            remain_q20_len=get_qual_num(qual_num,20,[max_index])
            remain_q30_len=get_qual_num(qual_num,30,[max_index])
            barcode_q20+=barcode_q20_len;UMI_q20+=UMI_q20_len;barcode_q30+=barcode_q30_len;UMI_q30+=UMI_q30_len
            remain_q20+=remain_q20_len;remain_q30+=remain_q30_len
        elif barcode:
            barcode_q20_len=get_qual_num(qual_num,20,barcode)
            barcode_q30_len=get_qual_num(qual_num,30,barcode)   
            max_index=max(barcode)
            remain_q20_len=get_qual_num(qual_num,20,[max_index])
            remain_q30_len=get_qual_num(qual_num,30,[max_index])
            barcode_q20+=barcode_q20_len;barcode_q30+=barcode_q30_len
            remain_q20+=remain_q20_len;remain_q30+=remain_q30_len
        elif UMI:
            UMI_q20_len=get_qual_num(qual_num,20,UMI)
            UMI_q30_len=get_qual_num(qual_num,30,UMI)
            max_index=max(UMI)
            remain_q20_len=get_qual_num(qual_num,20,[max_index])
            remain_q30_len=get_qual_num(qual_num,30,[max_index])
            UMI_q20+=UMI_q20_len;UMI_q30+=UMI_q30_len
            remain_q20+=remain_q20_len;remain_q30+=remain_q30_len
    if barcode and UMI:
        return read_flag,total_data,total_q20,total_q30,barcode_q20,barcode_q30,UMI_q20,UMI_q30,remain_q20,remain_q30
    elif barcode:
        return read_flag,total_data,total_q20,total_q30,barcode_q20,barcode_q30,remain_q20,remain_q30
    elif UMI:
        return read_flag,total_data,total_q20,total_q30,UMI_q20,UMI_q30,remain_q20,remain_q30
    else:
        return read_flag,total_data,total_q20,total_q30

def get_qualPercent_data(args):
    if args.file_type=='fastq':
        with open(args.fastq, "r") as handle:
            all_info=get_info(handle,barcode=args.barcode,UMI=args.umi)
    else:
        with gzip.open(args.fastq, "rt") as handle:
            all_info=get_info(handle,barcode=args.barcode,UMI=args.umi)
    return all_info

def main():
    parser = argparse.ArgumentParser(
        description="""stat raw fastq data info, including total data、(barcode、UMI)Q20 percent、Q30 percent \n
        python fastq_stat.py -f test.fastq -t fastq [-b 0 16 -u 16 25]""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--fastq',
        '-f',
        type=str,
        help='full path of fastq file')
    parser.add_argument('--file_type',
                        '-t',
                        type=str,
                        default="gzip",
                        choices=["gzip","fastq"],
                        help="""fastq file type, gzip of fastq""")
    parser.add_argument('--barcode',
                        '-b',
                        type=int,
                        nargs='*',
                        help='barcode index list')
    parser.add_argument('--umi', '-u', type=int,nargs='*',
                        help='umi index list')

    args = parser.parse_args()

    all_info=get_qualPercent_data(args)
    print ('\t'.join([str(i) for i in all_info]))

if __name__ == '__main__':
    main()