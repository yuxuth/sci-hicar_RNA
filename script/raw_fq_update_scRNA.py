import gzip


def update_fastq(r1,r2, out_r1,out_r2 ): ## process two files

    f_r1 = gzip.open(r1, 'rt')
    f_r2 = gzip.open(r2, 'rt')
    f_out_r1 = open(out_r1, 'w')
    f_out_r2 = open(out_r2, 'w')

    while True:
        cur_r1_name = f_r1.readline().strip()
        cur_r1_read = f_r1.readline().strip()
        cur_r1_plus = f_r1.readline().strip()
        cur_r1_qual = f_r1.readline().strip()
        
        cur_r2_name = f_r2.readline().strip()
        cur_r2_read = f_r2.readline().strip()
        cur_r2_plus = f_r2.readline().strip()
        cur_r2_qual = f_r2.readline().strip()
    
        if cur_r1_name == "" : break
        
        tn5_index_i5 = cur_r1_read[0:6] ## 6bp tn5 index sequence
        read1_umi = cur_r1_read[6:12]
        tn5_index_i7 = cur_r2_read[0:6] ## 6bp tn5 index sequence
        read2_umi = cur_r2_read[6:12]

        i5_i7 = (cur_r2_read[100:106])  ## current_110bp_formate, extract 100-106

        new_read2 = tn5_index_i5 + tn5_index_i7 + i5_i7  + read1_umi + read2_umi  # 6+6+6 18barcode + 6+6 12umi
        new_read2_qual = cur_r1_qual[0:6] + cur_r2_qual[0:6] + cur_r2_qual[100:106] + cur_r1_qual[6:12] +cur_r2_qual[6:12]

        
        cur_r1_read = cur_r1_read[37:] ## trim first 37 bp, remain is the RNA reads
        cur_r1_qual = cur_r1_qual[37:] # trim first 37 bp, remain is the RNA reads

        cur_r2_read = new_read2 
        cur_r2_qual = new_read2_qual
        
        cur_r2_name = cur_r2_name.replace("/2","/1")
                
        
        f_out_r1.write(cur_r1_name+"\n")
        f_out_r1.write(cur_r1_read+"\n")
        f_out_r1.write(cur_r1_plus+"\n")
        f_out_r1.write(cur_r1_qual+"\n")     
    
        f_out_r2.write(cur_r2_name+"\n")
        f_out_r2.write(cur_r2_read+"\n")
        f_out_r2.write(cur_r2_plus+"\n")
        f_out_r2.write(cur_r2_qual+"\n")
        
    f_r1.close()
    f_r2.close()
    f_out_r1.close()
    f_out_r2.close()
    
r1 = str(snakemake.input['r1'])
r2 = str(snakemake.input['r2'])
out_r1 = str(snakemake.output['r1'])
out_r2 = str(snakemake.output['r2'])

update_fastq(r1,r2, out_r1,out_r2 )
