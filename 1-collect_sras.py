from tbtools import collector, ISprofiler, tools

#To get SRR_Acc_List.txt file, navigate to sra database in NCBI by searching on bacteria name, then click the number beside the desire bacteria, then here are the tuberculosis available sras : https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=1763&lvl=3&p=sra&lin=f&keep=1&srchmode=1&unlock , click on a number in gray->at least equal to 20-> then "Send results to run selector" -> and then choose accession list


data=open('data/SRR_Acc_List.txt').read().split('\n')
outfolder='sequences'

#Download the SRAs in data from NCBI
c=collector.Collector(sras=data, outdir=outfolder)     #we can use bioprojects=[list of bio projects] or biosamples=[list of biosamples], the last parameter is outdir='sequences' by default->change it to sequences/sra         

c.get()

