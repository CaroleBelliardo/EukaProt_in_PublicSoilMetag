head ../../MetaNema_LRmg/10Metag/kraken_reads/aubergine_I_oct2__hifi_reads.fasta.krak -n1000 > aubergine_I_oct2__hifi_reads.fasta.krak.head1000
awk -F'\t' 'NR==FNR {h[] = ; next} {print bash,h[]}' /database/hub/TAXONOMY/online/taxonomy_lineage.txt aubergine_I_oct2__hifi_reads.fasta.krak.head1000 > aubergine_I_oct2__hifi_reads.fasta.krak.head1000_lineage 
