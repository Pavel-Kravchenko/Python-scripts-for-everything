cat p1.fasta > p.fasta
for x in 1 2 3 4 5 6 7 8 

do 
   let "y = $x + 1"
   msbar -count 7 -point 1 -block 0 p$x.fasta p$y.fasta
   descseq -name "Evolution_rate_is:"p$y p$y.fasta p$y.fasta 

done
cat p[2-9].fasta >> p.fasta