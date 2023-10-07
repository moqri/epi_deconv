wget -nc -i ../../ftps
for f in beta/*Blood*.hg38.beta; do  fn="${f%%.*}";  echo ${fn:5:-10}.bed;  /oak/stanford/scg/lab_mpsnyder/moqri/soft/wgbs_tools/src/python/wgbs_tools.py beta2bed --genome hg38 $f > bed/${fn:5:-10}.bed; done

for f in beta/*-WBC-WGBS-Rep1.beta; do  fn="${f%%.*}";  echo ${fn:5:-10}.bed;  /oak/stanford/scg/lab_mpsnyder/moqri/soft/wgbs_tools/src/python/wgbs_tools.py beta2bed --genome hg19 $f > bed/${fn:5:-10}.bed; done

for f in bed/*-WBC.bed
do 
 /labs/mpsnyder/moqri/soft/ucsc/liftOver $f /labs/mpsnyder/moqri/soft/ucsc/hg19ToHg38.over.chain.gz $f tmp
done
#mask
#for f in bed/*.bed; do  fn="${f%.*}";  echo ${fn:4};  bedtools subtract -a $f -b mask.hg38 > mask/${fn:4}; done
#bedtools subtract -a neo/meth_count.csv.bed -b mask.hg38 > mask/neo

for f in bed/*.bed; do  echo $f;  awk '{print $1,$2,"+","CpG",$4/$5,$5}' $f > meth/${f:5}; done

dnmtools merge meth/*_Blood-T-CD3.bed -o merge/cd3
dnmtools merge meth/*_Blood-T-CD4.bed -o merge/cd4
dnmtools merge meth/*_Blood-T-CD8.bed -o merge/cd8
dnmtools merge meth/*_Blood-NK.bed -o merge/nk
dnmtools merge meth/*_Blood-Monocytes.bed -o merge/mono
dnmtools merge meth/*_Blood-Granulocytes.bed -o merge/gran
dnmtools merge meth/*_Blood-B.bed -o merge/b

bedops -m cd3 cd4 cd8 nk mono gran b > merge
declare -a l=(cd3 cd4 cd8 nk mono gran b merge)

for f in ${l[@]}; do time dnmtools hmr-rep merge/$f -o hmr/$f; done
for f in ${l[@]}; do time bedtools subtract -a hmr/$f -b black -A > hmr/$f.sub; done
for f in ${l[@]}; do LC_ALL=C sort -k 1,1 -k 2,2n hmr/$f.sub -o hmr/$f.sub; done
dnmtools selectsites hmr/merge.sub merge/gran | awk {'print $1,$2'} > erm

for f in ${l[@]}; do time LC_ALL=C sort -k 1,1 -k 2,2n merge/$f -o merge/$f; done
for f in ${l[@]}; do time dnmtools roi hmr/$f merge/$f > hmr/$f.m; done
for f in ${l[@]}; do dnmtools roi hmr/$f merge/wbc > hmr/$f.b; done
for f in ${l[@]}; do dnmtools roi hmr/$f merge/neo > hmr/$f.n; done
for f in ${l[@]}; do time paste hmr/$f hmr/$f.n hmr/$f.b hmr/$f.m | awk '{printf "%s,%.0f,%.0f,%.0f,%.2f,%.2f,%.2f\n", $1,$2,$3,$5,$11,$17,$23}' > hmr/$f.a; done

time bedops -m cd3 cd4 cd8 nk mono gran b > merge
bedtools subtract -a merge -b ../black -A > merge.sub

############
awk '{print $1,$2-1,"+","CpG",$5,substr($4,5)}' neo/meth_count.csv.bed > meth/neo
awk '$6>0 {print $1,$2,"+","CpG",$5,$6}' meth/neo > merge/neo
for f in meth/*WBC.bed; do echo $f; LC_ALL=C sort -k 1,1 -k 2,2n $f -o $f; done
for f in meth/*WBC.bed; do echo $f; awk '!seen[$1,$2]++' $f > $f.unq; done
dnmtools merge meth/*WBC.unq -o merge/wbc

for f in ${l[@]}; do dnmtools roi hmr/merge merge/$f > hmr/$f.s; done
time paste hmr/*.s | awk '{printf "%s,%.0f,%.0f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n", $1,$2,$3,$5,$11,$17,$23,$29,$35,$41,$47,$53,$60}' > hmr/sim

