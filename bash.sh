for f in beta/blood/*.beta; do  fn="${f%%.*}";  echo ${fn:5:-10}.bed;  /oak/stanford/scg/lab_mpsnyder/moqri/soft/wgbs_tools/src/python/wgbs_tools.py beta2bed --genome hg38 $f > bed/${fn:11:-10}.bed; done

for f in bed/*.bed; do  fn="${f%.*}";  echo ${fn:4};  bedtools subtract -a $f -b mask > mask/${fn:4}; done
for f in mask/*; do  echo $f;  awk '{print $1,$2,"+","CpG",$4/$5,$5}' $f > meth/${f:5}; done

dnmtools merge meth/*_Blood-T-CD3 -o merge/cd3
dnmtools merge meth/*_Blood-T-CD4 -o merge/cd4
dnmtools merge meth/*_Blood-T-CD8 -o merge/cd8
dnmtools merge meth/*_Blood-NK -o merge/nk
dnmtools merge meth/*_Blood-Monocytes -o merge/mono
dnmtools merge meth/*_Blood-Granulocytes -o merge/gran
dnmtools merge meth/*_Blood-B -o merge/b
dnmtools merge meth/*_Blood-B-Mem -o merge/bmem

dnmtools merge merge/cd3 merge/cd4 merge/cd8 merge/nk -o merge/t
dnmtools merge merge/b merge/t -o merge/lym
dnmtools merge merge/mono merge/gran -o merge/mye

dnmtools hmr-rep merge/lym -o hmr/lym
dnmtools hmr-rep merge/mye -o hmr/mye

dnmtools hmr-rep merge/gram -o hmr/gram
dnmtools hmr-rep merge/t -o hmr/t
dnmtools hmr-rep merge/b -o hmr/b
dnmtools hmr-rep merge/nk -o hmr/nk

bedtools subtract -a lym -b mye -A > lym-mye
bedtools subtract -b lym -a mye -A > mye-lym
awk '{print $1,$2,$3,$3-$2, $4,$5}' hmr/lym-mye | sort -k 4 -rn |  head | sort -k 1,1 -k 2,2n > top/lym
awk '{print $1,$2,$3,$3-$2, $4,$5}' hmr/mye-lym | sort -k 4 -rn |  head | sort -k 1,1 -k 2,2n > top/mye

LC_ALL=C sort -k 1,1 -k 2,2n merge/lym -o merge/lym.sort
LC_ALL=C sort -k 1,1 -k 2,2n merge/mye -o merge/mye.sort

dnmtools roi top/lym merge/lym.sort | awk {'print $1,$2,$3,$5'} > top/lym.lym
dnmtools roi top/lym merge/mye.sort | awk {'print $1,$2,$3,$5'} > top/lym.mye
dnmtools roi top/mye merge/mye.sort | awk {'print $1,$2,$3,$5'} > top/mye.mye
dnmtools roi top/mye merge/lym.sort | awk {'print $1,$2,$3,$5'} > top/mye.lym

awk '{if (substr($0,4,1) != "X" && substr($0,4,1) != "Y" && substr($0,4,1) != "M") print }' h

