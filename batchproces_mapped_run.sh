#parallel -j 8 --ctag ./process_mapped_run.sh {1} ::: {01..04}
#hadd -j
instr=""
declare -i i=0
while [[ $i -lt 10 ]]; do
echo $i
i=$i+1
instr=" "+$i
done
echo $instr
