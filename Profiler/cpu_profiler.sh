#!/bin/bash
#CPU profiling with perf


echo -----Running profiler script...-----
FileName=$1
FileFormat=".log"
LogFileName="$FileName$FileFormat"
echo "Creating ""$LogFileName"" from ""$1"

rm -f $LogFileName

#basic data
perf stat -o $FileName"0"$FileFormat -e cycles,instructions ./$1
perf stat -o $FileName"1"$FileFormat -e r530110,r538010,r532010,r534010,r531111,r530111,r530211 ./$1
perf stat -o $FileName"2"$FileFormat -e cache-misses,cache-references,mem-loads,mem-stores ./$1
perf stat -o $FileName"3"$FileFormat -e L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores,L1-dcache-store-misses,L1-dcache-prefetches,L1-dcache-prefetch-misses ./$1
perf stat -o $FileName"4"$FileFormat -e L1-icache-loads,L1-icache-load-misses,L1-icache-prefetches,L1-icache-prefetch-misses ./$1
perf stat -o $FileName"5"$FileFormat -e LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,LLC-prefetches,LLC-prefetch-misses ./$1
perf stat -o $FileName"6"$FileFormat -e dTLB-loads,dTLB-load-misses,dTLB-stores,dTLB-store-misses,dTLB-prefetches,dTLB-prefetch-misses ./$1
perf stat -o $FileName"7"$FileFormat -e iTLB-loads,iTLB-load-misses,branches,branch-misses,branch-loads,branch-load-misses ./$1
perf stat -o $FileName"8"$FileFormat -e r53a3a3,r15301a3,r25302a3,r85308a3,r45304a3,r55305a3,rc530ca3,r65306a3 ./$1
perf stat -o $FileName"9"$FileFormat -e r531414,r530114,r531014 ./$1
perf stat -o $FileName"10"$FileFormat -e r53b1b1,r5302b1,r5301b1,r5318b1,r5310b1,r5320b1,r5330b1,r5340b1 ./$1
perf stat -o $FileName"11"$FileFormat -e r53a1a1,r5301a1,r5302a1,r530ca1,r5330a1,r5340a1,r5380a1 ./$1
perf stat -o $FileName"12"$FileFormat -e r53e10e,r53010e,r1f31a0e,r1d3180e,r53100e,r53200e,r53400e ./$1
perf stat -o $FileName"13"$FileFormat -e r538080,r530280,r530480,r530180 ./$1
perf stat -o $FileName"14"$FileFormat -e r538787,r530187,r530487 ./$1
perf stat -o $FileName"15"$FileFormat -e r53a2a2,r5301a2,r5304a2,r5308a2,r5310a2 ./$1


cat $FileName"0"$FileFormat $FileName"1"$FileFormat $FileName"2"$FileFormat $FileName"3"$FileFormat $FileName"4"$FileFormat $FileName"5"$FileFormat $FileName"6"$FileFormat $FileName"7"$FileFormat $FileName"8"$FileFormat $FileName"9"$FileFormat $FileName"10"$FileFormat $FileName"11"$FileFormat $FileName"12"$FileFormat $FileName"13"$FileFormat $FileName"14"$FileFormat $FileName"15"$FileFormat > $LogFileName

for i in {0..15}
do
	rm $FileName$i$FileFormat
done

sed -i '/^$/d' ./$LogFileName
sed -i '/started/d' $LogFileName 
sed -i '/Performance/d' $LogFileName

sed -i 's/cycle/Number of clock cycles/g' $LogFileName
sed -i 's/instructions/Number of instruction/g' $LogFileName
sed -i 's/r531010/Counts number of floating point events/g' $LogFileName
sed -i 's/r530110/Number of X87 uops executed/g' $LogFileName
sed -i 's/r531010/Number of SSE or AVX-128 double/g' $LogFileName
sed -i 's/r538010/Number of SSE or AVX-128 double precision FP scalar uops executed/g' $LogFileName
sed -i 's/r532010/Number of SSE or AVX-128 single precision FP scalar uops executed/g' $LogFileName
sed -i 's/r534010/Number of SSE or AVX-128 single precision FP packed uops executed/g' $LogFileName
sed -i 's/r531111/Counts 256-bit packed floating point/g' $LogFileName
sed -i 's/r530111/Counts 256-bit packed single-precision/g' $LogFileName
sed -i 's/r530211/Counts 256-bit packed double-precision/g' $LogFileName
sed -i 's/r53a3a3/Stalled cycles/g' $LogFileName
sed -i 's/cache-misses/Cache misses/g' $LogFileName
sed -i 's/cache-references/Cache references/g' $LogFileName
sed -i 's/r15301a3/CYCLES_L2_PENDING/g' $LogFileName
sed -i 's/r25302a3/CYCLES_LDM_PENDING/g' $LogFileName
sed -i 's/r85308a3/CYCLES_L1D_PENDING/g' $LogFileName
sed -i 's/r45304a3/CYCLES_NO_EXECUTE/g' $LogFileName
sed -i 's/r55305a3/STALLS_L2_PENDING /g' $LogFileName
sed -i 's/rc530ca3/STALLS_L1D_PENDING/g' $LogFileName
sed -i 's/r65306a3/STALLS_LDM_PENDING/g' $LogFileName
sed -i 's/r531414/Counts arithmetic multiply operations /g' $LogFileName
sed -i 's/r530114/Cycles that the divider is active, includes integer and floating point/g' $LogFileName
sed -i 's/r531014/Number of cycles the divider is activated, includes integer and floating point/g' $LogFileName
sed -i 's/r53b1b1/Uops executed/g' $LogFileName
sed -i 's/r5302b1/Counts total number of uops executed from any thread per cycle/g' $LogFileName
sed -i 's/r5301b1/Counts total number of uops executed per thread each cycle/g' $LogFileName
sed -i 's/r5318b1/Number of cycles with no uops executed/g' $LogFileName
sed -i 's/r5310b1/Cycles where at least 1 uop was executed per thread/g' $LogFileName
sed -i 's/r5320b1/Cycles where at least 2 uops were executed per thread/g' $LogFileName
sed -i 's/r5330b1/Cycles where at least 3 uops were executed per thread/g' $LogFileName
sed -i 's/r5340b1/Cycles where at least 4 uops were executed per thread/g' $LogFileName
sed -i 's/r53a1a1/Uops dispatch to specific ports/g' $LogFileName
sed -i 's/r5301a1/Cycles in which a uop is dispatched on port 0/g' $LogFileName
sed -i 's/r5302a1/Cycles in which a uop is dispatched on port 1/g' $LogFileName
sed -i 's/r530ca1/Cycles in which a uop is dispatched on port 2/g' $LogFileName
sed -i 's/r5330a1/Cycles in which a uop is dispatched on port 3/g' $LogFileName
sed -i 's/r5340a1/Cycles in which a uop is dispatched on port 4/g' $LogFileName
sed -i 's/r5380a1/Cycles in which a uop is dispatched on port 5/g' $LogFileName
sed -i 's/r53e10e/Uops issued/g' $LogFileName
sed -i 's/r53010e/Number of uops issued by the RAT to the Reservation Station (RS)/g' $LogFileName
sed -i 's/r1f31a0e/Alias to ANY:c=1:i=1:t=1/g' $LogFileName
sed -i 's/r1d3180e/Alias to ANY:c=1:i=1/g' $LogFileName
sed -i 's/r53100e/Number of flags-merge uops allocated. Such uops adds delay/g' $LogFileName
sed -i 's/r53200e/Number of slow LEA or similar uops allocated/g' $LogFileName
sed -i 's/r53400e/Number of multiply packed scalar single precision uops allocated/g' $LogFileName
sed -i 's/r538080/Ilnstruction Cache accesses/g' $LogFileName
sed -i 's/r530280/Number of Instruction Cache, Streaming Buffer and Victim Cache Misses. Includes UC accesses/g' $LogFileName
sed -i 's/r530480/Number of cycles wher a code-fetch stalled due to L1 instruction cache miss or iTLB miss/g' $LogFileName
sed -i 's/r530180/Number of Instruction Cache,Buffer/g' $LogFileName
sed -i 's/r538787/Instruction Length Decoder stalls/g' $LogFileName
sed -i 's/r530187/Stall caused by changing prefix length of the instruction/g' $LogFileName
sed -i 's/r530487/Stall cycles due to IQ full/g' $LogFileName
sed -i 's/r53a2a2/Resource related stall cycles/g' $LogFileName
sed -i 's/r5301a2/Cycles stalled due to Resource Related reason/g' $LogFileName
sed -i 's/r5304a2/Cycles stalled due to no eligible RS entry available/g' $LogFileName
sed -i 's/r5308a2/Cycles stalled due to no store buffers available/g' $LogFileName
sed -i 's/r5310a2/Cycles stalled due to re-order buffer full/g' $LogFileName

echo -----Script finished...-----