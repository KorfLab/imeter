import os

offset = []
slope = []
rate = []
for i in range(1001):
    if i % 100 == 0 and i != 0:
        offset.append(i)
        slope.append(1/i)
        rate.append('1/'+str(i))
#print(offset, slope, rate)

for off in offset:
    for i in range(len(slope)):
        trainer = f'python3 geometric_trainer.py at_ime_master.txt.gz {off} {slope[i]} > mod'
        os.system(trainer)
        decoder = f'python3 decoder.py mod db_IME_Rose_WT_introns.fa > introns'
        os.system(decoder)
        r2 = f'python3 r2_analysis.py introns {off} {rate[i]}'
        os.system(r2)
