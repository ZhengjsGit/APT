#compare 3 symplectic algorithms using 160 particles
#times 3e8 6e7 3e7
#dT	   0.01 0.05 0.1
#Init E 0.1 which means Regular
#./calc.sh 8 24 0.1 3e7 100 Boris_Regular
#./calc.sh 8 24 0.05 6e7 200 Boris_Regular
#./calc.sh 8 24 0.01 3e8 1000 Boris_Regular

./calc_para.sh 160 160 10 24 0.1 3e7 100 SymE_Regular_para
./calc_para.sh 160 160 10 24 0.05 6e7 200 SymE_Regular_para
./calc_para.sh 160 160 10 24 0.01 3e8 1000 SymE_Regular_para

./calc_para.sh 160 160 11 24 0.1 3e7 100 imE_Regular_para
./calc_para.sh 160 160 11 24 0.05 6e7 200 imE_Regular_para
./calc_para.sh 160 160 11 24 0.01 3e8 1000 imE_Regular_para

./calc_para.sh 160 160 12 24 0.1 3e7 100 SV_Regular_para
./calc_para.sh 160 160 12 24 0.05 6e7 200 SV_Regular_para
./calc_para.sh 160 160 12 24 0.01 3e8 1000 SV_Regular_para
