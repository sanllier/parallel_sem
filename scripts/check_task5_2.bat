@echo off
set HA=%1
set WA=%2
set WB=%3
set N=%4

call rand_matr.exe -f a.mtrx -h %HA% -w %WA%
call rand_matr.exe -f b.mtrx -h %WA% -w %WB%

call task1 -af a.mtrx -bf b.mtrx -resf c_.mtrx
call mpiexec -n %N% task5_2.exe -mata a.mtrx -matb b.mtrx -res c.mtrx
call FC /B c.mtrx c_.mtrx