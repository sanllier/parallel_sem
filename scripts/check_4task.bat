@echo off
set SIZE=%1
set POW=%2
set N=%3

call rand_matr.exe -f a.mtrx -h %SIZE% -w %SIZE%
call task4.exe -check true -n %POW% -f a.mtrx -res d_.mtrx
call mpiexec -n %N% task4.exe -n %POW% -f a.mtrx -res d.mtrx
call FC /B d.mtrx d_.mtrx