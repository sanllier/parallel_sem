@echo off
set N=%1

call mpiexec -n %N% task3.exe -mata a.mtrx -matb b.mtrx -res c.mtrx
call FC /B c.mtrx c_.mtrx