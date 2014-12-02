@echo off
set N=%1

for /l %%x in (100, 100, 2000) do (
call rand_matr.exe -f a.mtrx -h %x% -w %x%
call rand_matr.exe -f b.mtrx -h %x% -w %x%

call mpiexec -n %N% task5.exe -mata a.mtrx -matb b.mtrx -res c_1.mtrx
call mpiexec -n %N% task5_2.exe -mata a.mtrx -matb b.mtrx -res c_2.mtrx
call FC /B c_1.mtrx c_2.mtrx
)