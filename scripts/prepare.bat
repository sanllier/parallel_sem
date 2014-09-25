@echo off
set H=%1
set W=%2

call rand_matr.exe -f a.mtrx -h %H% -w %W%
call rand_matr.exe -f b.mtrx -h %W% -w %H%
call task1.exe -af a.mtrx -bf b.mtrx -resf c_.mtrx