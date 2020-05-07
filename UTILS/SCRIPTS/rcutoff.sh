#!/bin/bash
# Limpando arquivos
rm rcutoff.dat
# Atribuindo valores iniciais 
aini=11
afin=8
passo=-0.5
a=$aini
while [ $(echo "$a>=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*rcutoff.*/rcutoff $a 0.1/" HICOLM.in
# Rodando TBMC
./HICOLM
./properties
et=$(awk '/MED2/ {print $8}' HICOLM.dat)
at=$(awk '/Total of atoms/ {print $4}' HICOLM.out)
# Imprimindo em arquivo dat
echo $a $et $at >> rcutoff.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
