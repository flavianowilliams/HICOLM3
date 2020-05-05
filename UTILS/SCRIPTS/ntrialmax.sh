#!/bin/bash
# Limpando arquivos
rm ntrialmax.dat ntrialmax-*.out
# Atribuindo valores iniciais 
aini=10000
afin=100000
passo=10000
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*ntrialmax.*/ntrialmax $a/" TBMC.in
# Rodando TBMC
./TBMC
es=$(awk '/simula.:/ {print $3}' TBMC.out)
E=$(awk '/TOTL):/ {print $2}' TBMC.out)
dE=$(awk '/TOTL):/ {print $3}' TBMC.out)
v=$(awk '/F):/ {print $2}' TBMC.out)
dv=$(awk '/F):/ {print $3}' TBMC.out)
cp TBMC.out ntrialmax-$a.out
# Imprimindo em arquivo dat
echo $a $es $v $dv $E $dE >> ntrialmax.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
