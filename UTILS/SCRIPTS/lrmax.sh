#!/bin/bash
# Limpando arquivos
rm lrmax.dat lrmax-*.out
# Atribuindo valores iniciais 
aini=0.001
afin=0.0001
passo=-0.00005
a=$aini
while [ $(echo "$a>=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*lrmax.*/lrmax $a 2/" TBMC.in
# Rodando TBMC
./TBMC
es=$(awk '/simula.:/ {print $3}' TBMC.out)
E=$(awk '/TOTL):/ {print $2}' TBMC.out)
dE=$(awk '/TOTL):/ {print $3}' TBMC.out)
v=$(awk '/F):/ {print $2}' TBMC.out)
dv=$(awk '/F):/ {print $3}' TBMC.out)
cp TBMC.out lrmax-$a.out
# Imprimindo em arquivo dat
echo $a $es $v $dv $E $dE >> lrmax.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
