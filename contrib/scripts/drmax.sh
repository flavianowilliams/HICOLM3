#!/bin/bash
# Limpando arquivos
rm drmax.dat
# Atribuindo valores iniciais 
aini=0.1
afin=0.01
passo=-0.001
a=$aini
while [ $(echo "$a>=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*drmax.*/drmax $a/" TBMC.in
# Rodando TBMC
./TBMC
es=$(awk '/simula.:/ {print $3}' TBMC.out)
E=$(awk '/TOTL):/ {print $2}' TBMC.out)
dE=$(awk '/TOTL):/ {print $3}' TBMC.out)
# Imprimindo em arquivo dat
echo $a $es $E $dE >> drmax.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
