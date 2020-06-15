#!/bin/bash
# Limpando arquivos
rm densidade-press.dat
# Atribuindo valores iniciais 
aini=100
afin=1000
passo=100
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*preext.*/preext $a/" TBMD.in
# Rodando TBMC
./TBMD
./tbmc
d=$(awk '/MED2/ {print $3}' TBMD.dat)
ds=$(awk '/DSV2/ {print $3}' TBMD.dat)
# Imprimindo em arquivo dat
echo $a $d $ds >> densidade-press.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
