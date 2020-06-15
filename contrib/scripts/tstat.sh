#!/bin/bash
# Limpando arquivos
rm *.md TBMC.dat
# Atribuindo valores iniciais
aini=0.01
afin=0.001
passo=-0.001
a=$aini
while [ $(echo "$a>=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*ensemble.*/ensemble npt $a 0.1 4.9d-5/" TBMC.in
# Rodando TBMC
./TBMC
./TBMC2
t=$(awk '/DSV1/ {print $5}' TBMC.dat)
de=$(awk '/DSV1/ {print $8}' TBMC.dat)
#imprimindo em arquivo dat
echo $a $de $t >> tstat.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
