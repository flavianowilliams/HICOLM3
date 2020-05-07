#!/bin/bash
# Limpando arquivos
rm densidade.dat
# Atribuindo valores iniciais 
aini=73
afin=373
passo=25
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*text.*/text $a/" TBMC.in
# Rodando TBMC
./TBMC
d=$(awk '/DENSIDADE/ {print $2}' TBMC.out)
# Imprimindo em arquivo dat
echo $a $d >> densidade.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
