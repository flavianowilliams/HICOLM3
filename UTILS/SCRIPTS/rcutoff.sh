#!/bin/bash
# Limpando arquivos
rm rcutoff.dat
# Atribuindo valores iniciais 
aini=4
afin=10
passo=0.5
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*rcutoff.*/rcutoff $a 0.1/" TBMC.in
# Rodando TBMC
./TBMC
t=$(awk '/Tempo estimado/ {print $4}' TBMC.out)
c=$(awk '/energia de vdw/ {print $6}' TBMC.out)
# Imprimindo em arquivo dat
echo $a $t $c >> rcutoff.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
