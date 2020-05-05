# Limpando arquivos
rm *.dat *.out
# Atribuindo valores iniciais 
aini=1
afin=12
passo=1
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*monkhorst.*/monkhorst $a $a 1/" TB.in
# Rodando siesta
./TB
E=$(awk '/Energia/ {print $3}' TB.out)
mv TB.out TB-$a.out
# Imprimindo em arquivo dat
echo $a $E >> teste.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
