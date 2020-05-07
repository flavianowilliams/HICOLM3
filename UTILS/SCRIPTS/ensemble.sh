# Limpando arquivos
rm *.dat *.out
# Atribuindo valores iniciais 
aini=1
afin=0.4
passo=-0.1
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*ensemble.*/ensemble npt hoover $a $a/" HICOLM.in
# Rodando siesta
./HICOLM
./propriedade
tst=$(awk '/MED3/ {print $3}' HICOLM.dat)
# mv TB.out TB-$a.out
# Imprimindo em arquivo dat
echo $a $tst >> ensemble.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
