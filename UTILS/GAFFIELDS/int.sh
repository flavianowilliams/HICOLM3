# Limpando arquivos
rm int.dat ag.out
# Atribuindo valores iniciais
aini=0.1
afin=2.0
passo=0.25
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*#AQUI.*/$a #AQUI/" ag.in
# Rodando gaffields
./ag > ag.out
E=$(awk '/Erro estimado:/ {print $3}' ag.out)
#mv TB.out TB-$a.out
# Imprimindo em arquivo dat
echo $a $E >> int.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
