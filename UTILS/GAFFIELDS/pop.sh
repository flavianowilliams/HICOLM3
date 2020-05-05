# Limpando arquivos
rm *pop.dat ag.out
# Atribuindo valores iniciais
aini=1000
afin=1500
passo=50
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
sed -i -e "s/.*#AQUI.*/$a #AQUI/" ag.in
# Rodando gaffields
./ag > ag.out
E=$(awk '/Erro estimado:/ {print $3}' ag.out)
#mv TB.out TB-$a.out
# Imprimindo em arquivo dat
echo $a $E >> pop.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
