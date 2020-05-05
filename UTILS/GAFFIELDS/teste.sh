# Limpando arquivos
rm teste.dat ag.out
# Atribuindo valores iniciais
aini=1
afin=10
passo=1
a=$aini
while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
# Alterando input
#sed -i -e "s/.*#AQUI.*/$a 0.02 #AQUI/" ag.in
# Rodando gaffields
./ag > ag.out
E1=$(awk '/Individuo selecionado:/ {print $3}' ag.out)
E2=$(awk '/Individuo selecionado:/ {print $4}' ag.out)
E3=$(awk '/Individuo selecionado:/ {print $5}' ag.out)
E4=$(awk '/Individuo selecionado:/ {print $6}' ag.out)
E5=$(awk '/Individuo selecionado:/ {print $7}' ag.out)
#mv TB.out TB-$a.out
# Imprimindo em arquivo dat
echo $a $E1 $E2 $E3 $E4 $E5 >> teste.dat
# Atribuindo valores
a=$(echo "$a+$passo"|bc)
done
