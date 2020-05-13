# apagando arquivos
rm optimize.dat
#
#aini=1
#afin=50
#passo=1
#a=$aini
#while [ $(echo "$a<=$afin"|bc) -eq 1 ];do
#
i=$(awk '/SD/ {print $2}' HICOLM.out)
E=$(awk '/SD/ {print $4}' HICOLM.out)
# Imprimindo em arquivo dat
echo $i $E >> optimize.dat
#
#a=$(echo "$a+$passo"|bc)
#done

