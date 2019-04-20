which bwa &> checkPath.txt
which samtools &>> checkPath.txt
which bedtools &>> checkPath.txt
which R &>> checkPath.txt
i=0
check="T"
while read LINE
do
    let i=i+1
    if [ ! -f "$LINE" ]
    then
        if [ $i == 1 ]
        then
            echo "Error: bwa is not in system PATH"
            check="F"
        elif [ $i == 2 ]
        then
            echo "Error: samtools is not in system PATH"
            check="F"
        elif [ $i == 3 ]
        then
            echo "Error: bedtools is not in system PATH"
            check="F"
        elif [ $i == 4 ]
        then
            echo "Error: R is not in system PATH"
            check="F"
        fi
    fi
done < checkPath.txt
rm checkPath.txt
if [ $check = "T" ]
then
    java -cp program/ChIA_PET.jar process.Main parameter.txt
fi
