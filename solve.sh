for V in 50

do

    for l in 5 10 15

    do

        for d in 25 50 75

        do

            for c in 25 50 75 

            do

                for i in 1 2 3 4 5

                do

                    #./listcola instances/test/V$V.d$d.c$c.l$l.i$i
                    ./main instances/test/V$V.d$d.c$c.l$l.i$i

                done

            done

        done

    done

done
