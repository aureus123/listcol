for V in 50

do                                   

    for d in 25 50 75

    do

        ./gengraph V$V.d$d.graph $V $d

        for c in 25 50 75

        do

            for l in 5 10 15

            do

                for i in 1 2 3 4 5
                do

                    cp V$V.d$d.graph V$V.d$d.c$c.l$l.i$i.graph
                    ./genrandominst2 V$V.d$d.c$c.l$l.i$i $c $l 1 1

                done

            done

        done

    done

done
