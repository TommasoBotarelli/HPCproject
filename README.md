# TODO

- [x] domande 1/2/3/4
- [ ] domanda 10: codice per GPU
- [x] fare grafico per arithmetic analysis
- [x] implementazione di Openmp
- [ ] domanda 5: scalasca su differenti topologie e fare grafico per vedere ottimizzazione
    - [x] strong scaling
        - [x] hybrid
        - [x] mpi
        - [x] OpenMp -> OK
    - [x] weak scaling
        - [x] hybrid
        - [x] mpi
        - [x] OpenMp
- [ ] domanda 6: scalasca su collapse e non collapse
    - [x] collapse -> stessi risultati di openMP
    - [ ] no collapse -> da ricompilare con ottimizzazione
- [x] domanda 7
- [ ] domande 8 e 9: fare scalasca su tutte le configurazioni e riportare su Excel
- [x] per togliersi tutti i dubbi vedere che succede ai risultati lanciati con un solo processo mpi. Se sono corretti allora effetivamente i tempi registrati sono ok. -> Update: sono corretti cioè anche con un solo processo si genera il file giusto (grandezza 4000x4000)


# UPDATEs

## 20/12

Riguardando i risultati che sono anche nella cartella results non mi tornava alcune cose:
- Rispetto ai risultati sequenziali MPI con **un solo processo** (che nella pratica dovrebbe essere paragonabile al precedente) impiegava ben oltre i circa 20 secondi;
- OpenMP aveva degli andamenti ancora peggiori e lo speedup era osceno;

Ho quindi notato che:
- MPI stampava molte cose a console e quindi ho ripulito la stampa;
- Mentre per il sequenziale usavamo l'opzione _-O3_ in fase di compilazione (ottimizzazione di terzo livello) questo non veniva fatto negli altri eseguibili. 

Quindi ho ricompilato tutti i file e rilanciato tutti gli esempi. Ho anche modificato lo script bash per la generazione della tabella riprendendo anche l'effettivo tempo di esecuzione che abbiamo sempre preso per determinare la velocità del codice (ora la tabella è completa infatti).

I risultati sono riportati nei file _.out_ della cartella _results_. Ora però c'è una nuova cosa un pò strana... MPI con un singolo processo dovrebbe avere il tempo di esecuzione di un sequenziale ma sembrerebbe avere un tempo inferiore.

Non di molto ma sembra un pò strano che si guadagni quasi 3 secondi. C'é da considerare però che i tempi del sequenziale sono stati mediati su 50 esecuzioni e questa invece è una singola esecuzione. Inoltre è probabile che le impostazioni con la quale sono stati registrati i tempi del sequenziale non rispecchino le stesse impostazioni di esclusività del nodo con la quale sono stati runnati i test di strong scaling.

Ricontrollando sembrerebbe che le impostazioni del job con la quale sono stati lanciati siano le stesse. Tuttavia mi sono accorto che effettivamente i sequenziali potrebbero avere qualche tipo di overhead per il fatto che si registrano i tempi. Anche se continua a sembrare un pò strano...

Ho girato senza esclusività MPI e sembra che con un solo nodo comunque faccia tempi molto migliori. Questo continua a non tornarmi molto...

In ogni caso ho aggiunto i risultati con ottimizzazione dello strong scaling nella cartella results.

## 21/12

Fatti i test per weak scaling. I risultati sembrano avere senso. Lo scaling è stato effettuato nel seguente modo:

| n ranks | dx    | dy    | total points             |
|---------|-------|-------|--------------------------|
| 1       | 5     | 5     | $800\cdot800=640000$     |
| 2       | 2.5   | 5     | $1600\cdot800=1280000$   |
| 4       | 2.5   | 2.5   | $1600\cdot1600=2560000$  |
| 8       | 1.25  | 2.5   | $3200\cdot1600=5120000$  |
| 16      | 1.25  | 1.25  | $3200\cdot3200=10240000$ |
| 32      | 0.625 | 1.25  | $6400\cdot3200=20480000$ |
| 64      | 0.625 | 0.625 | $6400\cdot6400=40960000$ |

Non è stato controllato che effettivamente i risultati siano stabili (dal momento che sono stati mantenuti lo stesso numero di punti nello spazio di sempre e questo potebbe comportare per _dx, dy_ piccoli instabilità).

Per i risultati del tempo di computazione ho messo tutto nella cartella results.

Ho fatto anche la weak e lo strong scaling per il caso di openmp senza il collapse. I risultati li devo ancora mettere.