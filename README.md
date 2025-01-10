# TODO

- [x] domande 1/2/3/4
- [ ] domanda 10: codice per GPU
- [x] fare grafico per arithmetic analysis
- [x] implementazione di Openmp
- [x] domanda 5: scalasca su differenti topologie e fare grafico per vedere ottimizzazione
    - [x] strong scaling
        - [x] hybrid
        - [x] mpi
        - [x] OpenMp -> OK
    - [x] weak scaling
        - [x] hybrid
        - [x] mpi
        - [x] OpenMp
- [x] domanda 6: scalasca su collapse e non collapse
    - [x] collapse -> stessi risultati di openMP
    - [x] no collapse -> da ricompilare con ottimizzazione
- [x] domanda 7
- [ ] Scalasca per spiegare risultati del weak e strong scaling
- [ ] ultima domanda batimetria
    - [x] trovare dati
    - [x] descrizione del dataset (nx, ny, dx, dy)
    - [x] conversione nel formato che ci serve a noi
    - [x] simulazione
    - [x] aggiungere altri tipi di source
    - [ ] scaricare risultati e fare video delle 3 source


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

## 03/01

Ho aggiuto i risultati per lo scaling con no collapse.

## 04/01

I parametri del file di batimetria semplice che abbiamo sempre usato sono i seguenti:

| name     | value |
|----------|-------|
| nx       | 100   |
| ny       | 100   |
| dx       | 40.0  |
| dy       | 40.0  |

Dopo averci pensato un pò penso che la meglio cosa sia prendere un dataset della batimetria delle terre non emerse e mettere come condizione di riflesso totale i contorni delle terre emerse. Per questo motivo ho usato il sito [emodnet.ec.europa.eu](https://emodnet.ec.europa.eu/geoviewer/#) che permette di recuperare i dati della batimetria. Continuerei ad usare un dataset con la stessa shape (non numericamente ma almeno che sia un quadrato) giusto per non avere problemi nella parallelizzazione. Si può usare un pò più di punti. 

Ho recuperato quindi i dati di batimetria vicino all'isola di Pianosa. Ho ripulito il dataset e creato per il momento un file csv contenente i dati della batimetria. La shape dell'array è $150\cdot 150$.

Tramite questo [calcolatore](https://opendem.info/arc2meters.html) si può calcolare i due parametri _dx_ e _dy_ prendendo una latitudine di riferimento per l'intero dataset (assumendo quindi una griglia regolare anche se nella realtà non lo è). 

Per il momento ho fatto uno zip in cui ho messo i dati raw ottenuti dal sito e poi i miei dati processati nell'array $150\cdot 150$.

## 05/01

Ho calcolato alcuni parametri necessari per la descrizione del dataset della batimetria:

- A partire dal dataset _raw_ ho estratto un quadrato di $150$ punti in entrambe le direzioni ($150\cdot 150 = 22500$ punti);
- Le estremità del quadrato sono le seguenti:

    | index         | latitude       | longitude |
    |----------------|----------------|-----------|
    | $[0, 0]$       | 42.6651        | 9.9724    |
    | $[0, 149]$     | 42.6651        | 10.1286   |
    | $[149, 0]$     | 42.5089        | 9.9724    |
    | $[149, 149]$   | 42.5089        | 10.1286   |

- A partire da queste coordinate in gradi si può effettuare la conversione per ottenere il valore di _dx_ e _dy_ in metri. Per prima cosa si calcola la differenza in gradi fra le coordinate di inizio e di fine per la longitude e la latitude. Da questo si calcola _dx_ e _dy_ (in gradi) prendendo la differenza e dividendola per _n_ numero di punti lungo l'asse corrispondente (nel nostro caso 150 per tutte e due le direzioni):

    |                                    | latitude       | longitude |
    |------------------------------------|----------------|-----------|
    | absolute difference                | 0.1562         | 0.1563    |
    | discretization step (in degree)    | 0.0010417      | 0.0010417 |

    Si può assumere quindi che lo step sia lo stesso lungo entrambi gli assi e costante. Quindi si può calcolare la lunghezza dello step in metri usando il [calcolatore](https://opendem.info/arc2meters.html) considerando che $0.0010417° = 3.75012$ _seconds_ e che la latitude media (si assume che la griglia sia costante cosa che in realtà non sarebbe data la differenza di latitude ma che in questo caso è minima) è $42.587$. Otteniamo un valore di $85.22374$ _m_.

    |                                    | latitude       | longitude |
    |------------------------------------|----------------|-----------|
    | absolute difference between max and min (in degree)                | $0.1562   $      | $0.1563   $ |
    | discretization step (in degree)    | $0.0010417$      | $0.0010417$ |
    | discretization step (in arc seconds)                        | $3.75012  $      | $3.75012  $ |
    | discretization step based on average latitude (in meters)| $85.22\approx 85$ | $85.22\approx 85$ |

Riassumendo abbiamo le seguenti caratteristiche per il dataset:

| name     | value                      |
|----------|----------------------------|
| nx       | 150                        |
| ny       | 150                        |
| dx       | $85.22374 \approx 85$      |
| dy       | $85.22374 \approx 85$      |


## 06/01

Ho aggiustato l'eseguibile _pianosa.c_ per fare in modo che i contorni dell'isola fossero riflettenti. Ho runnato e tutto va bene ma nella soluzione ci sono alcuni valori abbastanza alti $\approx 300$ che sembrano strani. La simulazione sembra realistica.