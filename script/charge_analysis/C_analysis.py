import tarfile
import re
import numpy as np
import matplotlib.pyplot as plt

NOME_FILE_ARCHIVIO = "CASF-2016.tar.gz" 

cariche_totali = []
file_analizzati = 0

print(f"Apertura dell'archivio {NOME_FILE_ARCHIVIO} (senza estrarlo sul disco)...")
try:
    with tarfile.open(NOME_FILE_ARCHIVIO, "r:gz") as tar:
        for member in tar.getmembers():
            # Vogliamo solo i ligandi in formato mol2
            if member.name.endswith(".mol2") and "ligand" in member.name:
                
                f = tar.extractfile(member)
                if f is not None:
                    # Leggiamo il testo decodificandolo
                    righe = f.read().decode('utf-8').splitlines()
                    
                    dentro_blocco_atomi = False
                    
                    for riga in righe:
                        # Capiamo quando iniziano e finiscono gli atomi
                        if riga.startswith("@<TRIPOS>ATOM"):
                            dentro_blocco_atomi = True
                            continue
                        elif riga.startswith("@<TRIPOS>BOND"):
                            dentro_blocco_atomi = False
                            break
                        
                        # Se siamo dentro la lista degli atomi, estraiamo la carica
                        if dentro_blocco_atomi:
                            colonne = riga.split()
                            if len(colonne) >= 9: # Assicuriamoci che la riga sia completa
                                try:
                                    # L'ultima colonna è la carica
                                    carica = float(colonne[-1])
                                    cariche_totali.append(carica)
                                except ValueError:
                                    pass
                
                file_analizzati += 1
                if file_analizzati % 500 == 0:
                    print(f"Analizzati {file_analizzati} ligandi...")

    print("-" * 30)
    print("ANALISI COMPLETATA!")
    print(f"Totale ligandi analizzati: {file_analizzati}")
    print(f"Totale atomi (cariche) estratti: {len(cariche_totali)}")
    
    if len(cariche_totali) > 0:
        print(f"Carica MINIMA trovata: {min(cariche_totali):.4f}")
        print(f"Carica MASSIMA trovata: {max(cariche_totali):.4f}")

except FileNotFoundError:
    print("ERRORE: Il file compresso non è stato trovato nella cartella.")

if len(cariche_totali) > 0:
    cariche_array = np.array(cariche_totali)
    
    # 1. Calcoliamo i veri confini (Ignoriamo gli outlier estremi)
    p1 = np.percentile(cariche_array, 1)   # Il confine dell'1% più basso
    p99 = np.percentile(cariche_array, 99) # Il confine del 99% più alto
    
    print(f"\n--- ANALISI STATISTICA ---")
    print(f"Valor Medio: {np.mean(cariche_array):.4f}")
    print(f"Deviazione Standard: {np.std(cariche_array):.4f}")
    print(f"Il 98% degli atomi ha una carica compresa tra {p1:.4f} e {p99:.4f}")
    
    plt.figure(figsize=(10, 6))
    
    # Prova a cambiare 'bins=50' con 20 o 100 per vedere come cambia
    plt.hist(cariche_array, bins=100, color='skyblue', edgecolor='black')
    
    plt.axvline(p1, color='red', linestyle='dashed', linewidth=2, label=f'1° Percentile ({p1:.2f})')
    plt.axvline(p99, color='red', linestyle='dashed', linewidth=2, label=f'99° Percentile ({p99:.2f})')
    
    plt.title("Distribuzione delle Cariche Parziali (CASF-2016)")
    plt.xlabel("Carica Elettrica")
    plt.ylabel("Frequenza (Numero di Atomi)")
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    
    plt.show()