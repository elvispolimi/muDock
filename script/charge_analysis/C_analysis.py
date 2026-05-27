import tarfile
import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

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

NUM_BINS = 12 
    
#     print(f"\n--- GENERAZIONE THRESHOLD NON-LINEARI ({NUM_BINS} BIN) ---")
    
#     cariche_pulite = cariche_array[(cariche_array >= p1) & (cariche_array <= p99)]
#     percentuali_taglio = np.linspace(0, 100, NUM_BINS + 1)[1:-1]
    
#    # 3. Troviamo il valore esatto della carica a quelle percentuali
#     thresholds_grezzi = np.percentile(cariche_pulite, percentuali_taglio)
    
#     # Eliminiamo i duplicati, Se 5 bin cadono su 0.0, ne teniamo solo uno.
#     thresholds_unici = np.unique(thresholds_grezzi)
#     thresholds_str = ", ".join([f"{t:.5f}" for t in thresholds_unici])
#     print("\nCopia e incolla questa riga nel costruttore della tua classe C++:")
#     print("std::vector<fp_type> thresholds = {")
#     print(f"    {thresholds_str}")
#     print("};")
print(f"\n--- GENERAZIONE THRESHOLD K-MEANS ({NUM_BINS} BIN ESATTI) ---")

cariche_pulite = cariche_array[(cariche_array >= p1) & (cariche_array <= p99)]
cariche_reshaped = cariche_pulite.reshape(-1, 1)
kmeans = KMeans(n_clusters=NUM_BINS, random_state=42, n_init="auto")
kmeans.fit(cariche_reshaped)

centri_ordinati = np.sort(kmeans.cluster_centers_.flatten())

thresholds_esatti = (centri_ordinati[:-1] + centri_ordinati[1:]) / 2.0

thresholds_str = ", ".join([f"{t:.5f}" for t in thresholds_esatti])

print("\nCopia e incolla questo array nel tuo autodock_quant_protein.hpp:")
print(f"// Generato per avere ESATTAMENTE {NUM_BINS} mappe quantizzate")
print("static const std::vector<fp_type> thresh = {")
print(f"    {thresholds_str}")
print("};")