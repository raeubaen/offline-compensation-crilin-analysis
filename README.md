Plot da fare in ordine di priorità:
 - A calibrazione cherenlov: sum(cherenkov) vs EnergyTotal, fit lineare e estrazione coeff. calibrazione
 - B risoluzione con DCB (di (crilin + vd) / primary - 1 ) vs energia e singole DCB in slice in energia, con rete e senza rete [ANALYZER.PY], con e senza tagliare le MIP, smearando HCAL al 25%
 - B+, valutare la risoluzione fatta con la DCB e quella fatta fittando DCB e poi la FWHM/2.35 della curva della DCB
 - B++ smearare HCAL a 30, 35, 40, 45, 50%
 - C plot dello sciame 3D (con opzione Box2 in ROOT) a varie energie, profondità dello sciame, MIP/non-MIP, larghezza sciame, etc.
   
   tutti quelli sotto sono da fare tagliando le MIP
 - D RMS cluster vs CRILIN_RECO (somma cherenkov calibrata) / CRILIN_TRUE (Primary - VD energy), in slice in energia, con fit lineare al profilo, sovrapposto al TH2, tagliando MIP
 - D+ calcolo correzione lineare in base all'RMS, risoluzione in slice di energia (fatta con fit DCB e poi FWHM/2.35), e confronto con B
 - E D ma con la media in Z o un'altra variabile sensata in Z (per esempio la prima hit in Z), a slice in energia
 - E+ idem, ma con altre variabili che ti paiono sensate (la prima hit in Z, numero hit sopra una soglia, altre?)
 - F risoluzione con la rete in base a RMS/Z/le altre variabili, e risoluzione in base al layer dove sciama il pione
 - G andamento della correzione (crilin fatto dalla rete / crilin fatto sommando le hit) vs. RMS, Z, le altre variabili che hai scelto
 - G+ andamento della correzione plottandola come il peso (quindi il colore di un TH2) e mettendo RMS e Z in x e y

Dati in:
 https://rgargiul.web.cern.ch/dati_rete_per_vittoria/
