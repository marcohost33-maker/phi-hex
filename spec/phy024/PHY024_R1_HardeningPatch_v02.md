# PHY024-R1 Hardening Patch v0.2

Datum: 2026-06-01  
Status: angewandte Korrektur nach Audit der letzten Reaktion  
Geltung: ergänzt PHY024-R0 Spec/Runbook, ersetzt aber keine Main-Run-Ergebnisse.

## 1. Korrigierte Punkte

### 1.1 Smoke-Reports existieren, aber waren zu klein
Die vorhandenen R0-Smoke-Reports nutzten nur seeds=4 und radius=2. Sie sind geeignet als
Funktions- und Sanity-Check, nicht als Evidenz für einen Ordnungsmechanismus.

### 1.2 helicity_proxy ist kein echter Helicity Modulus
Im R0-Skript bezeichnet `helicity_analog` nur den Kantenmittelwert
`mean cos(theta_i - theta_j)`. Das ist ein lokaler Steifigkeits-Proxy, nicht die
thermodynamische Helicity-Modulus-Groesse aus der XY-Literatur.

Konsequenz: `helicity_proxy` darf nicht als PASS-Kriterium verwendet werden.
Es bleibt eine Hilfsmetrik.

### 1.3 Binder-like ist nur Hilfsdiagnostik
Der R0-Binder-Indikator ist auf kleinen Graphen stark finite-size-dominiert und darf
nicht als robustes Crossing-Kriterium interpretiert werden.

### 1.4 XY-Kopplung zeigt Synchronisation, aber keine Hex-Spezifitaet
Der Mini-Main-Run zeigt bereits: Shuffled-Nullmodelle synchronisieren ebenfalls.
Damit ist J_theta als Reparatur mechanistisch plausibel, aber noch nicht als
hex-geometrischer Ordnungsmechanismus validiert.

## 2. R1mini Resultatstatus

- Baseline J_theta=0: reproduziert PHY023-Richtung.
- J_theta > 0: starke lokale Phasenausrichtung.
- Shuffled Null: ebenfalls starke Synchronisation.
- Status: REVIEW-POSITIVE fuer lokale XY-Synchronisation, kein PASS fuer BKT/Hex-Ordnung.

## 3. Neue Pflichtverbesserungen fuer R2

### 3.1 Echter Helicity-Modulus
Implementiere eine Twist-basierte Steifigkeitsmessung:
- Variante A: numerischer Twist-Test: Energie/Funktionswert bei Rand-Twist +/- delta
- Variante B: XY-Formel mit Richtungskomponente, falls periodische Randbedingungen verfuegbar

Auf offenen Patches ist Helicity-Modulus nur approximativ interpretierbar. Fuer harte
BKT-Claims sind periodische Honeycomb-Patches vorzuziehen.

### 3.2 Periodische Honeycomb-Geometrie
R2 muss neben offenen Hex-Patches auch periodische Honeycomb-/Torus-Geometrie haben.
Nur dann sind BKT-/FSS-Vergleiche fachlich naeher am Standard.

### 3.3 Vortex-Pairing statt nur Vortex-Dichte
R2 muss Vortex-Antivortex-Paarabstaende explizit messen:
- Anzahl +1 / -1 Plaquettes
- minimaler Abstand zu Gegenladung
- Pairing-Fraction unter Distanzschwelle
- Abhaengigkeit von J_theta und sigma

### 3.4 Random-Graph-Nullmodell haerter machen
Neben degree-preserving shuffled graph:
- gleiche Nodezahl/Edgezahl random graph
- gleiche Degree-Sequenz
- geometric random graph
- rotated/perturbed geometry
- edge-weight shuffle

### 3.5 Fit-Validierung
Power-vs-exp Fit nur akzeptieren mit:
- Mindestanzahl positiver Distanzen
- Bootstrap ueber Seeds
- AIC/BIC Delta
- Residualdiagnostik
- Finite-size Trend ueber radius/torus size

## 4. R2-Gate

R2 darf nur PASS-Kandidat werden, wenn:
1. J_theta=0 Baseline wieder negativ ist.
2. Hex/Torus zeigt staerkeres oder qualitativ anderes g1(r)-Verhalten als alle Nullmodelle.
3. Vortex-Antivortex-Pairing steigt im selben Parameterfenster.
4. Fit-Diagnostik ist ueber mindestens drei Systemgroessen stabil.
5. Effekt bleibt bei mehr Seeds und Step-Halving erhalten.
6. Backreaction lambda verbessert echte Ordnungsmetriken, nicht nur Lifetime/Proxywerte.

## 5. Claim-Status nach R1mini

- "XY-Kopplung erzeugt lokale Synchronisation": ERLAUBT.
- "XY-Kopplung repariert Phi-Hex vollstaendig": NICHT ERLAUBT.
- "Hex-Geometrie ist kausal fuer Ordnung": NICHT BELEGT.
- "BKT-kompatible Indikatoren sind sichtbar": REVIEW, nicht PASS.
- "Echter BKT-Uebergang": VERBOTEN.
