# PHI HEX W4 — Honeycomb T_BKT Praeregistrierung v01

**Datum:** 2026-09-27  
**Status:** PRE-REGISTERED DESIGN / noch kein Produktionsresultat  
**Bezug:** GitHub Issue #45  
**Claim ceiling:** FINDING / CONSISTENCY CHECK; kein neuer T_BKT-Bestwert ohne
separates Cross-Family-Review und dokumentierte Finite-Size-/Modellwahl-Systematik.

## 1. Forschungsfrage

W4 versucht **nicht**, 0.573 gegen 0.576(3) als konkurrierende Wahrheiten
aufzuloesen. Diese Differenz ist kleiner als bzw. vergleichbar mit publizierten
Unsicherheiten und mit realistischen Finite-Size-/Modellwahl-Systematiken bei
Laptop-grossen Gittern.

Die belastbare Frage lautet:

> Ist eine unabhaengige PHI-Hex-Honeycomb-Schaetzung, unter explizit gemessener
> Sampler-Validitaet und Finite-Size-Systematik, mit den aktuellen
> Literaturergebnissen vereinbar oder zeigt sie reproduzierbare Spannung?

## 2. Vor dem Lauf fixierte Literaturbasis

Vertragsquelle: `spec/260703 PHI HEX honeycomb reference conventions audit v01.md`.

Primaerset:
- Okabe/Otsuka, arXiv:2501.07388v1 / J. Phys. A 2025: rough estimate 0.573.
- Jiang, PTEP 2024 103A02, publizierte Fassung: NN 0.572(3), Helicity 0.576(4).
- de Andrade/Jorge/DaSilva, arXiv:2406.12076v4: Upsilon SA 0.575(8),
  Upsilon WL 0.576(3).

Diagnostik:
- de Andrade/Jorge/DaSilva v4 Upsilon_4: SA 0.551(11), WL 0.568(1).
  Dieser Kanal wird wegen der im Paper selbst diskutierten Minimum-/FSS-
  Systematik nicht mit Y2 zu einem einzigen "Band" vermischt.

Supersedierte arXiv-v1/beta-Werte duerfen nur als Historie erscheinen.

## 3. Design-Gates vor jeder Grossproduktion

### W4-G0 — Walker-Identifizierbarkeit

Jede Gittergroesse, die in einen W4-Paar-Schaetzer eingeht, braucht
**mindestens 3 unabhaengige WL-Walker**. Ein oder zwei Walker duerfen als
Kalibrierung laufen, erzeugen aber **keine W4-Validitaetsfreigabe**.

Ein Ein-Walker-Spread wird weder als 0 noch als "voll gueltige Domaene"
serialisiert, sondern als `null / unmeasured`.

### W4-G1 — L=64 Kalibrierung vor Skalierung

Vor L=96/128:
- L=64, 2 Walker nur zur Kosten-/Konvergenzkalibrierung,
- reale `wl_sweeps`, Walltime, Binzahl, Produktionsdichte, Rand-Leak,
  unbesetzte kanonische Masse und Y2-Walker-Spread erfassen,
- daraus Kostenexponent **messen**, nicht aus L^p annehmen.

Die 2-Walker-Kalibrierung darf W4-G0 nicht ersetzen.

### W4-G2 — gemeinsame Validitaetsdomaene

Ein Paar (La,Lb) ist nur quotierbar, wenn:
1. beide L jeweils eine gemessene Walker-Domaene haben,
2. das Crossing innerhalb beider Domaenen liegt,
3. Rand-Leak und unbesetzte kanonische Masse die bestehenden 1e-3-Gates
   erfuellen,
4. der Root-Finder nur benachbarte gueltige T-Punkte verbindet.

Alternative Evidenz (z.B. PHY032-Drift-Guard) muss als eigene Basis
ausgewiesen werden und darf nicht als Walker-Spread umetikettiert werden.

### W4-G3 — Mindestinformation

Ein W4-Literaturvergleich ist `INCONCLUSIVE`, wenn mindestens eine Bedingung
gilt:
- weniger als 3 quotierbare L-Paare,
- kombinierte 1-sigma Unsicherheit der W4-Zusammenfassung > 0.010,
- kein Modellwahlvergleich zwischen mindestens zwei plausiblen
  BKT-FSS-Ansatzen,
- das groesste L liegt ausserhalb seiner gemessenen Sampler-Domaene am
  relevanten Crossing.

## 4. Finite-Size- und Modellwahl-Systematik

Mindestens zwei Auswertungen werden parallel berichtet:

A. bestehender C-eliminierter Weber-Minnhagen-Paarpfad;  
B. logarithmische BKT-FSS mit subleading-Korrektur / alternativer zulässiger
   Form, analog der Sensitivitaetslogik von Hsieh-Kao-Sandvik
   (arXiv:1302.2900).

Zu berichten:
- Fitfenster / L_min,
- Schaetzer pro Modell,
- statistische Unsicherheit,
- Modellwahl-Spread,
- Residuen und Stabilitaet bei L_min-Aenderung.

Die Modellwahl-Komponente wird **nicht** durch hoehere Walkerzahl ersetzt.

## 5. Entscheidungsregel gegen Literatur

Kein kuenstliches Min/Max-"Wahrheitsband".

Wenn W4-G0..G3 bestanden sind:

- **CONSISTENT_CURRENT_LITERATURE:** das W4-Unsicherheitsintervall ist mit
  mindestens einem aktuellen Ergebnis aus **Jiang/PTEP** und mindestens einem
  aktuellen Ergebnis aus **de Andrade et al. v4** vereinbar; Okabe/Otsuka 0.573
  wird als dritter Punkt ohne erfundene Fehlerbar berichtet.
- **TENSION_CURRENT_LITERATURE:** Informations-Gates bestanden, aber keine
  solche Zwei-Quellen-Kompatibilitaet.
- **INCONCLUSIVE:** Informations-/Sampler-/Modellwahl-Gates nicht bestanden.

Diese Kategorien bewerten Datenkonsistenz, nicht die "richtige" Publikation.

## 6. Stop-/Scale-Regel nach L=64

Nach der L=64-Kalibrierung:
- Produktionsbudget fuer L=96/128 nur freigeben, wenn gemessene Kosten und
  Konvergenz ein vollstaendiges W4-G0..G3-Design realistisch machen.
- Andernfalls NEGATIVE_RESULT/INCONCLUSIVE dokumentieren und auf effizientere
  Samplingstrategie oder HPC/GPU-Pfad wechseln, statt unterbudgetierte
  Gross-L-Werte zu produzieren.

## 7. Reproduzierbarkeit

Jeder Lauf protokolliert mindestens:
- Git-Commit,
- Python/NumPy-Version,
- master_seed und alle Stream-IDs,
- L, Walkerzahl, T-Gitter, Energie-Fenster, Binzahl,
- `lnf_final`, `wl_sweeps`, `prod_sweeps`, Walltime,
- Leak/Coverage-Diagnostik,
- Walker-Spread und Domaenenbasis,
- alle L-Paar-Crossings mit Quotierbarkeitsgrund,
- alle Negativ-Results.

Historische PHY041/PHY042-Ergebnisdateien bleiben unveraendert; Errata werden
append-only in Specs/CHANGELOG dokumentiert.
