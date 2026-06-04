# Phi-Hex Meta-Audit (PHY017-PHY025): Finale Code-Konsolidierung und Experten-Roadmap

## TL;DR

- Die im Projekt beobachtete BKT-Crossing-Temperatur driftet korrekt von oben (1.540 bei L=5, 1.496 bei L=7, 1.490 bei L=9) gegen den hochpraezisen Literaturwert T_BKT(triangular) = 1.418(2) J/k_B; der naive Helicity-Crossing-Wert ist systematisch zu hoch und MUSS via Nelson-Kosterlitz-Sprung plus Weber-Minnhagen-Subleading-Log und der Paar-Extrapolation (Hsieh/Sandvik, L1,L2) in den thermodynamischen Limes gebracht werden - ein blosses Polynomfit der Crossing-Punkte genuegt nicht.
- Die wirksamste Einzel-Verbesserung des gesamten Stacks ist der Wechsel vom lokalen Update auf den Wolff-Embedding-Cluster-Algorithmus (dynamischer kritischer Exponent z praktisch null statt z ca. 2 bei Metropolis), kombiniert mit deutlich groesseren Gittern (L bis 128/256); das ist die notwendige Voraussetzung jeder Praezisions-Bestimmung und macht die heutige Begrenzung auf L<=9 obsolet.
- Methodisch sind die zwei wirksamsten Haerten genau auf die zwei im Projekt entlarvten Scheinergebnisse zugeschnitten: Blind-Analyse mit zielwert-unabhaengiger Vorwaerts-Normierung (gegen die rueckwaerts aus dem Zielwert abgeleitete Bond-Normierung) und Method-of-Manufactured-Solutions plus Regressions-Gates gegen analytische Referenzwerte (gegen das Rausch-Flacker-Metrik-Artefakt).

## Key Findings

### A) Praezisions-Bestimmung der kritischen BKT-Temperatur

**Referenzwert und Einordnung des Coworker-Resultats.** Der derzeit praeziseste publizierte Monte-Carlo-Wert fuer das ferromagnetische XY-Modell auf dem Triangular-Lattice ist T_BKT = 1.418(2) (Einheiten J/k_B), bestimmt von A. O. Sorokin, “Critical density of topological defects upon a continuous phase transition”, Annals of Physics (Amsterdam) 411, 167952 (2019), DOI 10.1016/j.aop.2019.167952. Die Zuordnung ist secondhand belegt in Otsuka/Shiina/Okabe, J. Phys. A 56, 235001 (2023), arXiv:2305.00651, woertlich: “By using the MC study of the helicity modulus, the BKT temperature was estimated as 1.418(2).”  Sorokin bestimmte den Wert ueber die Helicity-Modulus-Methode (universeller Sprung) mit Finite-Size-Scaling. Damit ist der Phi-Hex-Befund qualitativ korrekt: die Crossing-Temperatur 1.540 -> 1.496 -> 1.490 (L=5,7,9) driftet erwartungsgemaess von oben gegen 1.418. Das Vorzeichen der Drift (von oben) ist genau das, was die Theorie verlangt, weil endliche Gitter den Sprung verschmieren und das Crossing nach oben verschieben.

**Verhaeltnis triangular/square.** Der Triangular-Wert ist “roughly 3/2” des Square-Wertes; woertlich in arXiv:2305.00651: “The values of the triangular lattice are roughly 3/2 of those of the square lattice. It is because the coordination number of the triangular lattice is six, whereas [the square lattice is four].”  Der praezise Square-Referenzwert ist T_BKT = 0.89294(8), also beta_KT = 1.1199(1), aus Hasenbusch, “The two dimensional XY model at the transition temperature: A high precision Monte Carlo study”, arXiv:cond-mat/0502556 (Single-Cluster bis L=2048).  Die in der Folgeliteratur (Hsieh et al. 2013) verwendete Zahl 0.8935(1) ist mit diesem Wert konsistent. Das Verhaeltnis 1.418/0.8929 ca. 1.59 liegt zwischen dem reinen Koordinationszahl-Verhaeltnis 6/4 = 1.5 und dem tatsaechlichen Wert - die “3/2”-Regel ist also nur eine grobe Heuristik, kein exaktes Gesetz.

**Cluster-Update-Algorithmen.** Das zentrale Problem nahe T_c ist das critical slowing down (Autokorrelationszeit tau ~ L^z mit z ca. 2 fuer lokale Metropolis-Updates). Loesung:

- **Wolff-Embedding (Single-Cluster).** Man waehlt eine zufaellige Einheits-Richtung r, projiziert jeden kontinuierlichen Spin auf das Vorzeichen seiner r-Komponente und erhaelt so eingebettete Ising-Variablen. Auf diesen baut man einen Fortuin-Kasteleyn-Cluster mit Aktivierungswahrscheinlichkeit p_ij = 1 - exp(min(0, -2 beta J (s_i.r)(s_j.r))) und spiegelt den gesamten Cluster an der zu r senkrechten Ebene. In 2D fuer das XY-Modell ist das critical slowing down praktisch vollstaendig eliminiert (z nahe null). Hasenbusch erreichte damit L bis 2048. 
- **Swendsen-Wang-Multi-Cluster.** Identifiziert in einem Sweep alle FK-Cluster der eingebetteten Ising-Spins und flippt jeden unabhaengig mit Wahrscheinlichkeit 1/2. In 2D XY ist die Embedded-SW-Variante etwas langsamer (z_chi ca. 1) als Single-Cluster; Caracciolo et al. zeigten, dass die SW-Performance an der Frustration der induzierten Ising-Cluster haengt.  Vorteil von SW: alle Cluster werden pro Sweep aktualisiert, was sich gut parallelisieren laesst und natuerlich mit Parallel-Tempering / Replica-Exchange kombinierbar ist (so machten es Okabe-Otsuka).
- **Empfehlung Phi-Hex:** Single-Cluster-Wolff als primaeres Update, optional SW-Multi-Cluster fuer GPU-Parallelisierung. Wichtig: Der XY-Kopplungsterm sin(theta_j - theta_i) bzw. cos(theta_i - theta_j) (Durchbruch in PHY023) ist die Voraussetzung dafuer, dass das Embedding ueberhaupt physikalisch sinnvolle Cluster erzeugt - amplitudenbasierte Kopplung erzeugt keine BKT-Physik.

**Korrekte FSS-Protokolle (triangular-spezifisch).**

- **Nelson-Kosterlitz-Kriterium:** rho_s(T_BKT) = (2/pi) T_BKT (universeller Sprung).  Bei endlichem L verschmiert der Sprung; das naive Crossing von Upsilon(T) mit der Geraden (2/pi)T liefert ein nach oben verzerrtes Pseudo-T_c (genau die beobachtete Drift).
- **Weber-Minnhagen-Methode:** Statt naivem FSS nutzt man die Kosterlitz-RG-Gleichung. Am kritischen Punkt gilt pi*Upsilon/(2T) = 1 + 1/(2 ln(L/L0)), also eine logarithmisch langsame Annaeherung mit fixiertem Vorfaktor A(T_BKT)=1.  Man bestimmt T_BKT als jene Temperatur, bei der dieser Log-Fit am besten passt (chi-Quadrat-Minimierung).
- **State-of-the-art-Paar-Methode (Hsieh et al. 2013; Sandvik-Variante, arXiv:1302.2900):** Man definiert ein groessenabhaengiges T_BKT(L1,L2) aus einem Gittergroessen-Paar (z.B. L2=2*L1), das die Single-Parameter-Log-Korrektur der Stiffness bereits enthaelt, und extrapoliert mit der naechst-hoeheren Log-Korrektur in den Limes.  Diese Methode ist nachweislich gut konditioniert und wurde fuer das Standard-XY-Modell auf L bis 512 (GPU) verifiziert;  die Subleading-Logs haben “significant effects on the extrapolation”. 
- **Vierter-Ordnung-Helicity-Modulus Upsilon_4:** Hat ein Minimum am BKT-Punkt, das mit L monoton (divergent) waechst;  nuetzlich als unabhaengiger zweiter Schaetzer (Minnhagen-Kim 2003).
- **Logarithmus-arme Alternativen:** Correlation-Ratio R(T) = <g(L/2)>/<g(L/4)> und die zweite-Momente-Korrelationslaenge xi_2nd/L sind weniger sensitiv auf multiplikative Log-Korrekturen als die Binder-Kumulante (die fuer BKT generell schwach ist). Okabe-Otsuka nutzen genau R(T) mit Data-Collapse fuer das triangular Gitter. 

**Welche Gittergroessen und Statistik fuer welche Praezision.** Aus der Literatur: Olsson/Schultka-Manousakis erreichten mit L bis 256-400 und nur leading Logs eine Praezision von ca. 0.1% in T_BKT, aber mit Restdiskrepanz wegen vernachlaessigter Subleading-Logs.  Hasenbusch musste auf L=2048 gehen, um diese Diskrepanz “by brute force” aufzuloesen  (Ziel: 5. Dezimalstelle). Sandvik (arXiv:1302.2900) erreichte mit L bis 512 und der Paar-Methode robuste Extrapolation.  **Konkrete Empfehlung Phi-Hex:** Fuer ca. 1% Genauigkeit reicht L = 16,24,32,48,64 mit n>=16 Seeds; fuer ca. 0.1% (vierte Stelle) braucht es L bis 128/256 und die Subleading-Log-Extrapolation. Die heutige Begrenzung L<=9 ist fuer eine seriose T_BKT-Bestimmung definitiv zu klein - hier liegt der groesste Hebel.

**Coulomb-Gas-Mapping als alternativer Schaetzer.** Das XY-Modell mappt exakt auf ein 2D-Coulomb-Gas (Vortices = Ladungen +-, logarithmische Wechselwirkung) bzw. das Sine-Gordon-Modell  (Benfatto, Lecture Notes Cond-Mat 2024). T_BKT laesst sich dann unabhaengig ueber den dielektrischen Konstanten / die Vortex-Paar-Korrelationen schaetzen: man misst die effektive Coulomb-Kopplung aus den Ladungs-Dichte-Dichte-Korrelationen (Linear-Response, Orkoulas-Panagiotopoulos)  bzw. die magnetische Permeabilitaet (arXiv:1808.07183). Vorteil: vollstaendig unabhaengiger systematischer Fehler gegenueber der Stiffness-Methode - ideal als Cross-Validation-Gate.

### B) Defekt-/Vortex-Dynamik im geordneten Regime

**Quench durch den BKT-Uebergang.** Die Referenzarbeit ist A. Jelic & L. F. Cugliandolo, “Quench dynamics of the 2d XY model”, J. Stat. Mech. (2011) P02032, arXiv:1012.0417. Kernbefunde: Nach einem schnellen Quench aus der ungeordneten Phase in die Quasi-Long-Range-Order-Region waechst die Korrelationslaenge mit einer logarithmischen Korrektur zum diffusiven Gesetz;  woertlich: “The functional form is consistent with a logarithmic correction to the diffusive law and it serves to validate dynamic scaling… This analysis clarifies the different dynamic roles played by bound and free vortices.”  Die Autoren liefern eine Theorie der Quench-Raten-Abhaengigkeit (Kibble-Zurek), die ueber die Gleichgewichts-Scaling-Argumente hinausgeht.

**Coarsening/Annihilation (Modell A, T=0).** Qian & Mazenko, “Vortex Dynamics in a Coarsening Two Dimensional XY Model”, arXiv:cond-mat/0304346: Im nicht-erhaltenen O(2)-TDGL-Modell dominieren +-1-Vortices (hoehere Ladungen sind instabil); die Vortex-Geschwindigkeitsverteilung skaliert mit der mittleren Speed ~ t^-1/2 und zeigt einen algebraischen Tail mit Exponent -3 (bzw. -4 in der Geschwindigkeitsverteilung).  Die wachsende Laengenskala traegt eine logarithmische Korrektur zum diffusiven Wachstum: L(t) ~ (t/ln t)^1/2 (Community-Konsens-Form, vgl. arXiv:1211.1462)  - der Spezialfall n=d=2, fuer den die Bray-Rutenberg-Methode “mute” ist. 

**Langevin vs Monte-Carlo - wann was korrekt ist.** Die korrekte dynamische Universalitaetsklasse fuer einen nicht-erhaltenen Vektor-Ordnungsparameter ohne Erhaltungssaetze ist Modell A (Halperin-Hohenberg), realisiert durch die time-dependent Ginzburg-Landau-Langevin-Gleichung. Entscheidungsregel:

- **MC-Dynamik (Glauber/Metropolis Single-Spin-Flip):** Liegt in derselben dynamischen Universalitaetsklasse (Modell A) fuer kritische Relaxation und Coarsening, ist also fuer universelle Exponenten und Scaling-Funktionen korrekt. ABER: MC-“Zeit” (Sweeps) ist keine physikalische Zeit und kein deterministischer Trajektorienverlauf - fuer Vortex-Bahnen, Driftgeschwindigkeiten oder Lebensdauern in physikalischen Einheiten ungeeignet. WICHTIG: Cluster-Updates (Wolff/SW) zerstoeren die physikalische Dynamik voellig (nichtlokale Spruenge) und duerfen NUR fuer Gleichgewichts-Sampling, NIE fuer Dynamik-Studien verwendet werden.
- **Langevin/TDGL:** Liefert physikalische Zeit, echte Vortex-Trajektorien und Geschwindigkeitsverteilungen. Erforderlich, sobald Driftgeschwindigkeiten, Lebensdauern oder Backreaction quantitativ in physikalischen Einheiten gefragt sind. Caveat: gitterabhaengige UV-Divergenzen des Rauschens erfordern Renormierung der Gitter-Langevin-Gleichung. 

**Statistisch saubere Messung im geordneten Regime.** Vortex-Identifikation: Plaquette-Wirbelstaerke (Summe der auf [-pi,pi) reduzierten Phasendifferenzen um eine Elementar-Plaquette, /2pi); auf dem Triangular-Lattice besteht jeder Vortex-Kern aus drei Spins pro Elementardreieck.  Im geordneten Regime ist die freie Vortex-Dichte exponentiell unterdrueckt (n_v ~ exp(-E_c/T)), gebundene Paare dominieren. Saubere Praxis: (i) gebundene von freien Vortices trennen (Paar-Abstands-Cutoff / Coarse-Graining wie bei der 2D-Bose-Gas-Analyse arXiv:0912.1675); (ii) Lebensdauern aus Survival-Funktionen / First-Passage-Statistik, nicht aus Momentaufnahmen; (iii) Paar-Korrelationen g_vv(r) ueber viele unabhaengige Seeds; (iv) explizit gegen Rausch-Flackern testen (siehe Methodik-Lehre). Genau hier entstand das entlarvte Metrik-Artefakt: eine scheinbare Defekt-Stabilisierung, die nur Rausch-Flackern war - die statistische Trennung gebunden/frei und die Survival-Analyse haetten das verhindert.

**Backreaction - was die Theorie ueber Defekt-Kontrolle sagt.** Es ist physikalisch moeglich, Vortex-Dichte/-Dynamik extern zu beeinflussen, aber die Mechanismen sind eng begrenzt: In getriebenen dissipativen Systemen lockt ein resonantes Treiben die Phase und unterdrueckt die Bildung topologischer Anregungen, waehrend ein “engineered driving”-Profil sie gezielt erzeugt  (Polariton-Superfluide, arXiv:1612.07028). In aktiven nematischen Fluiden lassen sich Defekte mit “active topological tweezers” gezielt transportieren und braiden  (PNAS 2021, doi:10.1073/pnas.2400933121). ABER: Im strikten Gleichgewichts-geordneten Regime unterhalb T_BKT ist die freie Vortex-Dichte exponentiell unterdrueckt und durch Vortex-Antivortex-Bindung geschuetzt. Eine externe Stell-Groesse, die behauptet, Defekte zu “stabilisieren”, muss daher zwingend gegen den Bindungs-/Unbinding-Mechanismus und gegen reine Rausch-Effekte abgeglichen werden - eine echte Backreaction muss die effektive Stiffness oder das Vortex-Kern-Potenzial messbar verschieben, nicht nur die momentane Zaehlung. Das ist der theoretische Hintergrund, vor dem das Phi-Hex-Metrik-Artefakt zu Recht verworfen wurde.

### C) Workflow- und Methodik-Haertung (modell-agnostisch)

**Gegen Bestaetigungsverzerrung.**

- **Blind-Analyse (Teilchenphysik-Standard):** Roodman, “Blind Analysis in Particle Physics”, arXiv:physics/0312102; Klein & Roodman, Annu. Rev. Nucl. Part. Sci. 55, 141 (2005); MacCoun & Perlmutter, Nature 526, 187 (2015). Kernidee: Das Ergebnis (und/oder die Daten) bleiben verborgen, bis die Analyse “frozen” ist  - so kann der Zielwert die Methodenwahl nicht beeinflussen. Direkt anwendbar auf die entlarvte rueckwaerts-abgeleitete Bond-Normierung: Normierungen MUSS man vorwaerts aus der Theorie ableiten und einfrieren, bevor man den Zielwert (J*sqrt(3), 2/pi-Sprung, T_BKT) ueberhaupt anschaut. Das dokumentierte “stopping bias”-Phaenomen (man sucht Fehler nur bei ueberraschenden Ergebnissen)  ist exakt das Phi-Hex-Risiko.
- **Vorab-Registrierung fuer Simulationen:** ADEMP-PreReg-Template (Siepe, Bartos, Morris, Boulesteix, Heck, Pawel, “Simulation Studies for Methodological Research”, Psychological Methods 2024, DOI 10.1037/met0000695; GitHub bsiepe/ADEMP-PreReg). ADEMP = Aims, Data-generating mechanism, Estimands, Methods, Performance measures. Enthaelt Formeln fuer Monte-Carlo-Standardfehler und die noetige Wiederholungszahl.  Mit Zeitstempel auf OSF/Zenodo hinterlegen. 

**Verifikation & Validierung (V&V).**

- **Method of Manufactured Solutions (MMS):** Roache, “Verification and Validation in Computational Science and Engineering” (Hermosa 1998); Roache, J. Fluids Eng. 124, 4 (2002). Man waehlt eine analytische “manufactured solution”, setzt sie in die Gleichungen ein, leitet den noetigen Quellterm symbolisch ab und prueft, ob der Code bei systematischer Gitterverfeinerung die theoretische Konvergenzordnung erreicht (“theorem-like quality with a clearly defined completion point”).  Fuer Phi-Hex: analytische Referenzen wie die erreichte Maschinengenauigkeit 2e-16 der Madelung-Bohm-Fisher-Variationsableitung, die O(N^-2)-Konvergenz von Hellinger-Link vs Wasserstein-Maas, und der J*sqrt(3)-T=0-Wert sind genau solche manufactured-solution-Gates.
- **Grid Convergence Index (GCI):** Roache; Oberkampf & Roy, “Verification and Validation in Scientific Computing” (Cambridge 2010); Roy & Oberkampf, Comput. Methods Appl. Mech. Engrg. 200, 2131 (2011). Liefert eine einheitliche, sicherheitsfaktor-gewichtete Fehlerschranke aus Richardson-Extrapolation. Klare Trennung beibehalten: Code-Verifikation (Order-of-Accuracy, gegen exakte Loesungen) vs Solution-Verifikation (Diskretisierungs-Unsicherheit der konkreten Rechnung) vs Validierung (gegen Experiment/Natur).

**Reproduzierbarkeit.**

- **FAIR4RS:** Barker et al., “Introducing the FAIR Principles for research software”, Scientific Data 9, 622 (2022), DOI 10.1038/s41597-022-01710-x; Lamprecht et al., Data Science 3 (2020). Jede Code-Version bekommt eine eigene persistente ID (PID); Git-Commit-Hashes sind ein guter Anfang, aber nicht global aufloesbar  - daher Zenodo-DOI pro Release.
- **Deterministisches Seed-Mgmt + Provenance + Container:** Seeds als explizite, geloggte Parameter (n>=16, dokumentierte Generator-Familie); Provenance-Tools (z.B. noWorkflow fuer Python) erfassen Ausfuehrungsmetadaten; Container (Docker/Apptainer) garantieren bit-identische Reproduktion ueber Maschinen hinweg (in der Literatur war nur die VM/Container-Variante maschinenunabhaengig reproduzierbar). 

**Statistik-Standards (passend zu n>=16, nicht-normal).**

- **Cliff’s delta / Vargha-Delaney A:** Verteilungsfrei, robust fuer kleine bis mittlere (n=10-50) nicht-normale Stichproben (Delaney-Vargha 2002; Feng-Cliff 2004).  Schwellen (Vargha-Delaney): |delta| = 0.11 / 0.28 / 0.43 fuer klein/mittel/gross.  Cliff’s delta ist eine lineare Transformation von A. 
- **BCa-Bootstrap:** Liefert bessere CI-Coverage als Cohen’s d gerade bei heterogenen Varianzen und Nicht-Normalitaet (konsistent 0.95-Coverage).  Fuer alle Effektstaerken-CIs verwenden.
- **Hedges’ g** ergaenzend, wenn Normalitaet plausibel; bei kleinen n die Bias-Korrektur (g statt d) zwingend.

**Selbstkorrektur-Infrastruktur.**

- **Architecture Decision Records (ADRs):** Pro nicht-trivialer Methoden-/Normierungs-Entscheidung ein kurzes, versioniertes Markdown-Dokument (Kontext, Entscheidung, Konsequenzen, Status). Direkt verknuepfbar mit der Vorwaerts-Normierungs-Disziplin.
- **Negativ-Ergebnis-Register:** Die entlarvten Scheinergebnisse (Metrik-Artefakt, rueckwaerts-Bond-Normierung) explizit als datierte Negativ-Eintraege fuehren - das ist die institutionalisierte Form der Phi-Hex-Selbstkorrektur und verhindert Re-Litigation.

**Codebasis-Struktur fuer mehrstufige Simulationen.**

- **TDD fuer wissenschaftlichen Code:** Tests vor Code; jeder analytische Referenzwert wird zu einem automatisierten Regressionstest (z.B. 2/pi-Sprung, J*sqrt(3), 2e-16-Variationsableitung, O(N^-2)-Konvergenz, Liouvillian-Gap-Cross-Check Žnidarič 2015 / Mori-Shirai 2020). Unit-Tests muessen unabhaengig und schnell sein;  Warnung: bestandene Unit-Tests beweisen keine Korrektheit (nur Abwesenheit bekannter Regressionen). 
- **Kontinuierliche Validierungs-Gates (CI):** Bei jedem Commit laufen die Regressionstests gegen die analytischen Referenzwerte; ein “code bisector” lokalisiert eingefuehrte Regressionen automatisch.  Modularitaet: klare Trennung von Kernel (PHY017-Validierung), Simulator (v1.9.x), Stell-Schleife (PHY020) und Analyse (PHY023/025), damit jede Stufe einzeln gegen ihre Referenz testbar ist.

## Details

Die drei Bereiche haengen ueber eine gemeinsame Diagnose zusammen: Die Phi-Hex-Pipeline hat die richtige Physik (BKT auf Triangular mit echtem Phasen-Kopplungsterm) und die richtige Selbstkorrektur-Kultur bereits etabliert; was fehlt, ist (1) algorithmische Reichweite (Cluster-Update + grosse Gitter), (2) das korrekte Subleading-Log-Extrapolationsprotokoll und (3) die Formalisierung der bereits intuitiv praktizierten Skepsis in pruefbare Gates.

Zur Interpretation der Crossing-Drift: Die Werte 1.540/1.496/1.490 fuer L=5/7/9 zeigen das fuer BKT typische, extrem langsame (logarithmische) Annaeherungsverhalten. Eine naive lineare oder polynomielle Extrapolation in 1/L unterschaetzt systematisch die Kruemmung, weil die fuehrende Korrektur ~1/ln(L) und nicht ~1/L ist. Das erklaert, warum man ohne den Weber-Minnhagen/Hsieh-Apparat nicht zuverlaessig bei 1.418 landet - und warum so kleine Gitter (L<=9) prinzipiell unzureichend sind: bei L=9 ist ln(L) ca. 2.2, die Log-Korrektur also riesig.

Zur Backreaction-Frage konkret: Bevor irgendeine Stell-Groesse als “Vortex-Kontrolle” interpretiert wird, sollte ein Null-Test laufen (Stell-Groesse aus, gleiches Rauschen): wenn die gemessene Defekt-Statistik sich nicht signifikant (Cliff’s delta mit BCa-CI, n>=16) unterscheidet, ist der Effekt ein Artefakt. Eine echte Backreaction muss zudem die renormierte Stiffness (Helicity-Modulus) oder das Vortex-Kern-Potenzial verschieben - eine reine Aenderung der Momentan-Zaehlung genuegt nicht.

## Recommendations

**Stufe 1 - sofort (Algorithmus & Skalierung, groesster Hebel):**

1. Wolff-Single-Cluster-Embedding implementieren und gegen das bestehende lokale Update benchmarken (Autokorrelationszeit vs L messen, z bestimmen). Erwartung: z faellt von ca. 2 auf nahe 0.
1. Gittergroessen auf L = 16,24,32,48,64 ausweiten (n>=16 Seeds). Benchmark-Gate: T=0-Helicity-Modulus muss J*sqrt(3) reproduzieren (flaechen-normiert) - als automatischer Regressionstest.

- Schwelle: Sobald die L=16..64-Daten den Nelson-Kosterlitz-Log-Fit mit chi^2/dof ~ 1 erfuellen, weiter zu Stufe 2.

**Stufe 2 - T_BKT-Praezision:**
3. Weber-Minnhagen-Log-Analyse plus Hsieh/Sandvik-Paar-Extrapolation T_BKT(L1,L2) mit Subleading-Log anwenden. Ziel: T_BKT mit Fehlerbalken, der 1.418(2) einschliesst.
4. Falls 0.1%-Praezision gefordert: L bis 128/256 (GPU/SW-Multi-Cluster + Parallel-Tempering).
5. Coulomb-Gas-Kreuzcheck (dielektrische Konstante / Permeabilitaet) als unabhaengiger zweiter T_BKT-Schaetzer - Konsistenz beider Methoden ist das Akzeptanzkriterium.

**Stufe 3 - Dynamik (separater Code-Pfad):**
6. Fuer Vortex-Dynamik/Backreaction strikt von Cluster-Sampling trennen: TDGL/Langevin (Modell A) fuer physikalische Zeit und Trajektorien; MC-Glauber nur fuer universelle Scaling-Checks. Vortex-Lebensdauern via Survival-/First-Passage-Statistik, gebundene/freie Vortices trennen.
7. Backreaction nur akzeptieren, wenn (a) Null-Test (Stell-Groesse aus) signifikanten Unterschied zeigt (Cliff’s delta + BCa-CI) UND (b) die renormierte Stiffness messbar verschoben ist.

**Stufe 4 - Methodik-Haertung (parallel, dauerhaft):**
8. Blind-Analyse-Disziplin: alle Normierungen vorwaerts aus der Theorie ableiten und vor dem Zielwert-Abgleich einfrieren (ADR pro Normierung).
9. ADEMP-PreReg vor jeder neuen Hypothesen-Kampagne; Zeitstempel auf Zenodo.
10. MMS-Gates + CI-Regressionstests gegen alle analytischen Referenzen (2/pi, J*sqrt(3), 2e-16, O(N^-2), Liouvillian-Gap). Negativ-Ergebnis-Register und ADRs als Pflicht-Artefakte. FAIR4RS: Zenodo-DOI pro Release, Container, geloggte Seeds.

## Caveats

- Der Triangular-Referenzwert T_BKT = 1.418(2) ist hier secondhand bestaetigt (woertliches Zitat aus Otsuka/Shiina/Okabe 2023, arXiv:2305.00651); die Primaerquelle ist Sorokin, Ann. Phys. 411, 167952 (2019), DOI 10.1016/j.aop.2019.167952, die ich nicht direkt einsehen konnte. Die genauen von Sorokin verwendeten Gittergroessen sind nicht verifiziert. Es besteht zudem eine kleine [28]-vs-[29]-Referenz-Label-Ambiguitaet in arXiv:2501.07388, die bei formaler Zitation direkt geprueft werden sollte.
- Die “triangular = 3/2 * square”-Regel ist eine grobe Heuristik (tatsaechliches Verhaeltnis ca. 1.59), kein exaktes Gesetz; nicht fuer quantitative Vorhersagen verwenden.
- Die flaechen-normierte J*sqrt(3)-T=0-Stiffness ist als geometrischer Normierungsfaktor plausibel und konsistent mit dem analogen Honeycomb-Faktor 4/(3*sqrt(3)) (arXiv:2406.12076),  aber ich habe keine Primaerquelle gefunden, die exakt “triangular T=0-Stiffness = J*sqrt(3)” woertlich aussagt - vor Nutzung als hartes Gate intern analytisch verifizieren.
- L(t) ~ (t/ln t)^1/2 und der Geschwindigkeits-Tail-Exponent stammen aus T=0-Coarsening (Qian-Mazenko); im endlich-temperierten geordneten Regime modifizieren Spin-Wellen und thermische Vortex-Paare das Bild.
- Die zitierten honeycomb-Werte (0.571-0.576) weichen signifikant vom analytischen 1/sqrt(2) ab - ein in der aktuellen Literatur offener Punkt; er illustriert, dass selbst “exakte” analytische BKT-Vorhersagen kritisch gegen Hochpraezisions-MC zu pruefen sind (genau die Phi-Hex-Lehre).