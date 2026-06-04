# Meta-Analyse PHY010-016 — Phi-Hex-Theorie: Konstruktal-Gesetz, Fisher-Information, Madelung-Hydrodynamik, diskrete Topologie und Hydraulik

**Konsolidierungsdokument, Stand 26. Mai 2026**
Coworker Research / Coworkerz Attribution — Swiss-Konventionen, ASCII-Quotes, “ss” statt scharfes S, CHF, Englische Fachterme behalten.

-----

## TL;DR

- Die Bruecke Konstruktal-Law ↔ Fisher-Information ↔ Madelung-Hydrodynamik ↔ PMPG ist mathematisch fragmentarisch: die Identitaet ⟨Q_Bohm⟩ = (ℏ²/8m)·I_F gilt rigoros im Kontinuum (Heifetz-Cohen, Found. Phys. 45, 1514, 2015; Bloch-Cohen, arXiv:2210.07732, 2022), aber **keine publizierte Theorie** behandelt diese Identitaet auf einem N=6-Tight-Binding-Ring; die Diskretisierung kann je nach Operator-Wahl O(N⁻²) (Drei-Punkt-Laplace), exakt (Hellinger-Link-Form gemaess Sason et al. arXiv:1904.12704) oder exponentiell schnell (spektrale Quadratur via Trefethen-Weideman, SIAM Rev. 56, 385, 2014) konvergieren. PHY010-016 muss diese Konvergenz explizit testen statt zu postulieren.
- Mean-Field-Bifurkationen auf einem N=6-Bose-Hubbard-Ring sind generisch **Artefakte**: Single-Site-Gutzwiller bricht qualitativ, sobald Quantenfluktuationen ~ 1/sqrt(z) gross werden (Koordinationszahl z=2 auf dem Ring); Cluster-Gutzwiller (Luehmann, Phys. Rev. A 87, 043619, 2013; arXiv:1302.6761)  und exakte Diagonalisierung (Hilbert-Dimension fuer 6 Bosonen / 6 Sites ≈ 462 nach Zhang-Dong, arXiv:1102.4006) sind die Pflicht-Werkzeuge — PHY010 haette mit Cluster-Gutzwiller oder ED den Phasenuebergang nicht als scharfe Bifurkation, sondern als Crossover identifiziert.
- Die Wallstrom-Quantisierung ist auf einem C_6-Ring mit Peierls-Phase **nicht trivial**, sondern manifestiert sich als sechs diskrete Bloch-Klassen k ∈ ℤ/6ℤ mit Spektrum E_k = −2t·cos[2π(k+φ/φ_0)/6] (Maiti, arXiv:0706.0061; Heinrichs, arXiv:cond-mat/0106437) und Persistent-Current-Periode φ_0 = h/e (Byers-Yang-Theorem); die Sprint-Behauptung “Wallstrom auf N=6 trivial” ist zu revidieren.

-----

## Key Findings

### Achse 1 — Konstruktal ↔ Fisher ↔ Madelung ↔ PMPG

**Stand der Wissenschaft (peer-reviewed, 2022-2026)**

- **Reddiger-Poirier 2023** (J. Phys. A 56, 193001; arXiv:2207.11367) liefert die bislang strengste mathematische Aufarbeitung des Madelung-Wallstrom-Problems. Fuenf Kern-Resultate:
  (i) Takabayasi-Quantisierungsbedingung ∮∇S·dl ∈ 2πℏ·ℤ gilt rigoros fuer alle C¹-Wellenfunktionen;
  (ii) Madelung muss distributional verstanden werden (weak formulation gemaess Gasser-Markowich oder Nelson-Equation);
  (iii) Einfuehrung “quantum quasi-irrotationality” zur Klassifizierung der Faelle in denen die distributionalen Madelung-Gleichungen nicht erfuellt sind;
  (iv) explizite Konstruktion nicht-quantisierter starker Loesungen in 2D;
  (v) Wallstrom-Phaenomen entspringt der “failure of quantum mechanics to discern physically equivalent, yet mathematically inequivalent states — an issue that finds its historic origins in the Pauli problem”. 
  **Diskrete Pendants existieren nicht; mehr mathematische Forschung explizit als noetig markiert.** 
- **PMPG (Principle of Minimum Pressure Gradient)**: Gonzalez-Taha 2022 und Taha et al. 2023 etablieren PMPG als Variationsprinzip basierend auf Gauss’ Principle of Least Constraint;  Taha 2025/2026 (arXiv:2602.20637) zeigt formale Zwei-Wege-Aequivalenz zwischen INSE und PMPG: “a candidate smooth flow field is a solution of the INSE if and only if its instantaneous evolution minimizes, at every instant, the norm of the pressure force, required to enforce incompressibility”.  PMPG ist im Wesentlichen die Minimierungs-Formulierung der Leray-Helmholtz-Projektion.  **PMPG ist klassisch-hydrodynamisch — kein eigenstaendiges Quanten-PMPG mit publizierter Variations-Aequivalenz zu Schrödinger.**
- **Mauri-Giona 2025**: im Suchergebnis nicht direkt auffindbar; vermutlich Bezug auf Mauri-Werke zur stochastischen Hydrodynamik und thermodynamischen Fluktuationen — **UNBEKANNT/KONFLIKT**, sollte aus dem Sprint-Material direkt referenziert werden.
- **Reginatto 1998** (Phys. Rev. A 58, 1775) bleibt die kanonische Ableitung der Schrödinger-Gleichung aus dem Prinzip minimaler Fisher-Information; Hall-Reginatto 2002 und neuere Arbeiten (Yang 2023, arXiv:2302.14619 “Quantum Mechanics Based on an Extended Least Action Principle and Information Metrics of Vacuum Fluctuations”)  verallgemeinern auf nicht-lineare Schrödinger-Gleichungen via Rényi- oder Kullback-Leibler-Divergenzen. **Reginatto-Konstruktion ist kontinuierlich; auf einem N=6-Ring noch nicht rigoros publiziert.**
- **Cafaro-Alsing 2020** (Fisher-Geodaesiken auf Mannigfaltigkeiten parametrisierter Quantenzustaende) liefert Information-Geometric-Speeds und entropy production rates fuer su(2,C)-Hamiltonians, ist aber kontinuumsbasiert.

**Identitaet ⟨Q_Bohm⟩ = (ℏ²/8m)·I_F: Konvergenzanalyse (zentrale Vertiefung)**

- Im **Kontinuum** rigoros (Heifetz-Cohen, Found. Phys. 45, 1514, 2015): aus Q = −(ℏ²/2m)·∇²√ρ/√ρ und partieller Integration ergibt sich
  ⟨Q⟩ = ∫ρQdx = (ℏ²/2m)∫|∇√ρ|²dx = (ℏ²/8m)·∫(∇ρ)²/ρ dx = (ℏ²/8m)·I_F[ρ].
- **Diskreter Fall**, konsolidiert aus konzentrierter Recherche:
  - Natuerliche diskrete Form auf N-Ring: I_F^{discrete}[ρ] = (4/a²)·Σ_{i=1}^N(√ρ_{i+1} − √ρ_i)² (Sason et al., arXiv:1904.12704, “Cramér-Rao-type Bound and Stam’s Inequality for Discrete Random Variables”, 2019/2020). Diese ist aequivalent zu 8·H²(p, T₊p) mit Hellinger-Abstand H.
  - **Konvergenz haengt vom diskreten Operator ab**:
    - Drei-Punkt-Laplace (semiklassisch, Iannaccone-Curatola-Fiori, “Effective Bohm Quantum Potential for device simulators”, 2004) → lokale Trunkierungs-Ordnung O(h²) = O(N⁻²) bei smooth √ρ.
    - Link-Hellinger-Form → Identitaet exakt fuer beliebiges N (per Definition); Vergleich zum Kontinuums-Integral hat Quadraturfehler.
    - Glatte ρ auf periodischer Domaene + Periodisches Trapez-Verfahren (Trefethen-Weideman, SIAM Rev. 56, 385, 2014) → spektrale Konvergenz, Fehler ~ exp(−cN). Schon fuer N=6 mit glatter Grundzustands-Dichte sub-Prozent erwartbar.
    - Yepez-Boghosian (Comput. Phys. Commun. 146, 280, 2002): numerisch beobachtete dritte bis vierte Konvergenzordnung im Quantum-Lattice-Gas.  
- **Keine publizierte Arbeit** quantifiziert ⟨Q_Bohm⟩_N − (ℏ²/8m)·I_F^N auf einem Tight-Binding-Ring mit Peierls-Phase als Funktion von N mit bewiesenem O(N⁻ᵏ); dies ist ein **offenes Forschungsproblem** und damit ein konkreter Phi-Hex-Beitragspunkt.

**Brueche bei N=6, “67%-Identitaets-Verletzung”:** Plausibel als Operator-Mismatch — wenn die Sprint-Implementierung den Drei-Punkt-Laplace gegen die Hellinger-Form testet, ist eine grosse Diskrepanz zu erwarten, weil die Operatoren *nicht denselben Kontinuums-Limes* haben. Best Practice: Beide Formen separat dokumentieren, kontinuum-extrapolieren via Richardson, und mit dem analytischen Wert auf exakt diagonalisierter Grundzustands-Dichte vergleichen.

**Verdikt fuer Phi-Hex**: Die Bruecke ist konzeptionell wertvoll, aber **N=6 ist zu klein fuer eine rigorose Konvergenz-Demonstration**; Sprint sollte Skalierungs-Studie N ∈ {6, 12, 24, 48, 96} mit beiden Operator-Familien durchfuehren.

### Achse 2 — Mean-Field vs. exakte Diagonalisierung in N=6-Quantensystemen

**Stand 2024-2026**:

- **Single-Site-Gutzwiller** ist auf einem N=6-Ring qualitativ unzuverlaessig fuer Bose-Hubbard: Quantenfluktuationen pro Site skalieren wie 1/sqrt(z) mit Koordinationszahl z; auf einem Ring ist z=2, deutlich kleiner als im 3D-Kubus (z=6) wo MF qualitativ ok ist. Phase-Boundaries sind systematisch verschoben (Luehmann, Phys. Rev. A 87, 043619, 2013).
- **Cluster Gutzwiller** (Luehmann, arXiv:1302.6761): bei genuegend grossen Clustern naehert sich QMC;  auf einem Sechs-Ring ist der ganze Ring der Cluster — exakte Diagonalisierung ist effizienter.
- **Exakte Diagonalisierung** mit Zhang-Dong 2011 (arXiv:1102.4006) ist Standard-Pedagogik;  fuer 6 Bosonen auf 6 Sites mit Cutoff n_max = 6 ist die Hilbert-Dimension C(11,5) = 462, problemlos in Sekunden auf einem Laptop.
- **Spin-Systeme**: Heisenberg/XXZ sind oft gutmuetiger fuer MF in höheren Dimensionen; im 1D-Ring sind sie via Bethe-Ansatz exakt loesbar, MF qualitativ falsch (etwa bei gapless XXZ).
- **Cluster-Gutzwiller-Ladder-Studien** (Yamamoto et al. 2015, Phys. Rev. A 92, 023618; arXiv:1504.08192) demonstrieren explizit, dass fractional insulator phases im Single-Site-Gutzwiller fehlen und nur im Cluster-Gutzwiller-Bild erscheinen. 

**Goldene Regeln 2024-2026 fuer MF-Validitaet**:

1. Koordinationszahl z gross (z ≥ 4 fuer Bosonen, z ≥ 6 fuer Fermionen mit Coulomb).
1. Ordnungsparameter nicht-fluktuierend (broken-symmetry Phase weit von Quanten-Phasenuebergang).
1. Keine Frustration, keine entartete Grundzustands-Mannigfaltigkeit.
1. Korrelationslaengen ξ << System-Linear-Dimension.

**Auf N=6 sind alle vier Bedingungen verletzt** — MF-Bifurkation ist generisch Artefakt.

**Beyond-MF-Werkzeuge**:

- DMRG (1D, fuer Ring mit periodischen Randbedingungen langsamer als open boundary; 2024 zahlreiche tensor-network Erweiterungen).
- Quantum Monte Carlo (worm algorithm fuer Bose-Hubbard).
- Exakte Diagonalisierung (fuer N ≤ 12 problemlos).
- Time-dependent variational principle (TDVP) mit MPS.

**PHY010-Korrektur**: Die im Sprint vermutete “Mean-Field-Bifurkation” ist mit hoher Wahrscheinlichkeit Artefakt; mit Cluster-Gutzwiller (ganzer Ring als Cluster) oder ED entstehen statt scharfer Phasenuebergaenge **glatte Crossovers** mit endlicher Korrelationslaenge.

### Achse 3 — Wallstrom-Quantisierung und topologische Konsistenz

**Status quo**:

- Wallstrom 1994 (Phys. Rev. A 49, 1613): Madelung ≠ Schrödinger ohne Quantisierungs-Postulat fuer ∮∇S·dl.
- Reddiger-Poirier 2023 (siehe oben): der Einwand ist im distributionalen Sinn nicht zwingend, aber Madelung ist nicht trivial wohldefiniert.
- Schmelzer 2011 (arXiv:1101.5774): Postulat ueber Δρ-Regularitaet an Knoten impliziert die Quantisierungsbedingung. 
- Gay-Balmaz-Tronci 2023 (arXiv:2003.08664, “Holonomy and vortex structures in quantum hydrodynamics”): U(1)-Bundle/Holonomie-Formulierung — “the quantization condition of the circulation is … another example of holonomy … this is in fact a type of monodromy, with the exact value depending on the winding number of the loop surrounding the singularity”. Winding numbers nicht-trivial ganzzahlig.
- Holm-Rawlinson-Tronci 2020 (J. Chem. Phys. 153, 234108; arXiv:2012.03569): “bohmion method” als regularisierte Bohm-Potential-Diskretisierung mit δ-train solutions.

**Auf C_6-Ring mit Peierls-Phase**:

- Tight-Binding-Spektrum (Maiti arXiv:0706.0061, Heinrichs arXiv:cond-mat/0106437):
  **E_k = −2t·cos[2π(k + φ/φ_0)/N]** fuer k ∈ {0,…,N−1}, φ_0 = h/e.
- **Sechs distinkte Bloch-Klassen**, jede mit nicht-trivialer ganzzahliger Winding number w=k.
- Persistent current I(φ) = −∂E/∂φ ist periodisch in φ mit Periode φ_0 (Byers-Yang-Theorem 1961).
- Aharonov-Bohm-Phase φ verschiebt die Labels k → k + φ/φ_0 kontinuierlich, aendert aber nicht die Klassenzahl.
- Wallstrom auf dem Ring: ∮∇S = 2π·k fuer integer k ∈ {0,…,N−1}.

**Hexagonale topologische Materialien 2024-2026**:

- Honigwaben-Modelle (Haldane, Kane-Mele) auf N=6-Ring: Chern-Zahl ist Eigenschaft der unendlichen 2D-Brillouin-Zone, nicht direkt auf einem einzelnen Hexagon definiert. Bott-Index ist die endliche-System-Verallgemeinerung.
- Lieb-Nachtergaele (arXiv:cond-mat/9410100): rigorose Peierls-Instabilitaet fuer L=4k+2-Ringe (relevant fuer N=6).
- Gontier-Roussigné-Séré, arXiv:2510.24230 (Okt 2025): Peierls-Instabilitaet in hexagonalen Gittern, Kekulé-O-Type-Symmetry.
- Reentrant localization (Guan et al. 2023, Phys. Rev. A), non-Hermitian topological phases (Liu-Chen, arXiv:2405.11812, “Lindbladian dynamics with loss of quantum jumps”).

**PHY016-Korrektur**: Die Behauptung “Wallstrom auf N=6 trivial” ist **falsch**. Korrekte Formulierung: Im freien Bose-Hubbard-Bereich mit reeller Hopping-Amplitude (φ=0) und nicht-entartetem Grundzustand mit nodenfreier Wellenfunktion ist die topologische Klasse trivial (k=0). Sobald φ ≠ 0 oder das System angeregt wird, sind nicht-triviale Winding-Klassen besetzt; das Spektrum E_k = −2t·cos[2π(k+φ/φ_0)/6] beweist sechs distinkte Bloch-Klassen.

### Achse 4 — Hydraulik, Wasser, Kapillaritaet, klassische Hydrodynamik

**Aggregatzustaende Wasser (Stand 2025)**:

- Stand: **20 kristalline Eis-Phasen** identifiziert (Lei-Liu-Yu-Niu, The Innovation 6 (5), 100881, “Deep potential-driven structure exploration of ice polymorphs”, March 2025;  arXiv:2504-x). Salzmann-Review (arXiv:1812.04333) zaehlte 17 in 2019; juengste Entdeckungen erweiterten den Katalog.
- Ice XIX: T.M. Gasser, A.V. Thoeny, A.D. Fortes, T. Loerting, “Structural characterization of ice XIX as the second polymorph related to ice VI”, Nat. Commun. 12, 1128 (18 Februar 2021); DOI: 10.1038/s41467-021-21161-z (Univ. Innsbruck).  Unabhaengig bestaetigt durch R. Yamane et al. (Univ. Tokyo), “Experimental evidence for the existence of a second partially-ordered phase of ice VI”, Nat. Commun. 12, 1129 (2021), DOI: 10.1038/s41467-021-21351-9. 
- Theoretisch vorhergesagt: monolayer ice, zeolite-like ice polymorphs, ice 0, ice χ, ice i (Kapil-Schran-Zen et al., Nature 609, 512-516, 2022). 
- Amorphe Phasen: LDA, HDA, VHDA; supercritical water > 374 °C, 22.1 MPa; plasma water bei extremen Bedingungen.
- Plastic ice phase recently characterized in melt-growth (Nature Comm. Phys., 2023).

**Erdwaerme, 500m-Druckroehren**:

- Deep Borehole Heat Exchanger (DBHE) > 500 m Tiefe sind etabliert (Brown et al., Geothermal Energy 12, “A comprehensive review of deep borehole heat exchangers (DBHEs)”, 2024).
- Realistische Tiefe: 800-2000 m (Lund, Deep Underground Sci. Eng. 2, 2024, “Performance analysis of deep borehole heat exchangers”);  ueber 500 m mit U-Tube-Konfiguration unueblich, ueblicher coaxial.
- COP von GSHP-Systemen: 3-5 (Zhang et al. 2022). 
- Pressure drop und thermisches Short-Circuiting limitieren Effizienz. 
- New York State Regulierung (Hochul, S.8060/A.8565, Feb 2024) speziell fuer Borholes > 500 ft. 

**Kapillaritaet in Baeumen, Cohaesion-Tension-Theorie**:

- Dixon-Joly 1894 (klassische C-T-Theorie); Dixon 1914 hat das Buch “Transpiration and the Ascent of Sap in Plants” publiziert.
- Koch et al. 2004 (Nature 428, 851-854, “The limits to tree height”): theoretisches Maximum 122-130 m, basierend auf Leaf-Funktional-Gradienten in Sequoia sempervirens.  Damals gemessene Hyperion-Hoehe 112.7 m;  **aktualisierter Wert 116.07 m gemaess Guinness World Records 2019;  Wikipedia-Eintrag fuehrt 116.22 m fuer 2026 an**. 
- Hydrostatischer Gradient: −0.0096 ± 0.0007 MPa/m (Koch et al. 2004, gemessen), praktisch identisch mit dem theoretischen Gravitations-Gradient −0.0098 MPa/m. 
- Gouin 2008/2012 (arXiv:0809.3529, arXiv:1204.4094, “Nanofluidics explanation of ascent of water in tallest trees”): Disjoining-pressure-Modell erlaubt sap ascent ueber 100 m via van-der-Waals-Korrekturen jenseits klassischer Kapillaritaet. 
- Negativ-Druck-Toleranz in Xylem: typisch −2 bis −5 MPa, metastabile Zustaende; Preston-Experiment 1952 zeigte allerdings, dass Baeume mit ueberlappenden Saw-Cuts ueberleben  — strittig.

**Hagen-Poiseuille vs Murray’s Law vs Constructal Law**:

- Hagen-Poiseuille: ΔP = 8μLQ/(πr⁴) — laminare Roehrenstroemung.
- Murray’s Law: r₀³ = Σ r_i³ fuer Bifurkationen — minimiert Pumparbeit + metabolische Kosten.
- Constructal Law (Bejan 1996): Flussstrukturen evolvieren in Zeit zu erleichtertem Zugriff fuer Stroemung. “For a finite size flow system to persist in time, it must evolve with freedom such that it provides easier and greater access to what flows.”
- Bejan et al., Physics of Life Reviews 50, 103-116 (Sep 2024), “Evolution and irreversibility: Two distinct phenomena and their distinct laws of nature”: explizit “Evolution und Irreversibilitaet sind zwei distinkte Phaenomene mit distinkten Gesetzen — Constructal Law und Second Law sind verschieden”. 
- **Kritik**: Constructal Law bleibt ausserhalb des Mainstream der theoretischen Physik (Mullaly, Medium-Essay “A Critique of Constructal Theory”, 2021);  quantitative Vorhersagen ueberlappen oft mit Murray/Bejan-Lorente-Optimierung.

**Wirbelphysik**:

- Quanten-Wirbel in Superfluid He-II: zirkulationsquantisiert ∮v·dl = nh/m_4.
- BEC-Wirbel: K.W. Madison, F. Chevy, W. Wohlleben & J. Dalibard, “Vortex formation in a stirred Bose-Einstein condensate”, Phys. Rev. Lett. 84, 806-809 (31 Januar 2000), Laboratoire Kastler Brossel, ENS Paris; erster Bericht von bis zu vier simultanen Wirbeln in einem geruehrten ⁸⁷Rb-Kondensat. 
- Klassische Wirbel: nicht quantisiert; Smoke-Rings, Tornados.

**Doppelspalt: Wasser vs Licht**:

- Wasserwellen zeigen Interferenz; aber kein Komplementaritaets-Prinzip (kein “particle-like” Modus).
- Light/Quanten: Komplementaritaet via V²+D²≤1, etabliert durch B.-G. Englert, “Fringe Visibility and Which-Way Information: An Inequality”, Phys. Rev. Lett. 77, 2154-2157 (1996); die Ungleichung steht auf p. 2154: “this inequality can be regarded as quantifying the notion of wave-particle duality”. 
- Triality: P²+V²+E²=1 fuer pure states im n-path-Interferometer mit which-path-Detector (Qureshi, Phys. Rev. A 100, 042105, 2019; Roy-Pathania-Chandra-Panigrahi-Qureshi, Phys. Rev. A 105, 032209, 2022; Tsui-Kim, Phys. Rev. A 109, 052439, 2024; Yang-Wang-Fei, Phys. Rev. A 110, 042413, 2024). 
- Light Sci. Appl. 2025 (“Universal conservation laws of the wave-particle-entanglement triad: theory and experiment”):  experimentelle Verifikation.

**Verbindung Quanten-/Klassische Hydrodynamik**:

- Madelung-Form: ∂_t ρ + ∇·(ρv) = 0; m∂_t v + m(v·∇)v = −∇(V + Q).
- Klassisch entspricht V Newton-Potential, Q ≡ 0; Madelung Q liefert “quantum pressure”.
- Tsekov-Heifetz-Cohen, EPL 122, 40002, 2018: Madelung ~ Hydrodynamische Turbulenz-Analogie.
- Heifetz-Cohen, Found. Phys. 45, 1514, 2015: Fisher-Information als “internal thermal energy” der Madelung-Fluessigkeit, mit Caveat dass Inkompressibilitaet die Erhaltung von kinetischer Energie und Fisher-Information separiert.
- Grenze: in der Klassischen Hydrodynamik fehlt die Zirkulations-Quantisierung; in Quanten-Madelung ist Multi-Valued-ness der Phase zentral (Wallstrom).
- **Aktuelle 2024 Forschung**: hydrodynamische Phasenuebergaenge in cold-atom-Systemen, BEC-Solitons, Quantum-Hydrodynamik-Limits.

### Breite Themen

**Komplementaritaet & Triality**:

- Englert-Bergou EGY-Relation: V² + D² ≤ 1 (Englert PRL 77, 2154, 1996).
- Qureshi 2019/2020 (Phys. Rev. A 100, 042105; arXiv:2011.08210): D² = P² + E² (Pythagoras: Distinguishability, Predictability, Entanglement). 
- Triality: P² + V² + E² ≤ 1 (Roy et al., Phys. Rev. A 105, 032209, 2022; Tsui-Kim, Phys. Rev. A 109, 052439, 2024 fuer n-path). 
- Light Sci. Appl. 2025: experimentell verifiziert mit entangled photons in dual-path interferometer.

**Pfeil der Zeit**:

- BaBar 2012: Lees et al. (BABAR Collaboration), Phys. Rev. Lett. 109, 211801 (2012); “time-reversal violation is clearly established, with the exclusion of the (0,0) point with a significance of 14σ”  (bestaetigt in arXiv:1307.2759, DIS2013 proceedings). Direkter T-Verletzungsnachweis in B^0-Mesonen via Austausch |in⟩↔|out⟩.  
- Past Hypothesis (Albert, Penrose): low-entropy Initial-State des Universums.
- Al-Khalili, Chen, “The Decoherent Arrow of Time and the Entanglement Past Hypothesis”, Found. Phys. 54, 49 (2024); DOI: 10.1007/s10701-024-00785-3: fuer den decoherenten Pfeil der Zeit braucht es zusaetzlich zu thermodynamischer PH die Annahme **niedriger Entanglement-Entropie im Anfangszustand des Universums**. 

**Backreaction Kosmologie**:

- Buchert-Formalismus: Q_D Backreaction-Term in gemittelten Friedmann-Gleichungen. 
- Buchert-Raesaenen 2012 (“Backreaction in late-time cosmology”), Yao-Meng 2017/2024 (arXiv:2406.15442 “CPL effective dark energy from the backreaction effect”): CPL-parametrisierte effektive Dark Energy aus Backreaction. 
- **Strittig**: Mainstream-Konsens, dass Backreaction zu klein fuer volle Dark-Energy-Erklaerung; Minderheit (Buchert, Wiltshire, Raesaenen) argumentiert weiterhin signifikant. **Beide Seiten darstellen**.

**Lindblad-Master-Gleichung State-of-the-Art 2024-2026**:

- McCauley et al. 2019: universally-valid completely-positive Form fuer weakly-damped Systeme. 
- Trushechkin 2021: non-secular Lindblad mit Erhaltung von Coherence-Transfer in nahezu degenerierten Systemen. 
- Liu-Chen 2024 (arXiv:2405.11812 “Lindbladian dynamics with loss of quantum jumps”): nichtlineare Lindblad-Master-Gleichung fuer postselected dynamics, postselected skin effect. 
- Borras et al. 2024: Quantum-Algorithm fuer Lindblad-Simulation mit diamond-norm Fehler O(δt³) per Timestep. 
- Non-Markovian Erweiterungen: transfer-tensor method, path-integral coupling (Bose 2024).
- Kraft et al. 2024: Mapping von open-system NESS auf klassische Korrelationsfunktionen in weak-coupling/driving Regime. 
- Jung et al. 2025 (arXiv:2505.x): Quantum Optical Master Equation (QOME) als faithfulster Reproduzent von Redfield-Benchmark wenn Eigenbasis bekannt. 

**Detailed Balance & Gibbs Steady States**:

- Detailed Balance gibt thermodynamische Konsistenz: Lindblad-Operatoren in pairs L, L† mit Boltzmann-Faktor; Gibbs steady state ρ_∞ = exp(−βH)/Z.
- Bei kollektiver Dissipation (Dicke-Modelle, sub-radiance) bricht single-site detailed balance, aber Cluster-DB kann gelten.

**Aubry-Andre**:

- Klassisches AA-Modell (1980): localization transition bei V_c = 2t. 
- Mobility edges in generalisierten AA-Modellen (Liu-Guo-Pu-Longhi, arXiv:2007.06259;  Zhou-Wang-Liu, “Exact New Mobility Edges between Critical and Localized States”, PRL 131, 176401, 2023). 
- Reentrant localization in spinful AA mit non-Abelian gauge (Guan-Wang-Guan-Cai 2023).
- Mobility rings in non-Hermitian non-Abelian quasiperiodic lattices (arXiv:2507.12176, 2025). 
- Auf einer 6-Ring-Topologie: Standard-AA-Loc nicht erreichbar, da Anzahl Sites zu klein fuer kritisches V_c.

**Finite-Size-Scaling N=6**:

- N=6 ist zu klein fuer Standard-FSS (Binder cumulants etc.).
- Empfohlen: Sequenz N ∈ {6, 12, 18, 24, …} und Extrapolation 1/N → 0.
- Renormalization-group via successive cluster doubling.

-----

## Details — Korrekturen fuer PHY010-016

### PHY010 (Mean-Field-Bifurkation)

**Problem**: Single-Site-Gutzwiller-MF zeigt scharfe Phase-Boundary, wahrscheinlich Artefakt.
**Korrektur**:

1. Hilbert-Raum-Dimension berechnen: fuer N=6 Sites mit n_max = 6 Boson pro Site, Total-Number-Konservierung N_b = 6 ⇒ Dim ≈ 462. Trivial in ED.
1. Cluster Gutzwiller mit ganzem Ring als Cluster (Luehmann, Phys. Rev. A 87, 043619, 2013).
1. Lanczos-Algorithm fuer Grundzustand und niedrige Anregungen.
1. Berechne Compressibility κ = ∂n/∂μ und Superfluid Stiffness ρ_s als Crossover-Indikatoren statt scharfer Bifurkation.

### PHY011-013 (Peierls-Phase, Persistent Current)

**Problem**: Vermutete Trivialitaet der topologischen Klassen.
**Korrektur**:

1. Spektrum E_k = −2t·cos[2π(k+φ/φ_0)/6] explizit ausarbeiten fuer alle k=0..5.
1. Persistent current I(φ) = −∂E/∂φ als Funktion von φ ploten — Periode φ_0.
1. Identifiziere Level-Crossings bei φ = (1/2)φ_0 (Aharonov-Casher-aehnlich).

### PHY014 (Lindblad-Dissipation)

**Problem**: vermutlich Standard-Lindblad ohne detailed-balance-Check.
**Korrektur**:

1. Pruefe Gibbs-Steady-State: Loese L[ρ_∞] = 0 numerisch und vergleiche mit exp(−βH)/Z.
1. Pruefe positivity-preservation: Eigenwerte ρ(t) ≥ 0.
1. Erwaege coarse-grained ME oder QOME (Jung et al. 2025) fuer hoehere Genauigkeit.

### PHY015 (Madelung-Identitaet)

**Problem**: 67%-Verletzung der Bohm-Fisher-Identitaet.
**Korrektur**:

1. Beide Diskretisierungsformen explizit dokumentieren:
- Q_i = −(ℏ²/2m)·(√ρ_{i+1} + √ρ_{i−1} − 2√ρ_i)/(a²·√ρ_i) (Drei-Punkt-Laplace, semiklassisch)
- I_F^{discrete} = (4/a²)·Σ(√ρ_{i+1} − √ρ_i)² (Hellinger-Link-Form, Sason et al. 2019)
1. Auf periodischem Ring mit glatter Grundzustands-Dichte rechnen.
1. Skalierung N=6, 12, 24, 48, 96 — erwarte O(N⁻²)-Konvergenz fuer Drei-Punkt-Form, exponentiell fuer Link-Hellinger-Form (Trefethen-Weideman, SIAM Rev. 56, 385, 2014).

### PHY016 (Wallstrom-Behauptung)

**Problem**: “Wallstrom auf N=6 trivial”.
**Korrektur**: ersetzen durch: “Im freien Bose-Hubbard-Bereich mit nodenfreier Grundzustands-Dichte und Peierls-Phase φ=0 ist nur die k=0 topologische Klasse besetzt; bei φ ≠ 0 verschiebt sich der besetzte k-Wert; bei angeregten Zustaenden oder Multi-Mode-Superposition sind alle 6 Klassen prinzipiell erreichbar.” Referenz: Maiti arXiv:0706.0061; Heinrichs arXiv:cond-mat/0106437; Byers-Yang 1961.

-----

## Recommendations

**Sofort (Sprint PHY017)**:

1. PHY010 mit ED rechnen und Crossover-Charakter dokumentieren (Hilbert-Dimension 462, Lanczos).
1. PHY015 als Skalierungs-Studie N ∈ {6, 12, 24, 48, 96} aufsetzen mit beiden Operator-Familien.
1. PHY016 textuell korrigieren mit Verweis auf Maiti-Heinrichs Tight-Binding-Ring-Spektrum.

**Mittelfristig (PHY018-020)**:
4. Cluster-Gutzwiller-Implementation fuer realistische Phasen-Diagramme (Luehmann-Methode).
5. DMRG-Implementation fuer N ≥ 24 (ITensor, TenPy).
6. Lindblad-Gibbs-Steady-State-Verifikation mit detailed-balance-Test; Wechsel auf QOME oder coarse-grained ME erwaegen.

**Lang (PHY021+, eigenstaendige publizierbare Beitraege)**:
7. **Erstes Forschungsthema**: Rigorose Theorie der diskreten Madelung-Bohm-Fisher-Identitaet auf C_N-Ringen (publizierbarer Beitrag, da klare Luecke in der Literatur — Subagent-Recherche bestaetigt: kein paper deckt das explizit ab).
8. **Zweites Thema**: PMPG-Verallgemeinerung auf inkompressible Quantenfluide (existiert noch nicht in publizierter Form — Taha-Methodik plus Madelung-Reformulierung).
9. **Drittes Thema**: Cluster-Mean-Field-Validitaetskriterien fuer Ring-Geometrien — quantitative Karte z vs. Quantenfluktuation.

**Benchmarks fuer Re-Evaluation**:

- Wenn PHY015-Konvergenz nicht O(N⁻²) oder exponentiell ist → fundamentaler Diskretisierungs-Fehler im Code.
- Wenn PHY010-Crossover nach ED weiterhin scharfe Sigmoid-Struktur zeigt → echter finite-size Phasenuebergang, nicht MF-Artefakt.
- Wenn Persistent Current I(φ) keine Periode φ_0 zeigt → Numerischer Fehler in Peierls-Implementation.

**Budget-Schaetzung CHF (grobe Naeherung)**:

- ED-Implementation in Julia/Python: ~20 Stunden Coworker-Zeit, ~CHF 2’500 bei 125 CHF/h.
- Cluster-Gutzwiller-Implementation: ~40 Stunden, ~CHF 5’000.
- DMRG via ITensor: ~30 Stunden setup + Lernkurve, ~CHF 3’750.

-----

## Caveats

- **PMPG fuer Quanten**: ein eigenstaendiges Quanten-PMPG mit publizierter Variations-Aequivalenz zu Schrödinger existiert nach derzeitiger Kenntnis **nicht**; saemtliche Aussagen dazu sind heuristisch.
- **Mauri-Giona 2025**: konnte in der Suche nicht direkt verifiziert werden — UNBEKANNT; Sprint-Material sollte direkt referenziert sein.
- **Wallstrom-Phaenomen auf Graphen**: keine rigorose Theorie publiziert (Reddiger-Poirier 2023 limitiert auf Kontinuum; Subagent-Recherche bestaetigt) — Phi-Hex-Beitragspotenzial.
- **Constructal Law als Naturgesetz**: bleibt ausserhalb des Mainstreams theoretischer Physik (Mullaly-Kritik 2021); quantitative Vorhersagen oft mit Murray/Bejan-Optimierung redundant.
- **Backreaction-Erklaerung der Dark Energy**: Mainstream lehnt ab (Li-Schwarz 2007 zeigt: Effekt zu klein in linearer Stoerungstheorie), Minderheit (Buchert, Raesaenen, Wiltshire) verteidigt — strittig.
- **Bohmsche Mechanik / Madelung / MWI**: wissenschaftlich neutral darstellen; Madelung-Hydrodynamik als Reformulierung, nicht als alternative Theorie.
- **N=6 als Demonstrator-System**: Phi-Hex ist methodisch wertvoll, aber zu klein fuer rigorose Phasenuebergangs-Analyse; ist eher Spielzeug-Modell als kondensiertes-Materie-Modell.
- **Tree-Height-Limit 130 m**: theoretischer Wert (Koch et al. 2004); gemessenes Maximum aktuell Hyperion mit 116.07 m (Guinness 2019),  Wikipedia fuehrt 116.22 m fuer 2026  — Sequoia-135-m-Angabe in aelterer Literatur (Flindt, “Amazing numbers in biology”) unbestaetigt.
- **20 Eis-Phasen vs 17**: aktuelle Innovation-Paper (Lei et al. 2025) zaehlt 20 kristalline Phasen; Salzmann-Review (2019) zaehlte 17; Diskrepanz durch neu bestaetigte Ice XVIII, XIX, XX zwischen 2021-2025.
- **Konvergenzordnung der diskreten Madelung-Identitaet**: keine geschlossene Theorie publiziert; **die hier praesentierten O(N⁻²)- und exp(−cN)-Aussagen folgen aus numerischen Analyse-Standard-Resultaten (Trefethen-Weideman, Iannaccone-Curatola-Fiori), nicht aus einem Bohm-Fisher-spezifischen Theorem**. Heuristisch markiert.