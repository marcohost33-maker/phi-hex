-----

## id: 260602_PHI_HEX_audit_4responses_v01
title: Audit der letzten vier Reaktionen - Korrektur der Extrapolation und verifizierte T_BKT-Bestimmung
attribution: Coworker Research / Coworkerz
date: 2026-06-02
status: R2_REFERENCE_EVIDENCE
tier: VERIFIZIERT
phase: Methodische Selbstkorrektur und unabhaengige Validierung
warning_level: GRUEN

# Audit der letzten vier Reaktionen

Coworker Research / Coworkerz, 02. Juni 2026

## Worum es geht

Dieses Dokument prueft die letzten vier Reaktionen kritisch und haelt fest, welcher Fehler dabei auftrat, wie er entdeckt wurde, und mit welcher ueberpruefbaren Vorkehrung er behoben wurde. Die vier Reaktionen umfassten zwei Meta-Audits des Gesamtchats, die Code-Konsolidierung mit der ersten Finite-Size-Scaling-Bestimmung, sowie die Wolff-Implementierung mit der darauf aufbauenden Praezisionsbestimmung. Der wichtigste Befund betrifft die vorletzte Reaktion und ist gravierend, doch seine Korrektur hat das wissenschaftliche Ergebnis des gesamten Projekts deutlich verbessert.

## Befund 1: Scheinpraezision in der Extrapolation (gravierend, behoben)

In der Reaktion zur Praezisionsbestimmung PHY027 wurde die kritische Temperatur durch eine Extrapolation der naiven Helicity-Crossings gegen die Form eins durch Logarithmus L zum Quadrat bestimmt, mit dem Ergebnis 1,481. Dieses Ergebnis wurde so dargestellt, als folge es zwingend aus den Daten. Die Nachpruefung in dieser Runde zeigt, dass das nicht stimmt. Drei gleichermassen plausible Extrapolationsformen, naemlich eins durch Logarithmus L, eins durch Logarithmus L zum Quadrat, und eins durch L, liefern Grenzwerte von 1,420, 1,482 und 1,496. Mit nur vier Datenpunkten und je zwei Anpassungsparametern sind diese Formen statistisch nicht zu unterscheiden, da alle aehnlich gut passen. Der praesentierte Wert war damit im Wesentlichen eine Funktion der willkuerlichen Modellwahl und nicht der Daten. Es wurden zudem weder Fehlerbalken auf den extrapolierten Wert gerechnet noch die Modellabhaengigkeit offengelegt. Methodisch ist das dieselbe Klasse von Fehler wie die Bestaetigungsverzerrung, die in einer frueheren Projektphase bereits als Negativ-Befund dokumentiert worden war, nur diesmal in der Form einer vorgetaeuschten Genauigkeit.

## Befund 2: Die korrekte Methode war bekannt, aber nicht umgesetzt

Die Web-Recherche bestaetigte, dass die etablierte Loesung dieses Problems nicht eine frei gewaehlte Extrapolation ist, sondern die Sandvik-Paar-Methode. Diese war im eigenen Meta-Audit zuvor sogar korrekt beschrieben, dann aber in der Implementierung nicht verwendet worden. Die Methode beruht auf einer eleganten Idee. Am kritischen Punkt erfuellt der Helicity-Modulus die Weber-Minnhagen-Beziehung, die eine unbekannte Konstante enthaelt. Statt diese Konstante zu schaetzen, betrachtet man ein Paar von Gittergroessen, bei dem die eine doppelt so gross ist wie die andere. Bildet man die Differenz der inversen Korrekturterme beider Groessen, so faellt die unbekannte Konstante exakt heraus, und es bleibt ein fester, theoretisch bekannter Wert, naemlich minus zwei mal der Logarithmus des Groessenverhaeltnisses. Die kritische Temperatur ist dann jene Temperatur, bei der diese konstantenfreie Bedingung erfuellt ist. Es muss keine Extrapolationsform mehr geraten werden.

## Befund 3: Unabhaengige Verifikation auf dem Quadratgitter

Bevor die korrigierte Methode auf das eigentliche Problem angewandt wurde, wurde sie auf dem Quadratgitter verifiziert, dessen kritische Temperatur mit hoher Praezision bekannt ist und 0,8935 betraegt. Das ist eine Code-Verifikation im Sinne der Verifikations- und Validierungs-Methodik, denn sie prueft die Korrektheit des Verfahrens an einem Fall mit bekannter Antwort. Das Ergebnis ist ueberzeugend: Das groesste Gittergroessen-Paar reproduziert den bekannten Wert mit 0,8950, einer Abweichung von nur 0,16 Prozent. Damit ist zweierlei unabhaengig bewiesen. Erstens arbeitet der Wolff-Cluster-Sampler korrekt. Zweitens funktioniert die konstantenfreie Paar-Methode. Das kleinere Paar zeigte erwartungsgemaess keinen Nulldurchgang, was die Literaturaussage bestaetigt, dass zu kleine Gitter den asymptotischen Bereich noch nicht erreichen.

## Befund 4: Die korrigierte Bestimmung auf dem Triangular-Gitter

Auf das eigentliche Triangular-Gitter angewandt, liefert die verifizierte Methode eine deutlich bessere Bestimmung als zuvor. Das kleinere Paar ergibt 1,425, das groessere 1,409. Beide Paare klammern den Literaturwert von 1,418 ein und liegen auf unter einem Prozent. Die Bestschaetzung als Mittel der beiden Paare betraegt 1,417 mit einer Spanne von plus minus 0,008, in praktisch exakter Uebereinstimmung mit dem Referenzwert von 1,418 mit der Literatur-Unsicherheit von zwei in der letzten Stelle. Das ersetzt die fragwuerdige fruehere Angabe von 1,481 vollstaendig und beruht erstmals auf einer sauberen, unabhaengig verifizierten Methode mit ausgewiesener Unsicherheit.

## Befund 5: Eine vorschnelle Fehlerdiagnose im eigenen Pruefprozess

Ein lehrreicher Nebenbefund betrifft den Pruefprozess selbst. Im Verlauf der Diagnose erschien der Helicity-Modulus bei tiefer Temperatur mit einem Wert von etwa 0,7 zunaechst als Implementierungsfehler, da der nackte Kopplungswert eins erwartet wurde. Die genauere Pruefung zeigte jedoch, dass die Formel auf einem perfekt ausgerichteten Zustand exakt eins liefert, der niedrigere Wert bei endlicher Temperatur also korrekt ist. Es handelt sich um die physikalisch reale Renormierung der Steifigkeit durch Spin-Wellen, die auf endlichen Gittern den gemessenen Wert unter den nackten Kopplungswert drueckt. Das urspruenglich formulierte Validierungs-Gate, das exakt den Kopplungswert forderte, war damit physikalisch naiv und wurde durch ein Gate ersetzt, das nur den zielwert-unabhaengig korrekten Fall des perfekt ausgerichteten Zustands prueft. Dieser Nebenbefund ist wichtig, weil er zeigt, dass die Disziplin des kritischen Nachpruefens in beide Richtungen wirken muss. Sie muss nicht nur ein vorschnelles Erfolgsurteil abfangen, sondern auch eine vorschnelle Fehlerdiagnose, die einen korrekten Befund faelschlich als Fehler verwirft.

## Staerken der geprueften Reaktionen, die bestehen bleiben

Mehrere Kernergebnisse der vier Reaktionen sind durch die Pruefung bestaetigt worden und bleiben gueltig. Die Wolff-Implementierung selbst ist korrekt, wie die Quadratgitter-Verifikation unabhaengig belegt. Der gemessene Geschwindigkeitsgewinn gegenueber dem lokalen Update ist real und betraegt im Mittel etwa einen Faktor zwoelf. Die qualitative Aussage, dass das XY-gekoppelte Modell in der korrekten Universalitaetsklasse liegt, ist nicht nur bestaetigt, sondern jetzt quantitativ auf unter ein Prozent abgesichert. Die beiden Meta-Audits waren inhaltlich solide und haben die Roadmap korrekt vorgezeichnet; der Fehler lag nicht in der Planung, sondern in der Umsetzung der Extrapolation.

## Die methodische Lehre als ueberpruefbare Vorkehrung

Der wertvollste Teil dieser Pruefung ist die dauerhafte Vorkehrung gegen eine Wiederholung. Aus dem Hauptbefund folgt eine konkrete Regel: Wann immer ein numerisches Ergebnis durch eine Extrapolation in einen Grenzfall gewonnen wird, ist die Wahl der Extrapolationsform zu begruenden und ihre Auswirkung auf das Ergebnis offenzulegen. Wenn mehrere plausible Formen stark verschiedene Grenzwerte liefern, ist das Ergebnis als unsicher auszuweisen, und es ist, wo immer moeglich, eine Methode zu verwenden, die ohne freie Formwahl auskommt, wie die konstantenfreie Paar-Methode. Diese Regel ist als zielwert-unabhaengiges Kriterium formuliert und damit selbst gegen Bestaetigungsverzerrung geschuetzt. Sie ergaenzt die bereits etablierte Regel, dass Normierungen vorwaerts aus der Theorie abzuleiten und vor dem Abgleich mit dem Zielwert einzufrieren sind.

## Ablage-Empfehlung

Abzulegen sind das Quadratgitter-Validierungsmodul PHY028, das korrigierte Triangular-Bestimmungsmodul PHY029, ihre beiden Reports sowie dieses Audit-Dokument. Ein Architecture Decision Record sollte die Ruecknahme der freien Extrapolation aus PHY027 dokumentieren und die konstantenfreie Sandvik-Paar-Methode als neuen Standard fuer alle Grenzwert-Bestimmungen festschreiben. Ergaenzend sollte im Negativ-Ergebnis-Register vermerkt werden, dass die freie Extrapolationsform ein vermeidbarer Fehler war, der durch die unabhaengige Quadratgitter-Verifikation aufgedeckt wurde, und dass das naive Tieftemperatur-Gate als physikalisch unzutreffend verworfen wurde.

ENDE des Audit-Dokuments. Status v01, GRUEN klassifiziert. Die kritische Temperatur des Triangular-Gitters ist mit einer verifizierten Methode auf unter ein Prozent bestimmt und stimmt mit dem Literaturwert ueberein.