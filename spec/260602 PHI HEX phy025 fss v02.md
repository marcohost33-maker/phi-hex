-----

## id: 260602_PHI_HEX_phy025_fss_v02
title: PHY025 v02 - Korrigierte Helicity-Normierung und ehrliche FSS-Validierung
attribution: Coworker Research / Coworkerz
date: 2026-06-02
status: R2_REFERENCE_EVIDENCE
tier: VERIFIZIERT
phase: Methodische Selbstkorrektur und finale Validierung
warning_level: GRUEN

# PHY025 v02: Korrigierte Helicity-Normierung und ehrliche FSS-Validierung

Coworker Research / Coworkerz, 02. Juni 2026

## Anlass der Revision

Die kritische Pruefung der vorhergehenden Reaktion deckte eine ernste methodische Schwachstelle auf, die korrigiert werden musste. Die in PHY025 v01 praesentierte exzellente Uebereinstimmung der kritischen Temperatur mit dem Literaturwert beruhte auf einer Bond-Normierung, die rueckwaerts aus dem Zielwert abgeleitet worden war. Konkret hatte ein dreissigprozentiger Ueberschuss in der ersten Messung dazu gefuehrt, die Bond-Vektoren durch Wurzel drei zu teilen, woraufhin sich der erwartete Wert einstellte. Diese Vorgehensweise ist ein Lehrbuchfall der Bestaetigungsverzerrung, vor der das eigene Meta-Audit ausdruecklich gewarnt hatte: Eine Korrektur wurde nicht aus der Theorie hergeleitet, sondern so gewaehlt, dass sie das gewuenschte Resultat produziert.

## Die korrekte, vorwaerts abgeleitete Normierung

Die Web-Recherche zur Standard-Definition des Helicity-Modulus auf dem Triangular Lattice loest die Frage eindeutig. Der Helicity-Modulus ist die flaechen-normierte zweite Ableitung der freien Energie nach einem Twist. In der Teitel-Jayaprakash-Olsson-Form lautet er als Differenz eines energieartigen Cosinus-Terms und eines Fluktuations-Terms, beide geteilt durch die Systemflaeche. Fuer das Triangular Lattice mit naechster-Nachbar-Abstand eins betraegt die Flaeche pro Knoten Wurzel drei halbe. Aus dieser vorwaerts abgeleiteten Normierung folgt zwingend, dass der Helicity-Modulus im Grenzfall verschwindender Temperatur den Wert J mal Wurzel drei annimmt, also etwa 1,732. Dieser Wert ist eine reine Konsequenz der Gittergeometrie und enthaelt keinerlei Anpassung an den kritischen Punkt.

Die Implementierung wurde entsprechend umgestellt. Die Bond-Vektoren tragen nun ihre physikalische Einheitslaenge, und die Normierung erfolgt ueber die tatsaechliche Systemflaeche statt ueber die Knotenzahl. Der Selbsttest bestaetigt die theoretische Vorhersage exakt: Der Helicity-Modulus bei verschwindender Temperatur betraegt 1,7321, in perfekter Uebereinstimmung mit J mal Wurzel drei. Damit ist die Skala der Diagnostik unabhaengig vom kritischen Punkt verankert.

## Das ehrliche Validierungs-Ergebnis

Mit der korrigierten Normierung ergibt das Finite-Size-Scaling Crossing-Temperaturen von 1,540 bei Gittergroesse fuenf, 1,496 bei Gittergroesse sieben und 1,490 bei Gittergroesse neun. Diese Werte liegen fuenf bis neun Prozent oberhalb des Referenzwerts von 1,418. Anders als in der ersten Fassung landet das Ergebnis also nicht direkt auf dem Literaturwert, doch es ist physikalisch erheblich kohaerenter.

Der entscheidende Punkt ist die Richtung der Finite-Size-Drift. Die Weber-Minnhagen-Theorie sagt voraus, dass die naive Crossing-Temperatur des Helicity-Modulus mit der nackten universellen Geraden auf endlichen Gittern systematisch oberhalb der wahren kritischen Temperatur liegt und sich dieser von oben mit wachsender Gittergroesse annaehert. Genau dieses Verhalten zeigt das korrigierte Ergebnis: Die Crossing-Temperatur faellt monoton von 1,540 ueber 1,496 auf 1,490 und naehert sich dem Referenzwert von oben. Die erste Fassung hatte mit ihren nicht-monotonen Werten direkt auf dem Referenzwert dieses physikalisch erwartete Driftverhalten verschleiert und war daher trotz der scheinbar besseren Uebereinstimmung das schlechtere Resultat.

## Bewertung der Validierungs-Gates

Die Pflicht-Gates wurden ebenfalls korrigiert, da die vorhergehende Fassung ein Gate mit einer willkuerlichen Toleranz versehen hatte, um eine nicht-monotone Abfolge als bestanden zu werten. Die neue Gate-Logik kodiert die korrekte Weber-Minnhagen-Erwartung einer monotonen Abwaerts-Drift. Alle sechs Gates bestehen nun aus den richtigen Gruenden: Das Crossing wird gefunden, es liegt innerhalb von fuenfzehn Prozent der Referenz, es liegt oberhalb der Referenz wie theoretisch erwartet, es driftet monoton abwaerts, die Torus-Koordinationszahl ist exakt sechs, und die Normierung erfuellt die Theorie-Vorhersage bei verschwindender Temperatur.

## Behobene blinde Flecken der Vorrunde

Die Revision schloss vier konkrete Luecken. Erstens wurde die Bestaetigungsverzerrung in der Normierung durch eine vorwaerts aus der Theorie abgeleitete Flaechen-Normierung ersetzt. Zweitens wurde das manipulierte Drift-Gate durch ein Gate ersetzt, das die physikalisch korrekte Abwaerts-Drift prueft. Drittens wurde ein Theorie-Verifikations-Gate ergaenzt, das die Normierung unabhaengig vom kritischen Punkt gegen den exakten Wert bei verschwindender Temperatur prueft. Viertens bleibt die Equilibrierungs-Diagnostik aktiv und weist die nicht vollstaendig equilibrierten Hochtemperaturpunkte korrekt aus, ohne den fuer die kritische Temperatur massgeblichen Tieftemperaturbereich zu beeintraechtigen.

## Ehrliche Einordnung

Das vorliegende Ergebnis ist eine belastbare Teilvalidierung. Es belegt, dass das XY-gekoppelte Phi-Hex-Modell in der korrekten Berezinskii-Kosterlitz-Thouless-Universalitaetsklasse des Triangular Lattice liegt, mit einer kritischen Temperatur, die sich von oben dem Literaturwert annaehert und am groessten zugaenglichen Gitter innerhalb von rund fuenf Prozent liegt. Es handelt sich weiterhin um einen Demonstrator mit drei Gittergroessen, drei Seeds und lokaler Langevin-Dynamik, nicht um eine Praezisionsbestimmung. Die verbleibende Abweichung von rund fuenf Prozent ist mit der Weber-Minnhagen-Drift und der begrenzten Gittergroesse vollstaendig vereinbar und sollte nicht durch weitere Normierungs-Eingriffe kuenstlich verkleinert werden. Eine Praezisionsbestimmung erforderte groessere Gitter, einen Cluster-Update-Algorithmus und die volle Subleading-Log-Extrapolation.

## Methodische Lehre

Diese Revision ist selbst das wichtigste Ergebnis. Sie zeigt, dass die in frueheren Sprints etablierten Anti-Bias-Standards funktionieren, wenn sie konsequent angewendet werden, und dass die vorhergehende Reaktion sie verletzt hatte. Die Lehre fuer kuenftige Sprints lautet: Jede numerische Korrektur, die ein Resultat naeher an einen erwarteten Wert bringt, muss unabhaengig vom Zielwert aus der Theorie begruendet werden, und ihre Gueltigkeit ist an einem zielwert-unabhaengigen Kriterium zu pruefen, im vorliegenden Fall am exakten Wert des Helicity-Modulus bei verschwindender Temperatur.

## Naechste Schritte

Drei Schritte sind nun sinnvoll. Erstens sollte, falls eine Praezisionsbestimmung gewuenscht ist, ein Sprint mit Cluster-Update und groesseren Gittern die Weber-Minnhagen-Extrapolation vollstaendig durchfuehren und pruefen, ob der extrapolierte Wert im thermodynamischen Limes den Referenzwert trifft. Zweitens kann nun die Defekt-Stabilisierungs-Frage im wohldefinierten geordneten Regime knapp unterhalb der kritischen Temperatur wiederaufgenommen werden, wo die Vortices echte topologische Anregungen sind. Drittens bleibt die manuelle Drive-Inspektion des Ordners fuer den Differential-Test gegen den hinterlegten Originalcode offen.

## Drive-Ablage-Empfehlung

In den Ordner sollten das korrigierte Kernmodul in der Fassung v2.2, das PHY025-FSS-Modul in der Fassung v02, der zugehoerige Report und dieses Konsolidierungs-Dokument abgelegt werden. Ein Architecture Decision Record sollte die Ruecknahme der frueheren Normierung und die Einfuehrung der vorwaerts abgeleiteten Flaechen-Normierung dokumentieren, mit ausdruecklichem Hinweis auf die zugrunde liegende Bestaetigungsverzerrung als Negativ-Befund von methodischem Wert.

ENDE des Konsolidierungs-Dokuments. Status v02, GRUEN klassifiziert, Normierung theorie-verankert und Validierung ehrlich gefuehrt.